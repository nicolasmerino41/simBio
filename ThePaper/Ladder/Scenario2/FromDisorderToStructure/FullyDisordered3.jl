using DifferentialEquations, Random, Statistics, LinearAlgebra
using QuadGK, NLsolve, Distributions
using CairoMakie

# --------------------------------------------------
# 1. Single-group analytic solver (Barbier et al. Gaussian approximation)
#    φ = ∫ Dz [1 - Φ((μ*φ + σ*√φ*z - K0)/ζ)]
# --------------------------------------------------
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    std_norm = Normal(0, 1)

    # RHS(φ) = ∫ Dz [1 - Φ((μ*φ + σ√φ z - K0)/ζ)]
    function RHS(φ)
        if φ <= 0
            return 0.0
        elseif φ >= 1
            φc = 1.0
            # still compute integral for φ=1
        else
            φc = φ
        end
        integrand(z) = begin
            arg = (μ*φc + σ*sqrt(φc)*z - K0) / ζ
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        val, _ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
        return val
    end

    f(φ) = RHS(φ) - φ

    # Check endpoints
    f0 = f(1e-8)
    if f0 < 0
        return 0.0
    end
    f1 = f(1 - 1e-8)
    if f1 > 0
        return 1.0
    end
    φ_sol = find_zero(f, (1e-8, 1 - 1e-8), Bisection(); xtol=tol)
    return φ_sol
end

# --------------------------------------------------
# 2. Three-group analytic solver (Gaussian approximation)
#    Groups: basal (b), intermediate (i), top (t).
#    Self-consistent equations:
#      φ_b = ∫ Dz [1 - Φ((μ_bb φ_b + μ_bi φ_i + μ_bt φ_t + sd_b z - K0_b)/ζ_b)]
#      φ_i = ∫ Dz [1 - Φ((μ_ib φ_b + μ_ii φ_i + μ_it φ_t + sd_i z - K0_i)/ζ_i)]
#      φ_t = ∫ Dz [1 - Φ((μ_tb φ_b + μ_ti φ_i + μ_tt φ_t + sd_t z - K0_t)/ζ_t)]
#    Here we use hierarchical no-omnivory: only nonzero blocks are:
#      μ_bi (effect of i on b, negative), μ_ib (effect of b on i, positive),
#      μ_it (effect of t on i, negative), μ_ti (effect of i on t, positive);
#    others zero.
#    Variances similarly σ_bb, σ_bi, etc.
# --------------------------------------------------
function analytical_phi_three_groups(K0_b, ζ_b, K0_i, ζ_i, K0_t, ζ_t,
                                     μ_bb, μ_bi, μ_bt,
                                     μ_ib, μ_ii, μ_it,
                                     μ_tb, μ_ti, μ_tt,
                                     σ_bb, σ_bi, σ_bt,
                                     σ_ib, σ_ii, σ_it,
                                     σ_tb, σ_ti, σ_tt;
                                     tol=1e-6)

    std_norm = Normal(0,1)

    function F!(Fout, φ_vec)
        φ_b, φ_i, φ_t = φ_vec
        # Clamp between 0 and 1 for integral calculation
        φb = clamp(φ_b, 0.0, 1.0)
        φi = clamp(φ_i, 0.0, 1.0)
        φt = clamp(φ_t, 0.0, 1.0)

        # Means for each group
        mean_b = μ_bb*φb + μ_bi*φi + μ_bt*φt
        mean_i = μ_ib*φb + μ_ii*φi + μ_it*φt
        mean_t = μ_tb*φb + μ_ti*φi + μ_tt*φt

        # Variances for each group
        var_b = σ_bb^2*φb + σ_bi^2*φi + σ_bt^2*φt
        var_i = σ_ib^2*φb + σ_ii^2*φi + σ_it^2*φt
        var_t = σ_tb^2*φb + σ_ti^2*φi + σ_tt^2*φt

        sd_b = sqrt(max(var_b, 0.0))
        sd_i = sqrt(max(var_i, 0.0))
        sd_t = sqrt(max(var_t, 0.0))

        # Integrals
        integrand_b(z) = begin
            arg = (mean_b + sd_b*z - K0_b) / ζ_b
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        RHS_b = quadgk(integrand_b, -Inf, Inf; rtol=1e-6)[1]

        integrand_i(z) = begin
            arg = (mean_i + sd_i*z - K0_i) / ζ_i
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        RHS_i = quadgk(integrand_i, -Inf, Inf; rtol=1e-6)[1]

        integrand_t(z) = begin
            arg = (mean_t + sd_t*z - K0_t) / ζ_t
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        RHS_t = quadgk(integrand_t, -Inf, Inf; rtol=1e-6)[1]

        Fout[1] = RHS_b - φ_b
        Fout[2] = RHS_i - φ_i
        Fout[3] = RHS_t - φ_t
    end

    # Initial guess: basal likely high, intermediate small, top smaller
    φ_init = [0.9, 0.3, 0.1]
    sol = nlsolve(F!, φ_init; xtol=tol, ftol=tol)
    φb = clamp(sol.zero[1], 0.0, 1.0)
    φi = clamp(sol.zero[2], 0.0, 1.0)
    φt = clamp(sol.zero[3], 0.0, 1.0)
    return φb, φi, φt
end

# --------------------------------------------------
# 3. Generate hierarchical trophic interaction matrix A
#    Groups: basal (size S_b), intermediate (S_i), top (S_t)
#    Fully hierarchical, no omnivory:
#      - Basal ↔ Intermediate links
#      - Intermediate ↔ Top links
#      - No direct Basal ↔ Top links
#    For each link: effect of resource→consumer positive, consumer→resource negative.
#    Draw Gaussian: mean and variance scaled by donor group size.
# --------------------------------------------------
function gen_alpha_three_trophic(S_b::Int, S_i::Int, S_t::Int,
                                 μ_bi::Float64, μ_ib::Float64,
                                 μ_it::Float64, μ_ti::Float64,
                                 σ_bi::Float64, σ_ib::Float64,
                                 σ_it::Float64, σ_ti::Float64,
                                 rng::AbstractRNG)
    S = S_b + S_i + S_t
    A = zeros(Float64, S, S)
    # Indices:
    # basal: 1:S_b
    # intermediate: S_b+1 : S_b+S_i
    # top: S_b+S_i+1 : S
    # Basal-Intermediate:
    for ib in 1:S_b, ii in (S_b+1):(S_b+S_i)
        # effect basal ib → intermediate ii: positive on intermediate
        # scale mean μ_ib / S_b, var (σ_ib^2)/S_b
        A[ii, ib] = randn(rng)*sqrt((σ_ib^2)/S_b) + (μ_ib / S_b)
        # effect intermediate ii → basal ib: negative on basal
        A[ib, ii] = randn(rng)*sqrt((σ_bi^2)/S_i) - (μ_bi / S_i)
    end
    # Intermediate-Top:
    for ii in (S_b+1):(S_b+S_i), it in (S_b+S_i+1):(S_b+S_i+S_t)
        # effect intermediate ii → top it: positive on top
        A[it, ii] = randn(rng)*sqrt((σ_ti^2)/S_i) + (μ_ti / S_i)
        # effect top it → intermediate ii: negative on intermediate
        A[ii, it] = randn(rng)*sqrt((σ_it^2)/S_t) - (μ_it / S_t)
    end
    # No direct basal-top links; all other entries zero or diagonal zero.
    return A
end

# --------------------------------------------------
# 4. gLV simulation functions
# --------------------------------------------------
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for idx in eachindex(u)
        du[idx] = u[idx] * (K[idx] - u[idx] - Au[idx])
    end
end

function simulate_equilibrium(K, A; t_end=200.0, reltol=1e-6, abstol=1e-6, rng=nothing)
    S = length(K)
    if rng === nothing
        u0 = rand(S) .* K
    else
        u0 = rand(rng, S) .* K
    end
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(), reltol=reltol, abstol=abstol)
    N_end = sol.u[end]
    N_end .= max.(N_end, 0.0)
    return N_end
end

function surviving_fraction(N_end; thresh=1e-3)
    return sum(N_end .> thresh) / length(N_end)
end

# --------------------------------------------------
# 5. Main: compare single-group vs. three-group approximations
# --------------------------------------------------
function compare_three_group_hierarchy()
    rng = MersenneTwister(2025)
    # Community sizes
    S = 300
    frac_b = 1/3
    frac_i = 1/3
    # Basal, intermediate, top sizes
    S_b = Int(round(frac_b * S))
    S_i = Int(round(frac_i * S))
    S_t = S - S_b - S_i

    # Intrinsic K0 and ζ for each group:
    K0_b, ζ_b = 5.0, 1.0      # basal: positive
    K0_i, ζ_i = -1.0, 1.0     # intermediate: negative intrinsic
    K0_t, ζ_t = -1.0, 1.0     # top: negative intrinsic

    # Variance parameters for interactions (we fix):
    σ_ib = 0.5   # basal→intermediate variance
    σ_bi = 0.5   # intermediate→basal variance
    σ_ti = 0.5   # intermediate→top variance
    σ_it = 0.5   # top→intermediate variance

    # Sweep a single mean μ for both links:
    mu_values = range(0.0, stop=5.0, length=21)
    nrep = 20

    # Storage
    μ_list = Float64[]
    φ_total_single_analytic = Float64[]
    φ_total_three_analytic = Float64[]
    φ_total_sim_mean = Float64[]
    φ_total_sim_std = Float64[]

    for μ in mu_values
        # Block means for hierarchical links:
        # basal→intermediate effect: μ_ib = μ (positive for intermediate)
        # intermediate→basal: μ_bi = μ (magnitude, negative on basal)
        μ_ib = μ
        μ_bi = μ
        # intermediate→top: μ_ti = μ (positive for top)
        # top→intermediate: μ_it = μ (magnitude, negative on intermediate)
        μ_ti = μ
        μ_it = μ

        # 5a. Single-group analytic:
        # Need global aggregate moments:
        # We can compute global mean μ_global = S * mean(off-diag A),
        # but in analytic we need μ and σ in form for the approximation:
        # Instead, approximate by average of block means weighted by block sizes.
        # However, simplest: simulate many draws? For clarity, we build a theoretical global mean:
        # For hierarchical blocks, the expected off-diagonal mean:
        #   There are S_b*S_i links basal->intermediate and vice versa,
        #   and S_i*S_t links intermediate->top and vice versa.
        # Effect signs: basal→intermediate (A[ii,ib] positive mean μ_ib/S_b),
        #               intermediate→basal (A[ib,ii] negative mean -μ_bi/S_i),
        # similar for intermediate-top.
        # Compute global mean_offdiag:
        total_links = S*(S-1)  # count off-diagonal entries
        # Sum of means over off-diagonals:
        sum_means = 0.0
        # basal→intermediate entries: count S_i*S_b, mean per entry = μ_ib/S_b
        sum_means += S_i*S_b * (μ_ib/S_b)
        # intermediate→basal: S_b*S_i entries, mean = -μ_bi/S_i
        sum_means += S_b*S_i * (-μ_bi/S_i)
        # intermediate→top: S_t*S_i entries, mean per = μ_ti/S_i
        sum_means += S_t*S_i * (μ_ti/S_i)
        # top→intermediate: S_i*S_t entries, mean = -μ_it/S_t
        sum_means += S_i*S_t * (-μ_it/S_t)
        # All other off-diagonals (within-group and basal-top) have mean zero.
        μ_global = sum_means / total_links * S  # μ_global = S * mean_offdiag

        # For σ_global, assume variances add similarly:
        # sum of variances over off-diagonals:
        sum_vars = 0.0
        # basal→intermediate: S_i*S_b entries, var = (σ_ib^2)/S_b
        sum_vars += S_i*S_b * ((σ_ib^2)/S_b)
        # intermediate→basal: S_b*S_i entries, var = (σ_bi^2)/S_i
        sum_vars += S_b*S_i * ((σ_bi^2)/S_i)
        # intermediate→top: S_t*S_i entries, var = (σ_ti^2)/S_i
        sum_vars += S_t*S_i * ((σ_ti^2)/S_i)
        # top→intermediate: S_i*S_t entries, var = (σ_it^2)/S_t
        sum_vars += S_i*S_t * ((σ_it^2)/S_t)
        # All other off-diagonals have var = 0.
        # Mean variance per off-diagonal:
        var_offdiag = sum_vars / total_links
        σ_global = sqrt(S * var_offdiag)

        # For global K0 and ζ: mix group intrinsic parameters:
        # assume fraction S_b/S have K0_b, S_i/S have K0_i, S_t/S have K0_t.
        K0_global = (S_b*K0_b + S_i*K0_i + S_t*K0_t) / S
        # Global ζ: std of mixture: sqrt( weighted variance + variance of means );
        # approximate by sampling variance: here assume each group has same ζ, but means differ:
        # Var(K) = (S_b*(ζ_b^2 + K0_b^2) + S_i*(ζ_i^2 + K0_i^2) + S_t*(ζ_t^2 + K0_t^2))/S - K0_global^2
        varK = (S_b*(ζ_b^2 + K0_b^2) + S_i*(ζ_i^2 + K0_i^2) + S_t*(ζ_t^2 + K0_t^2))/S - K0_global^2
        ζ_global = sqrt(max(varK, 1e-8))

        # Solve single-group analytic φ_total:
        φ_single = analytical_phi(K0_global, ζ_global, μ_global, σ_global)

        # 5b. Three-group analytic:
        φb_th, φi_th, φt_th = analytical_phi_three_groups(
            K0_b, ζ_b, K0_i, ζ_i, K0_t, ζ_t,
            # block means:
            0.0,             # μ_bb
            μ_bi, 0.0,       # μ_bi, μ_bt=0
            μ_ib, 0.0, μ_it, # μ_ib, μ_ii=0, μ_it
            0.0, μ_ti, 0.0,  # μ_tb=0, μ_ti, μ_tt=0
            # block variances:
            0.0, σ_bi, 0.0,
            σ_ib, 0.0, σ_it,
            0.0, σ_ti, 0.0;
            tol=1e-6
        )
        φ_three = (S_b*φb_th + S_i*φi_th + S_t*φt_th) / S

        # 5c. Simulations:
        totals = Float64[]
        for rep in 1:nrep
            # Generate A
            A = gen_alpha_three_trophic(S_b, S_i, S_t,
                                        μ_bi, μ_ib, μ_it, μ_ti,
                                        σ_bi, σ_ib, σ_it, σ_ti,
                                        rng)
            # Generate K vector
            K = vcat(randn(rng, S_b).*ζ_b .+ K0_b,
                     randn(rng, S_i).*ζ_i .+ K0_i,
                     randn(rng, S_t).*ζ_t .+ K0_t)
            K .= max.(K, 1e-6)
            # Simulate
            N_end = simulate_equilibrium(K, A; rng=rng)
            push!(totals, surviving_fraction(N_end))
        end
        mean_tot = mean(totals)
        std_tot = std(totals)

        # Store
        push!(μ_list, μ)
        push!(φ_total_single_analytic, φ_single)
        push!(φ_total_three_analytic, φ_three)
        push!(φ_total_sim_mean, mean_tot)
        push!(φ_total_sim_std, std_tot)

        @info "μ=$μ | single φ_an=$(round(φ_single,digits=3)), three φ_an=$(round(φ_three,digits=3)), sim=$(round(mean_tot,digits=3))±$(round(std_tot,digits=3))"
    end

    # --------------------------------------------------
    # 6. Plot side-by-side
    # --------------------------------------------------
    fig = Figure(; size = (900, 400))
    # Left: single-group vs simulation
    ax1 = Axis(fig[1, 1];
        xlabel = "Mean interaction magnitude μ",
        ylabel = "Total surviving fraction φ_total",
        title = "Single-group approx vs simulation"
    )
    # Analytical line
    lines!(ax1, μ_list, φ_total_single_analytic; label="Analytic (single-group)", color=:blue)
    # Simulation means with error bars
    for (i, μ) in enumerate(μ_list)
        y = φ_total_sim_mean[i]; err = φ_total_sim_std[i]
        lines!(ax1, [μ, μ], [y - err, y + err]; color=:red, linestyle=:dash)
    end
    scatter!(ax1, μ_list, φ_total_sim_mean; color=:red, marker=:circle, label="Simulated φ_total")
    axislegend(ax1; position = :rt)

    # Right: three-group vs simulation
    ax2 = Axis(fig[1, 2];
        xlabel = "Mean interaction magnitude μ",
        ylabel = "Total surviving fraction φ_total",
        title = "Three-group hierarchical approx vs simulation"
    )
    lines!(ax2, μ_list, φ_total_three_analytic; label="Analytic (three-group)", color=:green)
    for (i, μ) in enumerate(μ_list)
        y = φ_total_sim_mean[i]; err = φ_total_sim_std[i]
        lines!(ax2, [μ, μ], [y - err, y + err]; color=:red, linestyle=:dash)
    end
    scatter!(ax2, μ_list, φ_total_sim_mean; color=:red, marker=:circle, label="Simulated φ_total")
    axislegend(ax2; position = :rt)

    display(fig)
end

# Run the comparison
compare_three_group_hierarchy()
