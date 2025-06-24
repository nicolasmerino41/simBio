using DifferentialEquations, Random, Statistics, LinearAlgebra
using QuadGK, NLsolve, Distributions
using CairoMakie

# --------------------------------------------------
# Single-group analytical solver: φ = ∫ Dz [1 - Φ((μ*φ + σ*√φ*z - K0)/ζ)]
# --------------------------------------------------
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    std_norm = Normal(0, 1)
    function RHS(φ)
        if φ <= 0
            return 0.0
        elseif φ >= 1
            φ = 1.0
        end
        integrand(z) = begin
            arg = (μ*φ + σ*sqrt(φ)*z - K0) / ζ
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        val, _ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
        return val
    end
    f(φ) = RHS(φ) - φ
    # Check trivial cases
    if f(1.0 - 1e-8) > 0
        return 1.0
    elseif f(1e-8) < 0
        return 0.0
    end
    φ_sol = find_zero(f, (1e-8, 1.0 - 1e-8), Bisection(); xtol=tol)
    return φ_sol
end

# --------------------------------------------------
# Two-group analytical solver (as before)
# --------------------------------------------------
function analytical_phi_two_groups(K0_b, ζ_b, K0_c, ζ_c,
                                   μ_bb, μ_bc, μ_cb, μ_cc,
                                   σ_bb, σ_bc, σ_cb, σ_cc;
                                   tol=1e-6)

    std_norm = Normal(0,1)
    function F!(Fout, φ_vec)
        φ_b, φ_c = φ_vec
        φ_bc = clamp(φ_b, 0.0, 1.0)
        φ_cc = clamp(φ_c, 0.0, 1.0)

        mean_b = μ_bb * φ_bc + μ_bc * φ_cc
        mean_c = μ_cb * φ_bc + μ_cc * φ_cc
        var_b = σ_bb^2 * φ_bc + σ_bc^2 * φ_cc
        var_c = σ_cb^2 * φ_bc + σ_cc^2 * φ_cc
        sd_b = sqrt(max(var_b, 0.0))
        sd_c = sqrt(max(var_c, 0.0))

        integrand_b(z) = begin
            arg = (mean_b + sd_b*z - K0_b) / ζ_b
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        RHS_b = quadgk(integrand_b, -Inf, Inf; rtol=1e-6)[1]

        integrand_c(z) = begin
            arg = (mean_c + sd_c*z - K0_c) / ζ_c
            (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
        end
        RHS_c = quadgk(integrand_c, -Inf, Inf; rtol=1e-6)[1]

        Fout[1] = RHS_b - φ_b
        Fout[2] = RHS_c - φ_c
    end

    φ_init = [0.5, 0.1]
    sol = nlsolve(F!, φ_init; xtol=tol, ftol=tol)
    φb = clamp(sol.zero[1], 0.0, 1.0)
    φc = clamp(sol.zero[2], 0.0, 1.0)
    return φb, φc
end

# --------------------------------------------------
# Generate trophic interaction matrix (Gaussian blocks)
# --------------------------------------------------
function gen_alpha_trophic(S_b::Int, S_c::Int,
                           μ_bc::Float64, μ_cb::Float64,
                           σ_bc::Float64, σ_cb::Float64,
                           rng::AbstractRNG)
    S = S_b + S_c
    A = zeros(Float64, S, S)
    for i in 1:S_b, j in (S_b+1):(S_b+S_c)
        # basal i → consumer j: positive effect on consumer
        A[j, i] = randn(rng) * sqrt((σ_cb^2)/S_b) + (μ_cb / S_b)
        # consumer j → basal i: negative effect on basal
        A[i, j] = randn(rng) * sqrt((σ_bc^2)/S_c) - (μ_bc / S_c)
    end
    return A
end

# --------------------------------------------------
# gLV simulation functions
# --------------------------------------------------
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] - Au[i])
    end
end

function simulate_equilibrium(K, A; t_end=200.0, reltol=1e-6, abstol=1e-6, rng=nothing)
    S = length(K)
    u0 = rng===nothing ? rand(S).*K : rand(rng, S).*K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(), reltol=reltol, abstol=abstol)
    N_end = sol.u[end]
    N_end .= max.(N_end, 0.0)
    return N_end
end

# --------------------------------------------------
# Main: compare single-group vs two-group approximations, side by side
# --------------------------------------------------
function compare_single_two_group()
    rng = MersenneTwister(1234)
    S = 200
    frac_b = 0.5
    S_b = Int(round(frac_b * S))
    S_c = S - S_b

    # Intrinsic K0: basal positive, consumer negative
    K0_b, ζ_b = 5.0, 1.0
    K0_c, ζ_c = -1.0, 1.0

    # Block variance parameters
    σ_bc = 0.5
    σ_cb = 0.5

    # Sweep mean interaction magnitude μ for both directions
    mu_values = range(0.0, stop=5.0, length=21)
    nrep = 20

    μ_list = Float64[]
    φ_sg_analytic = Float64[]      # single-group analytic
    φ_tg_analytic = Float64[]      # two-group analytic total
    φ_sim_mean = Float64[]
    φ_sim_std = Float64[]

    # Precompute global K0 and ζ for single-group:
    # The global distribution of K_i is a mixture: basal ~ N(K0_b, ζ_b^2), consumer ~ N(K0_c, ζ_c^2)
    # Global mean:
    K0_global = frac_b*K0_b + (1-frac_b)*K0_c
    # Global variance = E[K^2] - K0_global^2
    # E[K^2] = frac_b*(ζ_b^2 + K0_b^2) + (1-frac_b)*(ζ_c^2 + K0_c^2)
    E_K2 = frac_b*(ζ_b^2 + K0_b^2) + (1-frac_b)*(ζ_c^2 + K0_c^2)
    ζ_global = sqrt(max(E_K2 - K0_global^2, 0.0))

    for μ in mu_values
        # Block means for two-group
        μ_bc = μ
        μ_cb = μ

        # --- Two-group analytic ---
        φb_th, φc_th = analytical_phi_two_groups(
            K0_b, ζ_b, K0_c, ζ_c,
            0.0, μ_bc, μ_cb, 0.0,
            0.0, σ_bc, σ_cb, 0.0;
            tol=1e-6
        )
        φ_tg = (S_b*φb_th + S_c*φc_th) / S

        # --- Single-group analytic ---
        # Compute global μ_global and σ_global from block parameters:
        # Off-diagonal pairs consist only of cross-group links; within-group zeros.
        # Total off-diagonal count ~ S*(S-1)
        # Sum of means over i≠j: sum_{i∈basal, j∈consumer} μ_cb/S_b  + sum_{i∈consumer, j∈basal} -μ_bc/S_c
        # Numerator_mean = S_b*S_c*(μ_cb/S_b) + S_b*S_c*(-μ_bc/S_c) = S_c*μ_cb - S_b*μ_bc
        # μ_global = Numerator_mean / (S*(S-1))
        numerator = S_c*μ_cb - S_b*μ_bc
        μ_global = numerator / (S*(S-1))

        # For variance: E[A^2] over i≠j:
        # sum_{basal→consumer} [(σ_cb^2/S_b + (μ_cb/S_b)^2)] * (S_b*S_c)
        # + sum_{consumer→basal} [(σ_bc^2/S_c + (μ_bc/S_c)^2)] * (S_b*S_c)
        # divided by (S*(S-1)), then σ_global = sqrt(S * Var[A]) in the gLV normalization?
        # But in single-group analytic we need μ_global and σ_global such that α_ij ~ N(μ_global, σ_global^2/S).
        # Actually we compute Var_off = E[A^2] - (E[A])^2:
        E_A2 = (S_b*S_c)*((σ_cb^2)/S_b + (μ_cb/S_b)^2 + (σ_bc^2)/S_c + (μ_bc/S_c)^2) / (S*(S-1))
        Var_A = E_A2 - μ_global^2
        # Then σ_global such that Var(A_ij) ≈ σ_global^2 / S  => σ_global = sqrt(S * Var_A)
        σ_global = sqrt(max(S * Var_A, 0.0))

        # Analytical single-group φ
        φ_sg = analytical_phi(K0_global, ζ_global, μ_global*S, σ_global; tol=1e-6)
        # Note: analytical_phi expects arguments K0, ζ, μ, σ where μ is the total mean load (=S * mean per-entry)
        # We computed μ_global as mean per-entry (~1/(S) scale). So pass μ_global * S.

        # --- Simulations ---
        totals = Float64[]
        for rep in 1:nrep
            A = gen_alpha_trophic(S_b, S_c, μ_bc, μ_cb, σ_bc, σ_cb, rng)
            K = vcat(randn(rng, S_b).*ζ_b .+ K0_b,
                     randn(rng, S_c).*ζ_c .+ K0_c)
            K .= max.(K, 1e-6)
            N_end = simulate_equilibrium(K, A; rng=rng)
            push!(totals, sum(N_end .> 1e-3) / S)
        end
        mean_sim = mean(totals)
        std_sim = std(totals)

        push!(μ_list, μ)
        push!(φ_sg_analytic, φ_sg)
        push!(φ_tg_analytic, φ_tg)
        push!(φ_sim_mean, mean_sim)
        push!(φ_sim_std, std_sim)

        @info "μ=$μ | SG_analytic=$(round(φ_sg,digits=3)), TG_analytic=$(round(φ_tg,digits=3)), sim=$(round(mean_sim,digits=3))"
    end

    # --------------------------------------------------
    # Plot side by side
    # --------------------------------------------------
    fig = Figure(; size = (1000, 400))

    # Single-group plot on left
    ax1 = Axis(fig[1, 1];
        xlabel = "Mean interaction μ",
        ylabel = "Overall surviving fraction φ",
        title = "Single-group approximation"
    )
    lines!(ax1, μ_list, φ_sg_analytic; color=:blue, label="Analytic (single-group)")
    for (i, μ) in enumerate(μ_list)
        y = φ_sim_mean[i]; err = φ_sim_std[i]
        lines!(ax1, [μ, μ], [y-err, y+err]; color=:red, linestyle=:dash)
    end
    scatter!(ax1, μ_list, φ_sim_mean; color=:red, marker=:circle, label="Simulation")
    axislegend(ax1; position=:rt)

    # Two-group plot on right
    ax2 = Axis(fig[1, 2];
        xlabel = "Mean interaction μ",
        ylabel = "Overall surviving fraction φ",
        title = "Two-group approximation"
    )
    lines!(ax2, μ_list, φ_tg_analytic; color=:blue, label="Analytic (two-group)")
    for (i, μ) in enumerate(μ_list)
        y = φ_sim_mean[i]; err = φ_sim_std[i]
        lines!(ax2, [μ, μ], [y-err, y+err]; color=:red, linestyle=:dash)
    end
    scatter!(ax2, μ_list, φ_sim_mean; color=:red, marker=:circle, label="Simulation")
    axislegend(ax2; position=:rt)

    display(fig)
end

# Run the comparison
compare_single_two_group()
