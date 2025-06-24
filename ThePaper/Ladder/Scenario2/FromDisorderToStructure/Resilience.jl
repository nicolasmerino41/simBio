using Random, Statistics, LinearAlgebra
using DifferentialEquations
using QuadGK, Roots, Distributions
using CairoMakie

# --------------------------------------------------
# 1. Analytical single-group φ and mean abundance
# --------------------------------------------------
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    std_norm = Normal(0, 1)
    function RHS(φ)
        if φ <= 0
            return 0.0
        elseif φ >= 1
            φc = 1.0
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
    f0 = f(1e-8)
    if f0 < 0
        return 0.0
    end
    f1 = f(1-1e-8)
    if f1 > 0
        return 1.0
    end
    φ_sol = find_zero(f, (1e-8, 1-1e-8), Bisection(); xtol=tol)
    return φ_sol
end

function approximate_mean_abundance(K0, μ, φ)
    val = K0 - μ*φ
    return max(val, 0.0)
end

# --------------------------------------------------
# 2. Two-group analytic φ_b, φ_c
# --------------------------------------------------
function analytical_phi_two_groups(K0_b, ζ_b, K0_c, ζ_c,
                                   μ_bb, μ_bc, μ_cb, μ_cc,
                                   σ_bb, σ_bc, σ_cb, σ_cc;
                                   tol=1e-6)
    std_norm = Normal(0,1)
    function F!(Fout, φ_vec)
        φ_b, φ_c = φ_vec
        φb = clamp(φ_b, 0.0, 1.0)
        φc = clamp(φ_c, 0.0, 1.0)
        mean_b = μ_bb*φb + μ_bc*φc
        mean_c = μ_cb*φb + μ_cc*φc
        var_b = σ_bb^2*φb + σ_bc^2*φc
        var_c = σ_cb^2*φb + σ_cc^2*φc
        sd_b = sqrt(max(var_b,0.0))
        sd_c = sqrt(max(var_c,0.0))
        integrand_b(z) = begin
            arg = (mean_b + sd_b*z - K0_b)/ζ_b
            (1 - cdf(std_norm, arg))*pdf(std_norm, z)
        end
        RHS_b = quadgk(integrand_b, -Inf, Inf; rtol=1e-6)[1]
        integrand_c(z) = begin
            arg = (mean_c + sd_c*z - K0_c)/ζ_c
            (1 - cdf(std_norm, arg))*pdf(std_norm, z)
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

function approximate_mean_abundances_two_groups(K0_b, μ_bc, φ_c, K0_c, μ_cb, φ_b)
    Nb = max(K0_b - μ_bc*φ_c, 0.0)
    Nc = max(K0_c + μ_cb*φ_b, 0.0)
    return Nb, Nc
end

# --------------------------------------------------
# 3. Build bipartite trophic interaction matrix
# --------------------------------------------------
function gen_alpha_trophic(S_b::Int, S_c::Int,
                           μ_bc::Float64, μ_cb::Float64,
                           σ_bc::Float64, σ_cb::Float64,
                           rng::AbstractRNG)
    S = S_b + S_c
    A = zeros(Float64, S, S)
    for i in 1:S_b, j in (S_b+1):(S_b+S_c)
        A[j, i] = randn(rng)*sqrt((σ_cb^2)/S_b) + (μ_cb / S_b)
        A[i, j] = randn(rng)*sqrt((σ_bc^2)/S_c) - (μ_bc / S_c)
    end
    return A
end

# --------------------------------------------------
# 4. gLV simulation + empirical resilience
# --------------------------------------------------
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for idx in eachindex(u)
        du[idx] = u[idx] * (K[idx] - u[idx] - Au[idx])
    end
end

"""
    gLV_survive_and_resilience(A; K0_b, ζ_b, K0_c, ζ_c, S_b, S_c, rng)

Simulate gLV, compute survivors, build Jacobian on survivors,
compute largest real part of eigenvalues of J, then resilience = -largest_real_part.
Returns (φ_total, resilience).
"""
function gLV_survive_and_resilience(A; K0_b, ζ_b, K0_c, ζ_c, S_b::Int, S_c::Int, rng=GLOBAL_RNG)
    S = S_b + S_c
    K = vcat(randn(rng, S_b).*ζ_b .+ K0_b,
             randn(rng, S_c).*ζ_c .+ K0_c)
    K .= max.(K, 1e-6)
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, 200.0), (K, A))
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    N_end = sol.u[end]
    N_end .= max.(N_end, 0.0)
    survivors = findall(x-> x>1e-3, N_end)
    φ_total = length(survivors)/S
    if length(survivors) <= 1
        # trivial: if ≤1 survivor, resilience = -(Jacobian eigenvalue) 
        # but we treat as large positive resilience or Inf; here return large positive
        return φ_total, Inf
    end
    idx = survivors
    n = length(idx)
    Jsub = zeros(Float64, n, n)
    for (ii, i) in enumerate(idx)
        if N_end[i] > 1e-8
            Jsub[ii, ii] = -N_end[i]
            for (jj, j) in enumerate(idx)
                if i != j
                    Jsub[ii, jj] = -N_end[i] * A[i, j]
                end
            end
        end
    end
    # full Jacobian Jsub; compute eigenvalues
    vals = eigen(Jsub).values
    λmax_real = maximum(real.(vals))
    # resilience = -λmax_real
    resilience = -λmax_real
    return φ_total, resilience
end

# --------------------------------------------------
# 5. Analytic resilience predictions
# --------------------------------------------------
"""
Single-group predicted resilience from matrix A:
- Compute global K0_glob, ζ_glob, global μ_glob, σ_glob from A and K0_b, K0_c.
- Solve φ_pred = analytical_phi(...)
- Approx mean abundance Nbar = K0_glob - μ_glob*φ_pred
- Predicted resilience r_pred = Nbar - Nbar * σ_glob * sqrt(φ_pred)
"""
function predicted_resilience_single_from_matrix(A, K0_b, ζ_b, K0_c, ζ_c, S_b, S_c)
    S = S_b + S_c
    # global K0, ζ
    K0_glob = (S_b*K0_b + S_c*K0_c)/S
    varK = (S_b*(ζ_b^2 + K0_b^2) + S_c*(ζ_c^2 + K0_c^2))/S - K0_glob^2
    ζ_glob = sqrt(max(varK, 1e-8))
    # global moments from A
    off = Float64[]
    for i in 1:S, j in 1:S
        if i != j
            push!(off, A[i,j])
        end
    end
    meanA = mean(off)
    varA = var(off)
    μ_glob = S * meanA
    σ_glob = sqrt(S * varA)
    φ_pred = analytical_phi(K0_glob, ζ_glob, μ_glob, σ_glob)
    Nbar = approximate_mean_abundance(K0_glob, μ_glob, φ_pred)
    if φ_pred < 1e-8 || Nbar < 1e-8
        return φ_pred, Inf
    end
    # predicted resilience
    rpred = Nbar - Nbar * σ_glob * sqrt(φ_pred)
    return φ_pred, rpred
end

"""
Two-group predicted resilience:
- Solve φ_b, φ_c analytically.
- Approx mean abundances Nbar_b, Nbar_c.
- n_b = S_b * φ_b, n_c = S_c * φ_c.
- Compute v_bc = 0.25*(Nbar_b^2*(σ_bc^2/S_c) + Nbar_c^2*(σ_cb^2/S_b)).
- Approx spectral radius of -N*A ≈ (n_b * n_c)^(1/4)*sqrt(v_bc).
- avg abundance Nbar_avg = (n_b*Nbar_b + n_c*Nbar_c)/(n_b+n_c).
- Predicted resilience rpred = Nbar_avg - (n_b*n_c)^(1/4)*sqrt(v_bc).
"""
function predicted_resilience_two_group(K0_b, ζ_b, K0_c, ζ_c,
                                         μ_bc, μ_cb, σ_bc, σ_cb,
                                         S_b::Int, S_c::Int)
    φ_b, φ_c = analytical_phi_two_groups(K0_b, ζ_b, K0_c, ζ_c,
                                         0.0, μ_bc, μ_cb, 0.0,
                                         0.0, σ_bc, σ_cb, 0.0)
    Nbar_b, Nbar_c = approximate_mean_abundances_two_groups(K0_b, μ_bc, φ_c,
                                                            K0_c, μ_cb, φ_b)
    n_b = S_b * φ_b
    n_c = S_c * φ_c
    if n_b < 1e-8 || n_c < 1e-8
        # one or both extinct: resilience dominated by surviving group: = -(-Nbar) = Nbar
        if n_b < 1e-8 && n_c < 1e-8
            return (S_b*φ_b + S_c*φ_c)/(S_b+S_c), Inf
        elseif n_b < 1e-8
            return (S_c*φ_c)/(S_b+S_c), Nbar_c
        else
            return (S_b*φ_b)/(S_b+S_c), Nbar_b
        end
    end
    varA_bc = σ_bc^2 / S_c
    varA_cb = σ_cb^2 / S_b
    v_bc = 0.25*(Nbar_b^2 * varA_bc + Nbar_c^2 * varA_cb)
    # approximate spectral radius ~ (n_b * n_c)^(1/4) * sqrt(v_bc)
    spec_rad = (n_b * n_c)^(1/4) * sqrt(v_bc)
    Nbar_avg = (n_b * Nbar_b + n_c * Nbar_c) / (n_b + n_c)
    rpred = Nbar_avg - spec_rad
    φ_total = (S_b*φ_b + S_c*φ_c)/(S_b+S_c)
    return φ_total, rpred
end

# --------------------------------------------------
# 6. Experiment: simulate & compare resilience predictions
# --------------------------------------------------
function experiment_resilience_two_group()
    rng = MersenneTwister(123)
    S = 200
    frac_b = 0.5
    S_b = Int(round(frac_b*S))
    S_c = S - S_b

    K0_b, ζ_b = 5.0, 1.0
    K0_c, ζ_c = -1.0, 1.0
    σ_bc = 0.5
    σ_cb = 0.5

    mu_values = range(0.0, stop=5.0, length=21)
    nrep = 10

    μ_list = Float64[]
    sim_vals = Float64[]
    sim_std = Float64[]
    pred1_vals = Float64[]
    pred2_vals = Float64[]

    for μ in mu_values
        sim_res = Float64[]
        pred1_res = Float64[]
        pred2_res = Float64[]
        for rep in 1:nrep
            A = gen_alpha_trophic(S_b, S_c, μ, μ, σ_bc, σ_cb, rng)
            φ_sim, r_sim = gLV_survive_and_resilience(A; K0_b=K0_b, ζ_b=ζ_b,
                                                      K0_c=K0_c, ζ_c=ζ_c,
                                                      S_b=S_b, S_c=S_c, rng=rng)
            push!(sim_res, r_sim)
            _, r1 = predicted_resilience_single_from_matrix(A, K0_b, ζ_b, K0_c, ζ_c, S_b, S_c)
            push!(pred1_res, r1)
            _, r2 = predicted_resilience_two_group(K0_b, ζ_b, K0_c, ζ_c,
                                                   μ, μ, σ_bc, σ_cb,
                                                   S_b, S_c)
            push!(pred2_res, r2)
        end
        push!(μ_list, μ)
        push!(sim_vals, mean(sim_res))
        push!(sim_std, std(sim_res))
        push!(pred1_vals, mean(pred1_res))
        push!(pred2_vals, mean(pred2_res))
    end

    # Plot: left single-group overlay, right two-group overlay
    fig = Figure(; size=(900,400))

    ax1 = Axis(fig[1,1];
        xlabel="Mean interaction μ",
        ylabel="Resilience (-Re(λ_max))",
        title="Simulated resilience with single-group prediction"
    )
    for i in eachindex(μ_list)
        μ = μ_list[i]
        lines!(ax1, [μ, μ], [sim_vals[i]-sim_std[i], sim_vals[i]+sim_std[i]]; color=:gray, linestyle=:dash)
    end
    scatter!(ax1, μ_list, sim_vals; color=:black, marker=:circle, label="Simulated resilience")
    lines!(ax1, μ_list, pred1_vals; color=:blue, linewidth=2, label="Single-group predicted")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1,2];
        xlabel="Mean interaction μ",
        ylabel="Resilience (-Re(λ_max))",
        title="Simulated resilience with two-group prediction"
    )
    for i in eachindex(μ_list)
        μ = μ_list[i]
        lines!(ax2, [μ, μ], [sim_vals[i]-sim_std[i], sim_vals[i]+sim_std[i]]; color=:gray, linestyle=:dash)
    end
    scatter!(ax2, μ_list, sim_vals; color=:black, marker=:circle, label="Simulated resilience")
    lines!(ax2, μ_list, pred2_vals; color=:green, linewidth=2, label="Two-group predicted")
    axislegend(ax2; position=:rt)

    display(fig)
end

# Run Script 2 adjusted
experiment_resilience_two_group()
