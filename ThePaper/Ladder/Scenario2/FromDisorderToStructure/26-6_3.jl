using Random, Statistics, LinearAlgebra
using Distributions
using DifferentialEquations
using QuadGK, Roots
using CairoMakie

# 1. gLV RHS
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] - Au[i])
    end
end

# 2. simulate φ_sim via ODE
function φ_sim(A; K0, ζ, t_end=200.0, rng=GLOBAL_RNG)
    S = size(A,1)
    K = max.(randn(rng, S).*ζ .+ K0, 1e-6)
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
    Nend = max.(sol.u[end], 0.0)
    return mean(Nend .> 1e-3)
end

# 3. single‐group analytical φ_pred
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    stdn = Normal(0,1)
    function RHS(φ)
        φc = clamp(φ,0.0,1.0)
        integrand(z) = (1 - cdf(stdn, (μ*φc + σ*sqrt(φc)*z - K0)/ζ)) * pdf(stdn, z)
        val,_ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
        return val
    end
    f(φ) = RHS(φ) - φ
    if f(1e-8) < 0
        return 0.0
    elseif f(1-1e-8) > 0
        return 1.0
    else
        return find_zero(f, (1e-8,1-1e-8), Bisection(); xtol=tol)
    end
end

function φ_pred(A; K0, ζ)
    S = size(A,1)
    off = [A[i,j] for i in 1:S, j in 1:S if i!=j]
    μ_glob = S*mean(off)
    σ_glob = sqrt(S*var(off))
    return analytical_phi(K0, ζ, μ_glob, σ_glob)
end

# 4. partial‐ordering trophic edges
function gen_edges(S, C, α, rng)
    edges = Tuple{Int,Int}[]
    for i in 1:S-1, j in i+1:S
        if rand(rng) < C
            if rand(rng) < α
                push!(edges,(i,j))      # cascade
            else
                rand(rng)<0.5 ? push!(edges,(i,j)) : push!(edges,(j,i))
            end
        end
    end
    return edges
end

# 5. build A
function build_A(S, edges, μ0, σ0, rng)
    A = zeros(S,S)
    for (i,j) in edges
        m = rand(rng, truncated(Normal(μ0/S, sqrt(σ0^2/S)), 0, Inf))
        A[j,i] = m
        A[i,j] = -m
    end
    return A
end

# 6. sweep and collect errors
function collect_errors()
    rng = MersenneTwister(2025)
    S = 60
    C = 0.1
    K0, ζ = 5.0, 1.0
    σ0 = 0.5

    alphas = [0.0, 0.25, 0.5, 0.75, 1.0]
    mu0s   = [1.0, 5.0, 10.0]
    nrep   = 5

    # preallocate
    err_mat = zeros(length(alphas), length(mu0s))

    for (iα, α) in enumerate(alphas)
        for (iμ, μ0) in enumerate(mu0s)
            errs = Float64[]
            for rep in 1:nrep
                edges = gen_edges(S, C, α, rng)
                A = build_A(S, edges, μ0, σ0, rng)
                φs = φ_sim(A; K0=K0, ζ=ζ, rng=rng)
                φp = φ_pred(A; K0=K0, ζ=ζ)
                push!(errs, abs(φp - φs))
            end
            err_mat[iα, iμ] = mean(errs)
        end
    end

    return alphas, mu0s, err_mat
end

# 7. plot two panels
function plot_errors()
    alphas, mu0s, errs = collect_errors()

    fig = Figure(; size=(800,350))

    # left: error vs μ for each α
    ax1 = Axis(fig[1,1];
        xlabel="μ₀",
        ylabel="Mean |φₚᵣₑd−φₛᵢₘ|",
        title="Error vs interaction strength for each α"
    )
    for (iα, α) in enumerate(alphas)
        lines!(ax1, mu0s, errs[iα, :],label="α=$(α)")
    end
    axislegend(ax1; position=:lt)

    # right: error vs α for each μ₀
    ax2 = Axis(fig[1,2];
        xlabel="α",
        ylabel="Mean |φₚᵣₑd−φₛᵢₘ|",
        title="Error vs ordering for each μ₀"
    )
    for (iμ, μ0) in enumerate(mu0s)
        lines!(ax2, alphas, errs[:, iμ], label="μ₀=$(μ0)")
    end
    axislegend(ax2; position=:lt)

    display(fig)
end

# Run it
plot_errors()

function plot_pred_vs_sim()
    rng = MersenneTwister(2025)
    # Reuse your parameter grid
    alphas = [0.0, 0.25, 0.5, 0.75, 1.0]
    mu0s   = [1.0, 5.0, 10.0]
    nrep   = 5

    # Storage matrices
    φ_sim_mat  = zeros(length(alphas), length(mu0s))
    φ_pred_mat = zeros(length(alphas), length(mu0s))

    # Re‐compute means of φ_sim and φ_pred over replicates
    for (iα, α) in enumerate(alphas), (iμ, μ0) in enumerate(mu0s)
        sims = Float64[]
        preds = Float64[]
        for rep in 1:nrep
            edges = gen_trophic_edges(60, 0.1, α, rng)
            A     = build_A(60, edges, μ0, 0.5, rng)
            push!(sims, simulate_equilibrium_phi(A; K0=5.0, ζ=1.0, rng))
            push!(preds, φ_pred_single(A; K0=5.0, ζ=1.0))
        end
        φ_sim_mat[iα, iμ]  = mean(sims)
        φ_pred_mat[iα, iμ] = mean(preds)
    end

    # Now plot
    fig = Figure(; size=(800,350))

    # Left: φ vs α for each mu0
    ax1 = Axis(fig[1,1],
        xlabel="Ordering α",
        ylabel="Surviving fraction φ",
        title="Predicted (lines) vs simulated (dots) φ vs α"
    )
    for (iμ, μ0) in enumerate(mu0s)
        lines!(ax1, alphas, φ_pred_mat[:,iμ], label="pred μ₀=$(μ0)", linewidth=2)
        scatter!(ax1, alphas, φ_sim_mat[:,iμ], label="sim μ₀=$(μ0)", marker=:circle)
    end
    axislegend(ax1; position=:lb)

    # Right: φ vs μ₀ for each α
    ax2 = Axis(fig[1,2],
        xlabel="Interaction strength μ₀",
        ylabel="Surviving fraction φ",
        title="Predicted (lines) vs simulated (dots) φ vs μ₀"
    )
    for (iα, α) in enumerate(alphas)
        lines!(ax2, mu0s, φ_pred_mat[iα,:], label="pred α=$(α)", linewidth=2)
        scatter!(ax2, mu0s, φ_sim_mat[iα,:], label="sim α=$(α)", marker=:diamond)
    end
    # axislegend(ax2; position=:lb)

    display(fig)
end

# Call it:
plot_pred_vs_sim()
