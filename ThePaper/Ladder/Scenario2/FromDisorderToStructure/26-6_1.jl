using Random, Statistics, LinearAlgebra
using Distributions
using DifferentialEquations
using QuadGK, Roots
using CairoMakie

# -------------------------
# 1. gLV RHS for ODE
#    dN_i/dt = N_i (K_i - N_i - ∑_j A[i,j] N_j)
# -------------------------
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] - Au[i])
    end
end

# -------------------------
# 2. Simulate to equilibrium
# -------------------------
function simulate_equilibrium_phi(A; K0_vec, ζ_vec, group_labels, t_end=200.0, rng=GLOBAL_RNG)
    S = size(A,1)
    # draw K_i by group
    K = similar(group_labels, Float64)
    for i in 1:S
        g = group_labels[i]
        Ki = rand(rng, Normal(K0_vec[g], ζ_vec[g]))
        K[i] = max(Ki, 1e-6)
    end
    # initial N(0)
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
    Nend = sol.u[end]
    Nend .= max.(Nend, 0.0)
    φ = mean(Nend .> 1e-3)
    return φ
end

# -------------------------
# 3. analytical_phi (single‐group prediction)
# -------------------------
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    stdn = Normal(0,1)
    function RHS(φ)
        φc = clamp(φ, 0.0, 1.0)
        integrand(z) = (1 - cdf(stdn, (μ*φc + σ*sqrt(φc)*z - K0)/ζ)) * pdf(stdn, z)
        val, _ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
        return val
    end
    f(φ) = RHS(φ) - φ
    if f(1e-8) < 0
        return 0.0
    elseif f(1-1e-8) > 0
        return 1.0
    else
        return find_zero(f, (1e-8, 1-1e-8), Bisection(); xtol=tol)
    end
end

# -------------------------
# 4. single-group φ_pred from empirical A
# -------------------------
function phi_pred_single(A, K0_vec, ζ_vec, group_labels)
    S = size(A,1); L = length(K0_vec)
    # global K0_glob, ζ_glob
    counts = zeros(Int,L)
    for g in group_labels
        counts[g] += 1
    end
    K0_glob = sum(counts[g]*K0_vec[g] for g in 1:L)/S
    E_K2 = sum(counts[g]*(ζ_vec[g]^2 + K0_vec[g]^2) for g in 1:L)/S
    ζ_glob = sqrt(max(E_K2 - K0_glob^2, 1e-8))
    # A moments
    off = [A[i,j] for i in 1:S for j in 1:S if i!=j]
    μ_glob = S * mean(off)
    σ_glob = sqrt(S * var(off))
    return analytical_phi(K0_glob, ζ_glob, μ_glob, σ_glob)
end

# -------------------------
# 5. power-law network generation
# -------------------------
function gen_PL_edges(S::Int, γ::Float64, kmin::Int, rng::AbstractRNG)
    kmax = S-1
    c1 = kmin^(1-γ)
    c2 = kmax^(1-γ)
    denom = c2 - c1   # note: negative for γ>1, but that’s fine

    # 1. sample degrees
    deg = zeros(Int, S)
    for i in 1:S
        u = rand(rng)
        # inverse‐transform for x^(1-γ) ∈ [c1, c2]
        kf = (u*denom + c1)^(1/(1-γ))
        # clamp to [kmin,kmax]
        deg[i] = clamp(floor(Int, kf), kmin, kmax)
    end

    # 2. make sum(deg) even
    isodd(sum(deg)) && (deg[rand(rng,1:S)] += 1)

    # 3. configuration‐model stubs
    stubs = Int[]
    for i in 1:S
        append!(stubs, fill(i, deg[i]))
    end
    shuffle!(stubs)

    # 4. pair stubs into undirected edges
    edges = Tuple{Int,Int}[]
    while length(stubs) >= 2
        i = pop!(stubs)
        j = pop!(stubs)
        if i != j
            push!(edges, (i, j))
        end
        # self‐loops are discarded
    end

    # 5. randomize direction
    directed = Tuple{Int,Int}[]
    for (i,j) in edges
        rand(rng) < 0.5 ? push!(directed, (i,j)) : push!(directed, (j,i))
    end

    return directed
end

# -------------------------
# 6. build trophic A from edges
# -------------------------
function build_A(S, edges, μ0, σ0, rng)
    A = zeros(S,S)
    for (i,j) in edges
        # sample positive magnitude
        mean_pos = μ0 / S
        sd_pos = sqrt(σ0^2 / S)
        mag = rand(rng, truncated(Normal(mean_pos, sd_pos), 0, Inf))
        A[j,i] = mag
        A[i,j] = -mag
    end
    return A
end

# -------------------------
# 7. main sweep
# -------------------------
function sweep_powerlaw_ODE()
    rng = MersenneTwister(42)
    S = 60
    # groups only for K-draw; no ordering
    groups = vcat(fill(1,20),fill(2,20),fill(3,20))
    K0_vec = [5.0, 3.0, 1.0]
    ζ_vec  = [1.0, 1.0, 1.0]
    σ0 = 0.5
    μ0 = 1.0
    gammas = [2.1,2.5,3.0,3.5,4.0]
    nrep = 50

    cvs = Float64[]; errs = Float64[]; gs = Float64[]
    for γ in gammas
        for rep in 1:nrep
            edges = gen_PL_edges(S, γ, 1, rng)
            deg = zeros(Int,S)
            for (i,j) in edges
                deg[i] += 1; deg[j] += 1
            end
            cv = std(deg)/mean(deg)
            A = build_A(S, edges, μ0, σ0, rng)
            φs = simulate_equilibrium_phi(A; K0_vec=K0_vec, ζ_vec=ζ_vec, group_labels=groups, rng=rng)
            φp = phi_pred_single(A, K0_vec, ζ_vec, groups)
            push!(cvs, cv); push!(errs, abs(φp-φs)); push!(gs, γ)
        end
    end

    # plot error vs CV
    fig = Figure(; size=(500,350))
    ax = Axis(fig[1,1];
        xlabel="Degree CV",
        ylabel="|φ_pred − φ_sim|",
        title="Power-law heterogeneity (ODE) → error vs CV"
    )
    scatter!(ax, cvs, errs; color=gs, colormap=:viridis, marker=:circle)
    display(fig)
end

# run
sweep_powerlaw_ODE()
