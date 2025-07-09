using Random, Statistics, LinearAlgebra
using Distributions, DifferentialEquations
using QuadGK, Roots
using CairoMakie

# 1. gLV RHS (same as before)
function gLV_rhs!(du,u,p,t)
    K,A = p
    Au = A*u
    @inbounds for i in eachindex(u)
        du[i] = u[i]*(K[i] - u[i] - Au[i])
    end
end

# 2. equilibrium φ via ODE
function φ_sim(A; K0,ζ,S,rng)
    # draw K_i iid ~ N(K0,ζ)
    K = max.(randn(rng,S).*ζ .+ K0, 1e-6)
    u0 = rand(rng,S).*K
    u0 .= max.(u0,1e-6)
    prob = ODEProblem(gLV_rhs!,u0,(0.0,200.0),(K,A))
    sol = solve(prob,Tsit5();reltol=1e-6,abstol=1e-6)
    Nend = max.(sol.u[end],0.0)
    return mean(Nend .> 1e-3)
end

# 3. single-group φ_pred
function φ_pred(A; K0,ζ)
    S = size(A,1)
    # global K0,ζ unchanged
    # A moments
    off = [A[i,j] for i in 1:S, j in 1:S if i!=j]
    μ_glob = S*mean(off)
    σ_glob = sqrt(S*var(off))
    return analytical_phi(K0,ζ,μ_glob,σ_glob)
end

# 4. network ensembles with partial ordering α
function gen_trophic_edges_cascade(S, C, α, rng)
    edges = Tuple{Int,Int}[]
    for i in 1:S-1, j in i+1:S
        if rand(rng) < C
            if rand(rng) < α
                # cascade: direct i→j (i lower rank → j higher)
                push!(edges,(i,j))
            else
                # random direction
                rand(rng)<0.5 ? push!(edges,(i,j)) : push!(edges,(j,i))
            end
        end
    end
    return edges
end

# 5. build A from edges
function build_A(S,edges,μ0,σ0,rng)
    A = zeros(S,S)
    for (i,j) in edges
        mean_pos = μ0/S
        sd_pos   = sqrt(σ0^2/S)
        mag = rand(rng, truncated(Normal(mean_pos,sd_pos),0,Inf))
        A[j,i]=mag; A[i,j]=-mag
    end
    return A
end

# 6. sweep α
function sweep_alpha()
    rng = MersenneTwister(2718)
    S = 60
    C = 0.1
    K0,ζ = 5.0,1.0
    μ0,σ0 = 1.0,0.5
    alphas = [0.0,0.25,0.5,0.75,1.0]
    nrep = 10

    errs = Float64[]
    xs = Float64[]
    for α in alphas
        for rep = 1:nrep
            edges = gen_trophic_edges_cascade(S,C,α,rng)
            A = build_A(S,edges,μ0,σ0,rng)
            φs = φ_sim(A;K0,ζ,S,rng=rng)
            φp = φ_pred(A;K0,ζ)
            push!(errs, abs(φp-φs))
            push!(xs, α)
        end
    end

    # plot
    fig = Figure(; size=(500,350))
    ax = Axis(fig[1,1],
        xlabel="Ordering α",
        ylabel="|φ_pred − φ_sim|",
        title="Single-group prediction error vs cascade ordering"
    )
    scatter!(ax, xs, errs; color=xs, colormap=:blues, marker=:circle)
    display(fig)
end

sweep_alpha()

