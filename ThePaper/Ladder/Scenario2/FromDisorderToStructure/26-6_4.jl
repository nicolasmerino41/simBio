using Random, Statistics, LinearAlgebra
using Distributions, DifferentialEquations
using QuadGK, Roots, NLsolve
using CairoMakie

# 1. gLV dynamics
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] - Au[i])
    end
end

# 2. Simulate to equilibrium → φ_sim (bipartite, group‐specific K)
function simulate_equilibrium_phi_bip(A; K0_b, ζ_b, K0_c, ζ_c, S_b, t_end=200.0, rng=GLOBAL_RNG)
    S = size(A,1)
    K = Vector{Float64}(undef, S)
    for i in 1:S
        if i ≤ S_b
            K[i] = max(randn(rng)*ζ_b + K0_b, 1e-6)
        else
            K[i] = max(randn(rng)*ζ_c + K0_c, 1e-6)
        end
    end
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol  = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
    Nend = max.(sol.u[end], 0.0)
    return mean(Nend .> 1e-3)
end

# 3. Single‐group analytic φ
function analytical_phi(K0, ζ, μ, σ; tol=1e-6)
    stdn = Normal(0,1)
    function RHS(φ)
        φc = clamp(φ,0,1)
        integrand(z) = (1 - cdf(stdn,(μ*φc + σ*sqrt(φc)*z - K0)/ζ)) * pdf(stdn,z)
        val,_ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
        return val
    end
    f(φ) = RHS(φ) - φ
    if f(1e-8) < 0    return 0.0
    elseif f(1-1e-8) > 0  return 1.0
    else                 return find_zero(f,(1e-8,1-1e-8),Bisection(); xtol=tol)
    end
end

function φ_pred_single(A; K0_b, ζ_b, K0_c, ζ_c, S_b)
    S = size(A,1)
    # global K0, ζ
    S_c = S - S_b
    K0_glob = (S_b*K0_b + S_c*K0_c)/S
    E_K2 = (S_b*(ζ_b^2 + K0_b^2) + S_c*(ζ_c^2 + K0_c^2)) / S
    ζ_glob = sqrt(max(E_K2 - K0_glob^2, 1e-8))
    # A moments
    off = [A[i,j] for i in 1:S, j in 1:S if i!=j]
    μ_glob = S*mean(off)
    σ_glob = sqrt(S*var(off))
    return analytical_phi(K0_glob, ζ_glob, μ_glob, σ_glob)
end

# 4. Two‐group analytic φ_b, φ_c
function analytical_phi_two_groups(K0_b, ζ_b, K0_c, ζ_c,
                                   μ_bb, μ_bc, μ_cb, μ_cc,
                                   σ_bb, σ_bc, σ_cb, σ_cc;
                                   tol=1e-6)
    stdn = Normal(0,1)
    function F!(Fout, φ)
        φb, φc = clamp.(φ, 0,1)
        mb = μ_bb*φb + μ_bc*φc
        mc = μ_cb*φb + μ_cc*φc
        vb = σ_bb^2*φb + σ_bc^2*φc
        vc = σ_cb^2*φb + σ_cc^2*φc
        sb, sc = sqrt(vb), sqrt(vc)
        integrand_b(z) = (1 - cdf(stdn,(mb + sb*z - K0_b)/ζ_b))*pdf(stdn,z)
        integrand_c(z) = (1 - cdf(stdn,(mc + sc*z - K0_c)/ζ_c))*pdf(stdn,z)
        RHSb,_ = quadgk(integrand_b, -Inf, Inf; rtol=1e-6)
        RHSc,_ = quadgk(integrand_c, -Inf, Inf; rtol=1e-6)
        Fout[1] = RHSb - φb
        Fout[2] = RHSc - φc
    end
    sol = nlsolve(F!, [0.5,0.5]; xtol=tol, ftol=tol)
    return clamp.(sol.zero, 0,1)
end

# 5. Bipartite edge generator
function gen_bipartite_edges(S_b, S_c, C, α, rng)
    edges = Tuple{Int,Int}[]
    for i in 1:S_b, j in (S_b+1):(S_b+S_c)
        if rand(rng) < C
            if rand(rng) < α
                push!(edges,(i,j))
            else
                rand(rng)<0.5 ? push!(edges,(i,j)) : push!(edges,(j,i))
            end
        end
    end
    return edges
end

# 6. Build A from edges
function build_A(S, edges, μ0, σ0, rng)
    A = zeros(S,S)
    for (i,j) in edges
        m = rand(rng, truncated(Normal(μ0/S, sqrt(σ0^2/S)), 0, Inf))
        A[j,i] = m; A[i,j] = -m
    end
    return A
end

# 7. Moments for two‐group
function group_moments(A, grp, g, h)
    vals = Float64[]
    for i in 1:size(A,1), j in 1:size(A,2)
        if i!=j && grp[i]==g && grp[j]==h
            push!(vals, A[i,j])
        end
    end
    if isempty(vals)
        return 0.0, 0.0
    else
        S = size(A,1)
        return mean(vals)*S, sqrt(var(vals)*S)
    end
end

# 8. Main sweep and plot
function main()
    rng = MersenneTwister(2025)
    # parameters
    S, S_b, S_c, C = 60, 30, 30, 0.1
    K0_b, ζ_b = 5.0, 1.0
    K0_c, ζ_c = 1.0, 1.0
    σ0 = 0.5

    alphas = [0.0,0.25,0.5,0.75,1.0]
    mu0s   = [1.0,5.0,10.0]
    nrep   = 5

    φ_sim_mat   = zeros(length(alphas), length(mu0s))
    φ_pred1_mat = zeros(length(alphas), length(mu0s))
    φ_pred2_mat = zeros(length(alphas), length(mu0s))

    for (iα, α) in enumerate(alphas), (iμ, μ0) in enumerate(mu0s)
        sims, p1s, p2s = Float64[], Float64[], Float64[]
        for rep in 1:nrep
            edges = gen_bipartite_edges(S_b, S_c, C, α, rng)
            A     = build_A(S, edges, μ0, σ0, rng)

            # φ_sim
            push!(sims, simulate_equilibrium_phi_bip(A;
                K0_b=K0_b,ζ_b=ζ_b, K0_c=K0_c,ζ_c=ζ_c,
                S_b=S_b, rng=rng))

            # single‐group pred
            push!(p1s, φ_pred_single(A;
                K0_b=K0_b,ζ_b=ζ_b, K0_c=K0_c,ζ_c=ζ_c,
                S_b=S_b))

            # two‐group pred
            # infer groups trivially:
            grp = vcat(fill(1,S_b),fill(2,S_c))
            μ_bb,σ_bb = group_moments(A,grp,1,1)
            μ_bc,σ_bc = group_moments(A,grp,1,2)
            μ_cb,σ_cb = group_moments(A,grp,2,1)
            μ_cc,σ_cc = group_moments(A,grp,2,2)
            φb,φc = analytical_phi_two_groups(
                K0_b,ζ_b, K0_c,ζ_c,
                μ_bb,μ_bc,μ_cb,μ_cc,
                σ_bb,σ_bc,σ_cb,σ_cc
            )
            push!(p2s, (S_b*φb + S_c*φc)/S)
        end

        φ_sim_mat[iα,iμ]   = mean(sims)
        φ_pred1_mat[iα,iμ] = mean(p1s)
        φ_pred2_mat[iα,iμ] = mean(p2s)
    end

    # plotting
    fig = Figure(; size=(900,350))

    ax1 = Axis(fig[1,1]; xlabel="α", ylabel="φ", title="φ vs ordering α")
    for iμ in 1:length(mu0s)
        μ0 = mu0s[iμ]
        lines!(ax1, alphas, φ_pred1_mat[:,iμ], color=:gray, linestyle=:dash, label="1‐group μ₀=$(μ0)")
        lines!(ax1, alphas, φ_pred2_mat[:,iμ], color=:black, linestyle=:solid, label="2‐group μ₀=$(μ0)")
        scatter!(ax1, alphas, φ_sim_mat[:,iμ], marker=:circle, label="sim μ₀=$(μ0)")
    end
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1,2]; xlabel="μ₀", ylabel="φ", title="φ vs interaction μ₀")
    for iα in 1:length(alphas)
        α = alphas[iα]
        lines!(ax2, mu0s, φ_pred1_mat[iα,:], color=:gray, linestyle=:dash, label="1‐group α=$(α)")
        lines!(ax2, mu0s, φ_pred2_mat[iα,:], color=:black, linestyle=:solid, label="2‐group α=$(α)")
        scatter!(ax2, mu0s, φ_sim_mat[iα,:], marker=:diamond, label="sim α=$(α)")
    end
    axislegend(ax2; position=:rt)

    display(fig)
end

# Run the experiment
main()
