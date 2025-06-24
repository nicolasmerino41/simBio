using Random, Statistics, LinearAlgebra
using DifferentialEquations
using QuadGK, Roots, Distributions
using CairoMakie

# --------------------------------------------------
# (Reuse or redefine the helper functions as before)
# --------------------------------------------------

# 1. analytical_phi
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

# 2. ER / BA (commented BA for future use)
function gen_ER_directed(S::Int, p::Float64, rng::AbstractRNG)
    edges = Tuple{Int,Int}[]
    for i in 1:S-1
        for j in i+1:S
            if rand(rng) < p
                if rand(rng) < 0.5
                    push!(edges, (i,j))
                else
                    push!(edges, (j,i))
                end
            end
        end
    end
    return edges
end

function gen_BA_directed(S::Int, m::Int, rng::AbstractRNG)
    @assert m >= 1 && m < S "m must be ≥1 and < S"
    adj = [Set{Int}() for _ in 1:S]
    for i in 1:m-1, j in i+1:m
        push!(adj[i], j); push!(adj[j], i)
    end
    deg = zeros(Int, S)
    for i in 1:m
        deg[i] = m-1
    end
    for k in m+1:S
        total = sum(deg[1:k-1])
        chosen = Set{Int}()
        while length(chosen) < m
            if total > 0
                r = rand(rng) * total
                cum = 0.0
                for i in 1:k-1
                    cum += deg[i]
                    if r < cum
                        push!(chosen, i)
                        break
                    end
                end
            else
                i = rand(rng, 1:k-1)
                push!(chosen, i)
            end
        end
        for i in chosen
            push!(adj[k], i); push!(adj[i], k)
            deg[k] += 1; deg[i] += 1
        end
    end
    directed_edges = Tuple{Int,Int}[]
    for i in 1:S
        for j in adj[i]
            if i < j
                if rand(rng) < 0.5
                    push!(directed_edges, (i,j))
                else
                    push!(directed_edges, (j,i))
                end
            end
        end
    end
    return directed_edges
end

# 3. Power-law directed via configuration model
function gen_powerlaw_directed(S::Int, γ::Float64, kmin::Int, rng::AbstractRNG)
    deg = zeros(Int, S)
    # continuous power-law sampling between kmin and (S-1)
    for i in 1:S
        u = rand(rng)
        # avoid γ=1 exactly; assume γ != 1
        # sample kf in [kmin, S-1] with density ∝ k^{-γ}
        # Inverse transform: CDF ~ (k^{1-γ} - kmin^{1-γ})/( (S-1)^{1-γ} - kmin^{1-γ} )
        # Solve for k:
        denom = (S-1)^(1-γ) - kmin^(1-γ)
        # if denom <= 0 (e.g. S-1 == kmin), just set deg = kmin
        if denom <= 0
            ki = kmin
        else
            kf = ( u*denom + kmin^(1-γ) )^(1/(1-γ))
            ki = clamp(floor(Int, kf), kmin, S-1)
        end
        deg[i] = ki
    end
    # ensure even sum
    if isodd(sum(deg))
        deg[rand(rng, 1:S)] += 1
    end
    # build stubs
    stubs = Int[]
    for i in 1:S
        append!(stubs, fill(i, deg[i]))
    end
    shuffle!(stubs)
    edges = Tuple{Int,Int}[]
    while length(stubs) >= 2
        i = pop!(stubs)
        j = pop!(stubs)
        if i != j
            push!(edges, (i,j))
        else
            # discard self-loop; could reinsert stub, but for simplicity discard
        end
    end
    # assign random direction
    directed = Tuple{Int,Int}[]
    for (i,j) in edges
        if rand(rng) < 0.5
            push!(directed, (i,j))
        else
            push!(directed, (j,i))
        end
    end
    return directed
end

# 4. Build interaction matrix
function build_interaction_matrix(S::Int, edges::Vector{Tuple{Int,Int}}, μ0::Float64, σ0::Float64)
    A = zeros(Float64, S, S)
    for (i,j) in edges
        # i->j means i prey, j predator
        A[j, i] = rand() * sqrt((σ0^2)/S) + (μ0 / S)
        A[i, j] = rand() * sqrt((σ0^2)/S) - (μ0 / S)
    end
    return A
end

# 5. gLV simulation + single-group prediction
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for idx in eachindex(u)
        du[idx] = u[idx] * (K[idx] - u[idx] - Au[idx])
    end
end

function simulate_and_predict(A; K0=5.0, ζ=1.0, t_end=200.0, rng=GLOBAL_RNG)
    S = size(A, 1)
    K = randn(rng, S) .* ζ .+ K0
    K .= max.(K, 1e-3)
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    N_end = sol.u[end]
    N_end .= max.(N_end, 0.0)
    φ_sim = sum(N_end .> 1e-3) / S
    # global moments
    off = Float64[]
    for i in 1:S, j in 1:S
        if i != j
            push!(off, A[i,j])
        end
    end
    meanA = mean(off)
    varA = var(off)
    μ = S * meanA
    σ = sqrt(S * varA)
    φ_pred = analytical_phi(K0, ζ, μ, σ)
    return φ_sim, φ_pred
end

# --------------------------------------------------
# 6. Experiment: now sweep PL exponents AND interaction strengths
# --------------------------------------------------
function experiment_PL_and_IS_sweep()
    rng = MersenneTwister(42)
    S = 200
    K0 = 5.0
    ζ = 1.0
    σ0 = 0.5   # keep σ0 fixed; sweep mean μ0
    nrep = 10

    # Define range of power-law exponents γ
    gammas = collect(2.0:0.3:3.4)  # example: 2.0, 2.3, 2.6, 2.9, 3.2, 3.4
    # You may extend below 2.0 or above 3.4 as you wish.

    # Define range of mean interaction strengths μ0 to sweep
    # mu_values = range(0.0, stop=2.0, length=11)  # e.g., 0,0.2,...,2.0
    mu_values = range(0.0, stop=50.0, length=11)  # adjust as needed

    # Storage: For each γ, store vectors of (mu_values, φ_sim_mean, φ_sim_std, φ_pred_mean)
    results = Dict{Float64, Dict{Symbol, Vector{Float64}}}()
    # results[γ] = Dict(:mu => [...], :φ_sim_mean => [...], :φ_sim_std => [...], :φ_pred_mean => [...])

    for γ in gammas
        # initialize
        results[γ] = Dict(:mu => Float64[], :φ_sim_mean => Float64[], :φ_sim_std => Float64[], :φ_pred_mean => Float64[])
        for μ0 in mu_values
            φ_sims = Float64[]
            φ_preds = Float64[]
            for rep in 1:nrep
                # generate PL network with exponent γ, kmin=2
                edges = gen_powerlaw_directed(S, γ, 2, rng)
                # compute degree heterogeneity optionally
                # deg = zeros(Int, S)
                # for (i,j) in edges
                #     deg[i] += 1; deg[j] += 1
                # end
                # cv = std(deg)/mean(deg)
                # build interaction matrix with mean μ0
                A = build_interaction_matrix(S, edges, μ0, σ0)
                φ_sim, φ_pred = simulate_and_predict(A; K0=K0, ζ=ζ, rng=rng)
                push!(φ_sims, φ_sim)
                push!(φ_preds, φ_pred)
            end
            push!(results[γ][:mu], μ0)
            push!(results[γ][:φ_sim_mean], mean(φ_sims))
            push!(results[γ][:φ_sim_std], std(φ_sims))
            push!(results[γ][:φ_pred_mean], mean(φ_preds))
        end
    end

    # Plot: for each γ, plot a curve φ_sim ± error vs μ0, and overlay φ_pred curve
    fig = Figure(; size=(800,600))
    nγ = length(gammas)
    # Choose a color scheme for different γ
    colors = Colors.distinguishable_colors(nγ)
    ax = Axis(fig[1,1];
        xlabel="Mean interaction strength μ0",
        ylabel="Surviving fraction φ",
        title="φ_sim ± std and φ_pred vs μ0 for different PL exponents γ",
    )
    for (idx, γ) in enumerate(gammas)
        dat = results[γ]
        mu_vec = dat[:mu]
        sim_mean = dat[:φ_sim_mean]
        sim_std = dat[:φ_sim_std]
        pred_mean = dat[:φ_pred_mean]
        c = colors[idx]
        # error bars for simulated
        for i in 1:length(mu_vec)
            x = mu_vec[i]
            y = sim_mean[i]
            yerr = sim_std[i]
            lines!(ax, [x, x], [y-yerr, y+yerr]; color=c, linestyle=:dash)
        end
        # scatter simulated means
        scatter!(ax, mu_vec, sim_mean; color=c, marker=:circle, label="sim γ=$(round(γ,digits=2))")
        # overlay predicted line
        lines!(ax, mu_vec, pred_mean; color=c, linewidth=2, linestyle=:solid, label="pred γ=$(round(γ,digits=2))")
    end
    axislegend(ax; position=:lb, nbanks=2)
    display(fig)
end

# Run the combined PL × IS sweep experiment
experiment_PL_and_IS_sweep()
