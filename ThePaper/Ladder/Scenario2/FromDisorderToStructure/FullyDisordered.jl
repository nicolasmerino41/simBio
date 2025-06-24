# Requires: DifferentialEquations.jl, Random, Statistics
using DifferentialEquations, Random, Statistics

# 1. Parameters for the reference model
S = 100                      # number of species
μ = 1.0                      # target μ = S * mean(alpha_ij)
σ = 0.5                      # target σ: sqrt(S * Var(alpha_ij))
ζ = 1.0                      # standard deviation of carrying capacities K_i around mean K0
K0 = 5.0                     # mean carrying capacity
rng = MersenneTwister(1234)

# Function to generate a random interaction matrix with specified moments:
# Here we do a Gaussian ensemble: alpha_ij ~ Normal(mean=μ/S, var=(σ^2)/S), independent,
# and set diagonal zero (self-regulation handled separately via K_i and negative N_i^2 term).
function gen_alpha_gaussian(S, μ, σ; rng=rng)
    m = μ / S
    v = σ^2 / S
    A = randn(rng, S, S) .* sqrt(v) .+ m
    # Optionally enforce zero diagonal:
    for i in 1:S
        A[i,i] = 0.0
    end
    return A
end

# Another ensemble: e.g., a binary ± interaction but scaled to match same mean/var:
# For demonstration, generate entries equal to a or b with probabilities p and 1-p,
# choose p, a, b so that mean = μ/S and var = σ^2/S.
function gen_alpha_binary(S, μ, σ; rng=rng)
    # Solve for p, a, b: we'll pick symmetrical ±c structure:
    # Let alpha_ij ∈ {+c, -c} with probability p, 1-p. Then mean = c*(2p-1) = μ/S;
    # Var = E[alpha^2] - mean^2 = c^2 - (μ/S)^2 = σ^2/S.
    # => c^2 = σ^2/S + (μ/S)^2. Then p = ( (μ/S)/c + 1 )/2.
    c = sqrt(σ^2/S + (μ/S)^2)
    p = ((μ/S)/c + 1)/2
    A = zeros(Float64, S, S)
    for i in 1:S, j in 1:S
        if i==j
            A[i,i] = 0.0
        else
            A[i,j] = rand(rng) < p ? c : -c
        end
    end
    return A
end

# 2. Carrying capacities K_i: normal around K0 with sd ζ
function gen_K(S, K0, ζ; rng=rng)
    return randn(rng, S) .* ζ .+ K0
end

# 3. Define gLV ODE: dN_i/dt = N_i*(K_i - N_i - sum_j alpha_ij * N_j)
function make_gLV!(du, u, p, t)
    # u: vector of abundances N_i
    # p: tuple (K, A)
    K, A = p
    # compute interaction term: A * u
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] - Au[i])
    end
end

# 4. Simulate dynamics to equilibrium given initial conditions
function simulate_equilibrium(K, A; tspan=(0.0, 200.0), reltol=1e-6, abstol=1e-6)
    S = length(K)
    u0 = rand(rng, S) .* K  # random positive initial abundances, e.g., uniform [0, K_i]
    prob = ODEProblem(make_gLV!, u0, tspan, (K, A))
    sol = solve(prob, Tsit5(), reltol=reltol, abstol=abstol)
    N_end = sol(tspan[end])
    # truncate negatives or tiny values:
    N_end .= max.(N_end, 0.0)
    return N_end
end

# 5. Compute assembly outcomes: surviving fraction, total biomass, mean abundance, etc.
function summarize_outcome(N_end; thresh=1e-3)
    survivors = sum(N_end .> thresh)
    frac = survivors / length(N_end)
    total_biomass = sum(N_end)
    mean_abund = survivors > 0 ? mean(N_end[N_end .> thresh]) : 0.0
    return (frac=frac, total=total_biomass, mean_survivor=mean_abund)
end

# 6. Run two ensembles with same moments but different distributions
# Gaussian
A1 = gen_alpha_gaussian(S, μ, σ, rng=rng)
K1 = gen_K(S, K0, ζ, rng=rng)
N1 = simulate_equilibrium(K1, A1)
out1 = summarize_outcome(N1)

# Binary ±
A2 = gen_alpha_binary(S, μ, σ, rng=rng)
K2 = gen_K(S, K0, ζ, rng=rng)   # use same distribution for K_i; or fix K2=K1 to isolate A differences
# To isolate effect of interaction structure only, reuse the same K:
K2 = K1
N2 = simulate_equilibrium(K2, A2)
out2 = summarize_outcome(N2)

println("Outcome Gaussian ensemble: ", out1)
println("Outcome Binary ± ensemble: ", out2)

# 7. Compare: they should be similar if Barbier’s disordered-limit holds.
# You can repeat multiple replicates (with different RNG seeds) and gather statistics:
function replicate_outcomes(nrep)
    results = Vector{NamedTuple{(:frac1,:total1,:frac2,:total2),Tuple{Float64,Float64,Float64,Float64}}}(undef, nrep)
    for rep in 1:nrep
        # optionally reseed or use different rng seeds
        A1 = gen_alpha_gaussian(S, μ, σ, rng=rng)
        A2 = gen_alpha_binary(S, μ, σ, rng=rng)
        K = gen_K(S, K0, ζ, rng=rng)
        N1 = simulate_equilibrium(K, A1)
        N2 = simulate_equilibrium(K, A2)
        s1 = summarize_outcome(N1)
        s2 = summarize_outcome(N2)
        results[rep] = (s1.frac, s1.total, s2.frac, s2.total)
    end
    return results
end

# Example: run 10 replicates
res = replicate_outcomes(10)
println("Replicates (frac1, total1, frac2, total2):")
for r in res
    println(r)
end

# 8. (Optional) Compute aggregate moments from each generated A to confirm they match:
function compute_moments(A, K)
    S = size(A,1)
    # Collect off-diagonal entries via comprehension
    offdiag = [A[i,j] for i in 1:S for j in 1:S if i != j]
    meanA = mean(offdiag)
    varA  = var(offdiag)
    μ_est = S * meanA
    σ_est = sqrt(S * varA)
    ζ_est = std(K)
    # Reciprocity γ: correlation of A[i,j] and A[j,i] over i<j
    pairs = [(A[i,j], A[j,i]) for i in 1:S for j in i+1:S]
    if !isempty(pairs)
        v1 = [p[1] for p in pairs]
        v2 = [p[2] for p in pairs]
        γ_est = cor(v1, v2)
    else
        γ_est = NaN
    end
    return (μ_est=μ_est, σ_est=σ_est, ζ_est=ζ_est, γ_est=γ_est)
end


println("Moments A1: ", compute_moments(A1, K1))
println("Moments A2: ", compute_moments(A2, K2))

