# network_vs_demography.jl
# Demonstrate:
# 1. Bulk stability (resilience) ≈ min(d_i) independent of off-diagonals when σ << min(d)
# 2. SL_i requires full network to compute
# 3. Threshold σ ≈ min(d) when network begins to matter

using LinearAlgebra, Random

# Function to build J = D * M with M = -I + A_scaled
function build_J(d::Vector{Float64}, S::Int; σ=0.1)
    # D: diagonal rates
    D = diagm(d)
    # random off-diagonal M: zero-mean normal entries
    M = -Matrix{Float64}(I, S, S)
    G = randn(S, S) * σ
    for i in 1:S, j in 1:S
        if i != j
            M[i, j] += G[i, j]
        end
    end
    return D * M
end

# Compute resilience = - max real part of eigenvalues
function measure_resilience(J::Matrix{Float64})
    vals = eigvals(J)
    return -maximum(real(vals))
end

# Compute SL for species i: B_i / K_i, with B = (I - A)^(-1) * K
function compute_SL(A::Matrix{Float64}, K::Vector{Float64})
    B = (I - A) \ K
    return B ./ K
end

# Main experiment
function main()
    S = 50
    Random.seed!(42)

    # Sample demographic diagonals d_i from uniform
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)

    # Range of σ values
    sigmas = vcat(0.01 * min_d, 0.1 * min_d, 0.5 * min_d, 1.0 * min_d, 2.0 * min_d)
    println("σ/mind | resilience(J) | min(d)")
    for σ in sigmas
        J = build_J(d, S; σ=σ)
        ρ = measure_resilience(J)
        "%.4f   |  %.4f      | %.4f, $(σ/min_d), $ρ, $min_d"
    end

    # SL demonstration
    # Build a random A with controlled σ on off-diagonals
    σ = 0.5 * min_d
    A = zeros(S, S)
    for i in 1:S, j in 1:S
        if i != j && rand() < 0.1
            A[i, j] = randn() * σ
        end
    end
    # Choose K randomly
    K = abs.(randn(S) .+ 1.0)
    sl = compute_SL(A, K)
    @show sl[1:5]

    # Show that shuffling A changes SL but not resilience when σ << min(d)
    println("\nEffect of randomizing A on SL vs resilience:")
    # Original J & resilience
    J0 = build_J(d, S; σ=0.1*min_d)
    r0 = measure_resilience(J0)
    println("Resilience initial: ", r0)
    # SL initial
    sl0 = compute_SL(A, K)

    # Randomize A
    A_rand = copy(A)
    for i in 1:S, j in 1:S
        if i != j && A[i,j] != 0.0
            A_rand[i,j] = randn() * σ
        end
    end
    sl_rand = compute_SL(A_rand, K)

    # New J resilience
    J1 = build_J(d, S; σ=0.1*min_d)
    r1 = measure_resilience(J1)
    println("Resilience after randomization: ", r1)
    println("SL change (first 5 species):", sl_rand[1:5] - sl0[1:5])
end

main()
