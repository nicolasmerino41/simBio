# generalize_visualize_network_vs_demography_makie.jl
# 1) Generalize across replicates and σ values
# 2) Visualize with Makie: resilience vs σ/min_d and SL-change vs σ/min_d
# 3) Interpret implications in comments

using Random, LinearAlgebra, Statistics, CairoMakie

# Function to build J = D * M with M = -I + random G scaled by σ
function build_J(d::Vector{Float64}, S::Int; σ=0.1)
    D = diagm(d)
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

function measure_reactivity(J::Matrix{Float64})
    J_sym = (J + J') / 2
    ev_sym = eigvals(J_sym)
    return maximum(real.(ev_sym))
end

# Compute SL = B ./ K with B = (I - A) \ K
function compute_SL(A::Matrix{Float64}, K::Vector{Float64})
    B = (I - A) \ K
    return B ./ K
end

function experiment(S=50, replicates=100)
    Random.seed!(42)
    # Sample diagonals and find min
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)

    # σ grid as multiples of min_d
    ratios = [0.01, 0.1, 0.5, 1.0, 2.0]
    σs = ratios .* min_d

    # Storage
    mean_res = Float64[]
    mean_sl_change = Float64[]

    # Random K and base A mask
    sparsity = 0.1
    K = abs.(randn(S) .+ 1.0)
    mask = rand(S, S) .< sparsity
    mask .= mask .& .!I(S)

    for σ in σs
        # Bulk resilience
        res_vals = [measure_resilience(build_J(d, S; σ=σ)) for _ in 1:replicates]
        push!(mean_res, mean(res_vals))

        # SL-change: average abs diff
        sl_diffs = Float64[]
        for _ in 1:replicates
            A = zeros(S, S)
            A[mask] .= randn(sum(mask)) .* σ
            sl0 = compute_SL(A, K)
            A_rand = copy(A)
            A_rand[mask] .= randn(sum(mask)) .* σ
            sl1 = compute_SL(A_rand, K)
            push!(sl_diffs, mean(abs.(sl1 - sl0)))
        end
        push!(mean_sl_change, mean(sl_diffs))
    end

    # Plot 1: resilience vs σ/min_d
    begin
        fig = Figure(; size = (600, 400))
        ax = Axis(fig[1,1]; xlabel = "σ / min(d)", ylabel = "Mean resilience",
                  title = "Bulk Stability vs Interaction Strength")
        lines!(ax, ratios, mean_res; color = :blue, linewidth = 2)
        scatter!(ax, ratios, mean_res; color = :blue)
        display(fig)
    end

    # Plot 2: SL-change vs σ/min_d
    begin
        fig = Figure(; size = (600, 400))
        ax = Axis(fig[1,1]; xlabel = "σ / min(d)", ylabel = "Mean |Δ SL|",
                  title = "Network Sensitivity of SL vs Strength")
        lines!(ax, ratios, mean_sl_change; color = :red, linewidth = 2)
        scatter!(ax, ratios, mean_sl_change; color = :red)
        display(fig)
    end
end

experiment()
