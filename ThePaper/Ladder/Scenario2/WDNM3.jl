using Random, LinearAlgebra, Statistics, CairoMakie

# (same build_J, measure_resilience, compute_SL from before) …

function experiment_grainy(S=50, replicates=500)
    Random.seed!(123)
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)

    # Log-spaced ratios from 0.01 to 10
    ratios = 10 .^ range(-2, 1; length=15)
    σs = ratios .* min_d

    mean_sl_change = Float64[]

    # fixed sparsity mask & K
    sparsity = 0.1
    K = abs.(randn(S) .+ 1.0)
    mask = rand(S, S) .< sparsity
    mask .= mask .& .!I(S)

    for σ in σs
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

    # Grainier SL‐sensitivity plot
    begin
        fig = Figure(; size=(700,400))
        ax = Axis(fig[1,1];
                  xlabel = "σ / min(d) (log scale)",
                  ylabel = "Mean |Δ SL|",
                  xscale = log10,
                  title = "Network Sensitivity of SL vs Interaction Strength")
        lines!(ax, ratios, mean_sl_change; color=:red, linewidth=2)
        scatter!(ax, ratios, mean_sl_change; color=:red, markersize=8)
        fig |> display
    end
end

experiment_grainy()
