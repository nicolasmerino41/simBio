using Random, LinearAlgebra, Statistics, CairoMakie

# Build J directly from a given A and diagonal d
function build_J_from_A(d::Vector{Float64}, A::Matrix{Float64})
    D = diagm(d)
    M = -I(size(A,1)) + A
    return D * M
end

# Compute collectivity from a given A
function measure_collectivity_from_A(_d, _S, A::Matrix{Float64})
    return maximum(abs, eigvals(A))
end

function test_all(S=50, reps=200)
    Random.seed!(123)
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)
    ratios = 10 .^ range(-2, 1; length=15)
    σs = ratios .* min_d

    # Pre-generate K and mask
    sparsity = 0.9
    K = abs.(randn(S) .+ 1.0)
    mask = rand(S, S) .< sparsity
    mask .= mask .& .!I(S)

    # Storage
    Δρ=Float64[]; Δrea=Float64[]; ΔRmed=Float64[]
    Δprt=Float64[]; Δpers=Float64[]; Δphi=Float64[]; Δsl=Float64[]

    for σ in σs
        # build a shared A for this σ
        A_shared = zeros(S, S)
        A_shared[mask] .= randn(sum(mask)) .* σ

        bufρ=Float64[]; bufrea=Float64[]; bufRmed=Float64[]
        bufprt=Float64[]; bufpers=Float64[]; bufphi=Float64[]; bufsl=Float64[]

        for _ in 1:reps
            # perturb A_shared to A0 and A1
            A0 = copy(A_shared)
            A1 = copy(A_shared)
            A1[mask] .= randn(sum(mask)) .* σ

            # build Jacobians
            J0 = build_J_from_A(d, A0)
            J1 = build_J_from_A(d, A1)

            # bulk metrics
            push!(bufρ,    abs(abs(measure_resilience(J1))   - abs(measure_resilience(J0))))
            push!(bufrea,  abs(measure_reactivity(J1)    - measure_reactivity(J0)))
            push!(bufRmed, abs(analytical_Rmed(J1)       - analytical_Rmed(J0)))
            push!(bufprt,  abs(pulse_return_time(J1)     - pulse_return_time(J0)))
            push!(bufpers, abs(measure_persistence(J1)   - measure_persistence(J0)))
            push!(bufphi,  abs(measure_collectivity_from_A(d,S,A1) - measure_collectivity_from_A(d,S,A0)))

            # SL change
            sl0 = compute_SL(A0, K)
            sl1 = compute_SL(A1, K)
            push!(bufsl, mean(abs.(sl1 - sl0)))
        end

        push!(Δρ,    mean(bufρ))
        push!(Δrea,  mean(bufrea))
        push!(ΔRmed, mean(bufRmed))
        push!(Δprt,  mean(bufprt))
        push!(Δpers, mean(bufpers))
        push!(Δphi,  mean(bufphi))
        push!(Δsl,   mean(bufsl))
    end

    # Compute non‐normality sensitivity
    Δnn = test_nonnormality_sensitivity(d, σs; replicates=reps)

    # Plot
    fig = Figure(; size=(1200, 450))
    titles = ["Δ Resilience", "Δ Reactivity", "Δ Rₘₑ𝒹",
              "Δ Pulse RT", "Δ Persistence", "Δ Collectivity", "Δnonnormality", "Δ SL"]
    data   = [Δρ, Δrea, ΔRmed, Δprt, Δpers, Δnn, Δphi, Δsl]
    colors = [:blue, :green, :purple, :orange, :brown, :gray, :yellow, :darkgreen]
    Label(fig[0, 1:3], "Random Case"; fontsize=18)
    for (i, (title, Δ, c)) in enumerate(zip(titles, data, colors))
        row = (i-1) ÷ 4 + 1
        col = (i-1) % 4 + 1
        ax = Axis(fig[row, col];
                  xlabel="σ / min(d) (log)", xscale=log10,
                  yscale  = title == "Δ Pulse RT" ? log10 : identity,
                  ylabel=title, title=title, titlesize=10)
        lines!(ax, ratios, Δ; color=c, linewidth=2)
        scatter!(ax, ratios, Δ; color=c, markersize=6)
    end

    display(fig)
end

test_all()
