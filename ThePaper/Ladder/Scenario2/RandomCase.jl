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

# 3) analytical Rmed
function analytical_Rmed(J; t=0.01)
    E = exp(t*J)
    M = E*E'
    rates = diag(M) .- 1 .|> x-> -x/(2t)
    mean(rates)
end

# 4) pulse return‐time approx = 1/resilience
pulse_return_time(J) = 1/measure_resilience(J)

# 5) persistence approx = fraction of negative real‐parts
measure_persistence(J) = count(λ->real(λ)<0, eigvals(J))/length(eigvals(J))

# 6) collectivity = spectral radius of G part
measure_collectivity(d,S;σ=0.1) = begin
    # reconstruct A'=G (zero mean)
    G = randn(S,S)*σ
    for i in 1:S; G[i,i]=0; end
    maximum(abs, eigvals(G))
end

# Henrici non-normality
function nonnormality(J::AbstractMatrix)
    C = J*J' - J'*J
    return norm(C)
end

function test_nonnormality_sensitivity(d, σs; replicates=200)
    Δnn = Float64[]
    for σ in σs
        ds = Float64[]
        for _ in 1:replicates
            J0 = build_J(d, length(d); σ=σ)
            J1 = build_J(d, length(d); σ=σ)
            push!(ds, abs(nonnormality(J1) - nonnormality(J0)))
        end
        push!(Δnn, mean(ds))
    end
    return Δnn
end

function test_all(S=50, reps=200)
    Random.seed!(123)
    # Sample diagonals
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)
    ratios = 10 .^ range(-2, 1; length=15)
    σs = ratios .* min_d

    # Storage for Δ‐metrics
    Δρ = Float64[]; Δrea = Float64[]; ΔRmed = Float64[]
    Δprt = Float64[]; Δpers = Float64[]; Δphi = Float64[]

    for σ in σs
        bufρ=Float64[]; bufrea=Float64[]; bufRmed=Float64[] 
        bufprt=Float64[]; bufpers=Float64[]; bufphi=Float64[]

        for _ in 1:reps
            J0 = build_J(d, S; σ=σ)
            J1 = build_J(d, S; σ=σ)

            push!(bufρ,    abs(abs(measure_resilience(J1))   - abs(measure_resilience(J0))))
            push!(bufrea,  abs(measure_reactivity(J1)    - measure_reactivity(J0)))
            push!(bufRmed, abs(analytical_Rmed(J1)       - analytical_Rmed(J0)))
            push!(bufprt,  abs(pulse_return_time(J1)     - pulse_return_time(J0)))
            push!(bufpers, abs(measure_persistence(J1)   - measure_persistence(J0)))
            push!(bufphi,  abs(measure_collectivity(d,S;σ=σ) - measure_collectivity(d,S;σ=σ)))
        end

        push!(Δρ,    mean(bufρ))
        push!(Δrea,  mean(bufrea))
        push!(ΔRmed, mean(bufRmed))
        push!(Δprt,  mean(bufprt))
        push!(Δpers, mean(bufpers))
        push!(Δphi,  mean(bufphi))
    end

    # Compute non‐normality sensitivity
    Δnn = test_nonnormality_sensitivity(d, σs; replicates=reps)
    
    ###################### SL change ######################
    mean_sl_change = Float64[]
    # Random K and base A mask
    sparsity = 0.1
    K = abs.(randn(S) .+ 1.0)
    K = abs.(rand(Pareto(2.0), S) .+ 1.0)
    mask = rand(S, S) .< sparsity
    mask .= mask .& .!I(S)

    for σ in σs

        # SL-change: average abs diff
        sl_diffs = Float64[]
        for _ in 1:100
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
    ######################################################################
    
    # Plot all eight metrics in a 2×4 grid
    fig = Figure(; size = (1200, 450))
    titles = ["Δ Resilience", "Δ Reactivity", "Δ Rₘₑ𝒹",
              "Δ Pulse RT", "Δ Persistence", "Δ Collectivity", "Δ Nonnormality",
              "Δ SL"]
    data   = [Δρ, Δrea, ΔRmed, Δprt, Δpers, Δphi, Δnn, mean_sl_change]
    colors = [:blue, :green, :purple, :orange, :brown, :gray, :teal, :yellow]
    Label(fig[0, 1:3], "Random Case"; fontsize=18)
    for (i, (title, Δ, c)) in enumerate(zip(titles, data, colors))
        row = (i - 1) ÷ 4 + 1
        col = (i - 1) % 4 + 1
        ax = Axis(fig[row, col];
                  xlabel = "σ / min(d) (log)",
                  xscale  = log10,
                  yscale  = title == "Δ Pulse RT" ? log10 : identity,
                  ylabel = title,
                  title = title,
                  titlesize = 10)
        lines!(ax, ratios, Δ; color = c, linewidth = 2)
        scatter!(ax, ratios, Δ; color = c, markersize = 6)
    end

    display(fig)
end

test_all()
