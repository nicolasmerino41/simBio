using Random, LinearAlgebra, Statistics, CairoMakie, Distributions

# Build J from a given interaction matrix A and diagonal rates d
function build_J_from_A(d::Vector{Float64}, A::Matrix{Float64})
    D = diagm(d)
    M = -I(size(A,1)) + A
    return D * M
end

# Build bipartite interaction A
function make_bipartite_A(S, R; σ=0.1, pareto_exponent=2.5)
    A = zeros(S,S)
    C = S - R
    for i in 1:R, j in R+1:S
        s = rand(Pareto(1, pareto_exponent)) * σ
        A[j,i] += s   # consumer gain
        A[i,j] -= s   # resource loss
    end
    return A
end

function test_all_bipartite_fromA(S=50, R=20, reps=200)
    Random.seed!(123)
    # sample self-regulation rates
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)
    ratios = 10 .^ range(-2,1; length=15)
    σs = ratios .* min_d

    # pre-draw K for SL
    K = abs.(randn(S) .+ 1)
    sparsity = 0.1
    mask = rand(S,S) .< sparsity
    mask .= mask .& .!I(S)

    # storage
    Δρ=Float64[]; Δrea=Float64[]; ΔRmed=Float64[]
    Δprt=Float64[]; Δpers=Float64[]; Δsl=Float64[]; Δphi=Float64[]

    for σ in σs
        bufρ=Float64[]; bufrea=Float64[]; bufRmed=Float64[]
        bufprt=Float64[]; bufpers=Float64[]; bufsl=Float64[]; bufphi=Float64[]

        # shared bipartite A
        A_shared = make_bipartite_A(S, R; σ=σ)

        for _ in 1:reps
            # perturb A_shared to A0 and A1
            A0 = copy(A_shared)
            A1 = copy(A_shared)
            # random small rewiring on non‐bipartite links
            for i in 1:S, j in 1:S
                if mask[i,j]
                    A1[i,j] = randn()*σ
                end
            end

            # build Jacobians
            J0 = build_J_from_A(d, A0)
            J1 = build_J_from_A(d, A1)

            # bulk metrics
            push!(bufρ,    abs(measure_resilience(J1)   - measure_resilience(J0)))
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
        push!(Δsl,   mean(bufsl))
        push!(Δphi,  mean(bufphi))
    end

    # nonnormality
    Δnn = test_nonnormality_sensitivity(d, σs; replicates=reps)

    # plot
    fig = Figure(; size=(1200,500))
    titles = ["Δ Resilience","Δ Reactivity","Δ Rₘₑ𝒹",
              "Δ Pulse RT","Δ Persistence","Δ Collectivity","Δ Nonnormality", "Δ SL"]
    data   = [Δρ,Δrea,ΔRmed,Δprt,Δpers,Δphi,Δnn,Δsl]
    colors = [:blue,:green,:purple,:orange,:grey,:brown,:yellow,:teal]
    Label(fig[0, 1:3], "Bipartite Case"; fontsize=18)
    for (i,(t,Δ,c)) in enumerate(zip(titles,data,colors))
        row=(i-1)÷4+1; col=(i-1)%4+1
        ax = Axis(fig[row,col];
                  xlabel="σ/min(d) (log)", xscale=log10,
                  yscale=t in ["Δ SL", "Δ Collectivity"] ? log10 : identity,
                  ylabel=t, title=t, titlesize=10)
        lines!(ax, ratios, Δ; color=c, linewidth=2)
        scatter!(ax, ratios, Δ; color=c, markersize=6)
    end

    display(fig)
end

test_all_bipartite_fromA()
