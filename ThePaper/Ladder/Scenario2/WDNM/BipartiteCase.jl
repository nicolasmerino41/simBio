# Function to build J = D * M with M = -I + random G scaled by σ
# Build a bipartite Jacobian with power-law off-diagonals
function build_bipartite_J(d::Vector{Float64}, S::Int, R::Int; σ=0.1, pareto_exponent=2.5)
    # S total species, first R are resources, rest are consumers
    C = S - R
    D = diagm(d)
    # Initialize interaction matrix M = -I
    M = -I(S)
    # Bipartite: consumers (j>R) eat resources (i<=R)
    for i in 1:R, j in R+1:S
        # sample power-law magnitude
        strength = rand(Pareto(1.0, pareto_exponent)) * σ
        M[j,i] += strength   # consumer gain
        M[i,j] -= strength   # resource loss
    end
    # scale other entries (non-bipartite) with small Gaussian noise
    for i in 1:S, j in 1:S
        if i!=j && !(i<=R && j>R) && !(i>R && j<=R)
            M[i,j] += randn()*σ*0.1
        end
    end
    return D * M
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

##################################################################################
##################################################################################
##################################################################################
####################### TEST_ALL_BIPARTITE #######################################
# generalize_visualize_bipartite_vs_demography_makie.jl
# Like test_all, but using bipartite heavy-tailed Jacobians

# Build a bipartite Jacobian with power-law off-diagonals
function build_bipartite_J(d::Vector{Float64}, S::Int, R::Int; σ=0.1, pareto_exponent=2.5)
    C = S - R
    D = diagm(d)
    M = -Matrix{Float64}(I, S, S)
    # Consumer-resource links
    for i in 1:R, j in R+1:S
        strength = rand(Pareto(1.0, pareto_exponent)) * σ
        M[j,i] += strength   # consumer gain
        M[i,j] -= strength   # resource loss
    end
    # small noise elsewhere
    for i in 1:S, j in 1:S
        if i!=j && !((i<=R && j>R) || (i>R && j<=R))
            M[i,j] += randn() * σ * 0.1
        end
    end
    return D * M
end

# Metrics
function measure_collectivity_bip(d,S,R;σ=0.1,pareto_exponent=2.5)
    # spectral radius of the underlying bipartite A
    A = zeros(S,S)
    for i in 1:R, j in R+1:S
        A[j,i] = rand(Pareto(1.0,pareto_exponent))*σ
    end
    return maximum(abs, eigvals(A))
end

# Sensitivity of non-normality
function test_nonnormality_sensitivity_bip(d,σs,S,R;reps=200)
    Δnn=Float64[]
    for σ in σs
        buf=Float64[]
        for _ in 1:reps
            J0 = build_bipartite_J(d,S,R;σ=σ)
            J1 = build_bipartite_J(d,S,R;σ=σ)
            push!(buf, abs(nonnormality(J1)-nonnormality(J0)))
        end
        push!(Δnn, mean(buf))
    end
    return Δnn
end

function test_all_bipartite(S=50, R=20, reps=200)
    Random.seed!(123)
    d = rand(0.5:0.1:2.0, S)
    min_d = minimum(d)
    ratios = 10 .^ range(-2,1; length=15)
    σs = ratios .* min_d

    Δρ=Float64[]; Δrea=Float64[]; ΔRmed=Float64[]
    Δprt=Float64[]; Δpers=Float64[]; Δphi=Float64[]

    for σ in σs
        bufρ=Float64[]; bufrea=Float64[]; bufRmed=Float64[]
        bufprt=Float64[]; bufpers=Float64[]; bufphi=Float64[]

        for _ in 1:reps
            J0 = build_bipartite_J(d,S,R;σ=σ)
            J1 = build_bipartite_J(d,S,R;σ=σ)
            push!(bufρ, abs(measure_resilience(J1)-measure_resilience(J0)))
            push!(bufrea, abs(measure_reactivity(J1)-measure_reactivity(J0)))
            push!(bufRmed,abs(analytical_Rmed(J1)-analytical_Rmed(J0)))
            push!(bufprt, abs(pulse_return_time(J1)-pulse_return_time(J0)))
            push!(bufpers,abs(measure_persistence(J1)-measure_persistence(J0)))
            push!(bufphi, abs(measure_collectivity_bip(d,S,R;σ=σ) - measure_collectivity_bip(d,S,R;σ=σ)))
        end
        push!(Δρ,    mean(bufρ))
        push!(Δrea,  mean(bufrea))
        push!(ΔRmed, mean(bufRmed))
        push!(Δprt,  mean(bufprt))
        push!(Δpers, mean(bufpers))
        push!(Δphi,  mean(bufphi))
    end

    Δnn = test_nonnormality_sensitivity_bip(d,σs,S,R;reps=reps)

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
    #####################################################

    fig = Figure(; size=(1000,460))
    titles = ["Δ Resilience","Δ Reactivity","Δ Rₘₑ𝒹",
              "Δ Pulse RT","Δ Persistence","Δ Collectivity","Δ Nonnormality"]
    data   = [Δρ,Δrea,ΔRmed,Δprt,Δpers,Δphi,Δnn]
    colors = [:blue,:green,:purple,:orange,:brown,:gray,:teal]
    Label(fig[0, 1:3], "Bipartite Case"; fontsize=18)
    for (i,(t,Δ,c)) in enumerate(zip(titles,data,colors))
        row=(i-1)÷4+1; col=(i-1)%4+1
        ax = Axis(
            fig[row,col];
            xlabel="σ / min(d) (log)", xscale=log10,
            ylabel=t, title=t, titlesize=10
        )
        lines!(ax, ratios, Δ; color=c, linewidth=2)
        scatter!(ax, ratios, Δ; color=c, markersize=6)
    end
    display(fig)
end

test_all_bipartite()

