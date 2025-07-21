function plot_error_vs_structure(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:7,
    remove_unstable::Bool=false
)
    # 1) optional filter on stability
    if remove_unstable
        rescols = Symbol.(string(var) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # 2) extract structural props
    N = nrow(df)
    conn = Vector{Float64}(undef, N)
    mdeg = Vector{Float64}(undef, N)
    modu = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[8]            # S×S interaction matrix
        # build unweighted adjacency
        Adj = A .!= 0.0
        G = SimpleGraph(Adj)

        # connectance = edges / [S*(S-1)]
        S = size(A,1)
        conn[i] = ne(G) / (S*(S-1))

        # mean degree
        mdeg[i] = mean(degree(G))

        # —— spectral modularity —— 
        # build modularity matrix B = A_uw - (k kᵀ)/(2m)
        k = sum(Adj, dims=2)[:]              # degree vector
        m = sum(k) / 2                       # number of edges
        B = zeros(Float64, S, S)
        for u in 1:S, v in 1:S
            B[u,v] = (Adj[u,v] ? 1.0 : 0.0) - (k[u]*k[v])/(2m)
        end

        # leading eigenvector of B (symmetric)
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        # partition by sign of v1 → s ∈ {+1, -1}
        s = map(x-> x >= 0 ?  1.0 : -1.0, v1)

        # Q = (1/(4m)) sᵀ B s
        modu[i] = (s' * (B * s)) / (4m)
    end

    df[!,:connectance] = conn
    df[!,:mean_degree] = mdeg
    df[!,:modularity]  = modu

    # 3) build long form of errors
    long = DataFrame(prop=String[], step=Int[], error=Float64[], value=Float64[])
    fullcol = Symbol(string(var)*"_full")

    for s in steps
        col = Symbol(string(var)*"_S$(s)")
        for row in eachrow(df)
            err = abs(row[col] - row[fullcol])
            push!(long, ("connectance", s, err, row.connectance))
            push!(long, ("mean_degree", s, err, row.mean_degree))
            push!(long, ("modularity", s, err, row.modularity))
        end
    end

    # 4) plot 3×1 scatter + trend
    fig = Figure(size=(1000,450))
    props  = ["connectance"]#,"mean_degree","modularity"]
    titles = ["Connectance"]#,"Mean degree","Modularity"]

    for (i, prop) in enumerate(props)
        ax = Axis(fig[i,1]; title=titles[i], xlabel=titles[i], ylabel="Absolute error")
        sub = filter(r -> r.prop==prop, long)

        for s in steps
            ss = filter(r -> r.step==s, sub)
            scatter!(ax, ss.value, ss.error;
                label="S$s", markersize=4, alpha=0.6)
            xs = sort(unique(ss.value))
            ys = [mean(filter(r->r.value==x, ss).error) for x in xs]
            lines!(ax, xs, ys; linestyle=:dash, color=:black)
        end

        axislegend(ax, position=:rt)
    end

    display(fig)
    return fig
end


plot_error_vs_structure(
    G, :resilience;
    # steps=1:7,
    remove_unstable=true
)

function plot_error_vs_structure_box(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:7,
    n_bins::Int=5,
    remove_unstable::Bool=false
)
    # 1) optional filter
    if remove_unstable
        rescols = Symbol.(string(var) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # 2) compute structural properties
    N = nrow(df)
    conn = similar(Vector{Float64}(undef, N))
    mdeg = similar(conn)
    modu = similar(conn)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[8]                  # interaction matrix
        Adj = A .!= 0.0
        G = SimpleGraph(Adj)

        S = size(A,1)
        conn[i] = ne(G) / (S*(S-1))         # connectance
        mdeg[i] = mean(degree(G))           # mean degree

        # spectral modularity
        k = sum(Adj, dims=2)[:]             # degrees
        m = sum(k) / 2
        B = Adj .- (k * k')/(2m)
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        signs = map(x -> x >= 0 ? 1.0 : -1.0, v1)
        modu[i] = (signs' * (B * signs)) / (4m)
    end

    df[!,:connectance] = conn
    df[!,:mean_degree] = mdeg
    df[!,:modularity]  = modu

    # 3) bin each property into quantile bins
    props = [:connectance, :mean_degree, :modularity]
    edges = Dict(p => quantile(df[!,p], range(0,1,length=n_bins+1)) for p in props)
    bins  = Dict{Symbol, Vector{Int}}()
    for p in props
        e = edges[p]
        bins[p] = [ min(searchsortedlast(e, v), n_bins) for v in df[!,p] ]
    end

    # 4) prepare figure: 7 rows × 3 cols
    fig = Figure(size=(900, 1100))
    fullcol = Symbol(string(var)*"_full")

    for (r, s) in enumerate(steps), (c, p) in enumerate(props)
        ax = Axis(fig[r, c];
            title = "S$(s) vs Full | $(p)",
            xlabel = "Bin",
            ylabel = "|error|"
        )

        # compute errors for this step
        scol = Symbol(string(var)*"_S$(s)")
        errs = abs.(df[!, scol] .- df[!, fullcol])
        bx   = bins[p]  # bin index per row

        # draw boxplot via scatter+stat summary
        # flatten xs, ys vectors
        boxplot!(ax, bx, errs)
        
        # label x-axis
        # ax.xticks = (1:n_bins,
        #              [ @sprintf("%.3f–%.3f", edges[p][i], edges[p][i+1]) for i in 1:n_bins ])
    end

    display(fig)
    return fig
end


fig = plot_error_vs_structure_box(
    G, :resilience;
    # steps=1:7,
    remove_unstable=true
)

filename = "errorVsStructureAndSteps.png"
save(filename, fig; px_per_unit=6.0)