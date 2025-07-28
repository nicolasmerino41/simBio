function plot_error_vs_structure(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:7,
    remove_unstable::Bool=false
)
    # 1) optional filter on stability
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # 2) extract structural props
    N = nrow(df)
    conn = Vector{Float64}(undef, N)
    mdeg = Vector{Float64}(undef, N)
    modu = Vector{Float64}(undef, N)
    nest = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]            # S×S interaction matrix
        # build unweighted adjacency
        Adj = A .!= 0.0
        G = SimpleGraph(Adj)

        # connectance = edges / [S*(S-1)]
        S = size(A,1)
        conn[i] = ne(G) / (S*(S-1))

        # mean degree
        mdeg[i] = mean(degree(G))
        cv = std(degree(G)) / mean(degree(G))

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

        nested_sum = 0.0
        nested_count = 0
        for u in 1:S, v in u+1:S
            du = sum(Adj[u, :])
            dv = sum(Adj[v, :])
            denom = max(du, dv)
            if denom > 0
                overlap = sum(Adj[u, :] .& Adj[v, :]) / denom
                nested_sum += overlap
                nested_count += 1
            end
        end
        nest[i] = nested_count > 0 ? nested_sum / nested_count : NaN

        # degree CV
        degree_cv[i] = cv
    end

    df[!,:connectance] = conn
    df[!,:mean_degree] = mdeg
    df[!,:modularity]  = modu
    df[!,:nestedness]  = nest
    df[!,:degree_cv]   = degree_cv

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
            push!(long, ("nestedness", s, err, row.nestedness))
            push!(long, ("degree_cv", s, err, row.degree_cv))
        end
    end

    # 4) plot 3×1 scatter + trend
    fig = Figure(size=(1000,450))
    props  = ["connectance", "modularity", "degree_cv"]#,"mean_degree","modularity"]
    titles = ["Connectance", "Modularity", "Degree CV"]#,"Mean degree","Modularity"]

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
    G, :reactivity;
    steps=[1, 2, 3, 5],
    remove_unstable=true
)

function plot_error_vs_structure_box(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:5,
    n_bins::Int=5,
    remove_unstable::Bool=false,
    save_plot::Bool=false
)
    # 1) optional filter
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # 2) extract structural props
    N = nrow(df)
    conn = Vector{Float64}(undef, N)
    mdeg = Vector{Float64}(undef, N)
    modu = Vector{Float64}(undef, N)
    nest = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]                  # interaction matrix
        Adj = A .!= 0.0
        G = SimpleGraph(Adj)

        S = size(A,1)
        conn[i] = ne(G) / (S*(S-1))         # connectance
        mdeg[i] = mean(degree(G))           # mean degree
        cv = std(degree(G)) / mean(degree(G))

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

        nested_sum = 0.0
        nested_count = 0
        for u in 1:S, v in u+1:S
            du = sum(Adj[u, :])
            dv = sum(Adj[v, :])
            denom = max(du, dv)
            if denom > 0
                overlap = sum(Adj[u, :] .& Adj[v, :]) / denom
                nested_sum += overlap
                nested_count += 1
            end
        end
        nest[i] = nested_count > 0 ? nested_sum / nested_count : NaN

        # degree CV
        degree_cv[i] = cv
    end

    df[!,:connectance] = conn
    df[!,:mean_degree] = mdeg
    df[!,:modularity]  = modu
    df[!,:nestedness]  = nest
    df[!,:degree_cv]   = degree_cv

    # 3) bin each property into quantile bins and print bin edges
    props = [:connectance, :modularity, :nestedness, :degree_cv]
    edges = Dict(p => quantile(df[!,p], range(0,1,length=n_bins+1)) for p in props)
    bins  = Dict{Symbol, Vector{Int}}()

    println("Bin ranges for each property:")
    for p in props
        e = edges[p]
        println("  $(p):")
        for i in 1:n_bins
            println("    Bin $i: $(round(e[i], digits=4)) – $(round(e[i+1], digits=4))")
        end
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
        # draw boxplot
        boxplot!(ax, bx, errs)

        # label x-axis
        # ax.xticks = (1:n_bins,
        #              [ @sprintf("%.3f–%.3f", edges[p][i], edges[p][i+1]) for i in 1:n_bins ])
    end

    display(fig)
    if save_plot
        filename = "errorVsStructureAndSteps_$(var).png"
        save(filename, fig; px_per_unit=6.0)
    end
    
end


fig = plot_error_vs_structure_box(
    G, :reactivity;
    steps=[1,2,3,5],
    remove_unstable=true,
    save_plot=false
)

function plot_error_vs_structure_points_continuous(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:5,
    remove_unstable::Bool=false,
    save_plot::Bool=false
)
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    N = nrow(df)
    conn = Vector{Float64}(undef, N)
    mdeg = Vector{Float64}(undef, N)
    modu = Vector{Float64}(undef, N)
    nest = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]
        Adj = A .!= 0.0
        G = SimpleGraph(Adj)
        S = size(A, 1)

        conn[i] = ne(G) / (S*(S-1))
        mdeg[i] = mean(degree(G))
        degree_cv[i] = std(degree(G)) / mdeg[i]

        k = sum(Adj, dims=2)[:]
        m = sum(k) / 2
        B = zeros(Float64, S, S)
        for u in 1:S, v in 1:S
            B[u,v] = (Adj[u,v] ? 1.0 : 0.0) - (k[u]*k[v]) / (2m)
        end
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        s = map(x -> x >= 0 ? 1.0 : -1.0, v1)
        modu[i] = (s' * (B * s)) / (4m)

        nested_sum = 0.0
        nested_count = 0
        for u in 1:S, v in u+1:S
            du, dv = sum(Adj[u,:]), sum(Adj[v,:])
            denom = max(du, dv)
            if denom > 0
                overlap = sum(Adj[u,:] .& Adj[v,:]) / denom
                nested_sum += overlap
                nested_count += 1
            end
        end
        nest[i] = nested_count > 0 ? nested_sum / nested_count : NaN
    end

    df[!,:connectance] = conn
    df[!,:mean_degree] = mdeg
    df[!,:modularity]  = modu
    df[!,:nestedness]  = nest
    df[!,:degree_cv]   = degree_cv

    props = [:connectance, :modularity, :nestedness, :degree_cv]
    fullcol = Symbol(string(var)*"_full")

    fig = Figure(size=(1000, 1000))
    for (r, s) in enumerate(steps), (c, p) in enumerate(props)
        ax = Axis(fig[r, c];
            title = "$var S$(s) vs Full | $(p)",
            xlabel = string(p),
            ylabel = "|error|"
        )

        scol = Symbol(string(var)*"_S$(s)")
        errs = abs.(df[!, scol] .- df[!, fullcol])
        xs = df[!, p]

        scatter!(ax, xs, errs; markersize=3, color=:black, alpha=0.3)
    end

    for (s) in steps
        scol = Symbol(string(var)*"_S$(s)")
        fullcol = Symbol(string(var)*"_full")

        raw_errs = df[!, scol] .- df[!, fullcol]
        rel_errs = raw_errs
        # rel_errs = abs.(raw_errs) ./ (abs.(df[!, fullcol]) .+ 1e-6)
        # rel_errs = (rel_errs .+ 1e-6) ./ (1 + 2e-6)  # squeeze into (0,1)
        
        println("=== Step S$s ===")
        for p in props
            df_tmp = DataFrame(y=rel_errs, x=df[!, p])
            model = lm(@formula(y ~ x), df_tmp)
            println("Metric: ", p)
            display(coeftable(model))
        end
    end

    display(fig)
    if save_plot
        filename = "errorVsStructureAndSteps_continuous_$(var).png"
        save(filename, fig; px_per_unit=6.0)
    end

    return df
end


T = G
G = T
G = filter(row -> row.IS > 0.01, G)

for met in [:resilience, :reactivity, :rt_pulse, :after_press]
    println("=== $(met) ===")
    fig = plot_error_vs_structure_points_continuous(
        G, 
        :resilience;
        steps=[1, 2, 3, 5],
        remove_unstable=true,
        save_plot=true
    )
end

df=fig
count(<(0.05), df.degree_cv)

using Plots

histogram(
    df.degree_cv;
    bins=100,
    xlabel="Degree CV",
    ylabel="Frequency",
    # title="Distribution of Degree CV across Networks"
)
