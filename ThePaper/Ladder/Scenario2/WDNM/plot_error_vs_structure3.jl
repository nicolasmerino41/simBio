using GLM, StatsModels

function plot_lm_error_vs_structure_metrics(
    df::DataFrame;
    steps=1:5,
    remove_unstable::Bool=false,
    save_plot::Bool=false,
    outlier_quantile=nothing,
    resolution = (900, 650)
)
    # 1) Filter out runs whose resilience at any step is unstable
    df_local = deepcopy(df)
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df_local = filter(row -> all(row[c] < 0 for c in rescols), df_local)
    end

    # 2) Compute structural properties
    N = nrow(df_local)
    conn      = Vector{Float64}(undef, N)
    mdeg      = Vector{Float64}(undef, N)
    modu      = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)
    nest      = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df_local))
        A   = row.p_final[2]
        Adj = A .!= 0.0
        S   = size(A, 1)

        g = SimpleGraph(A .!= 0)
        degs = degree(g)
        degree_cv[i] = std(degs) / mean(degs)

        conn[i]      = sum(Adj) / (S*(S-1))
        mdeg[i]      = mean(sum(Adj, dims=2))

        # modularity
        k = sum(Adj, dims=2)[:]
        m = sum(k) / 2
        B = zeros(Float64, S, S)
        for u in 1:S, v in 1:S
            B[u,v] = (Adj[u,v] ? 1.0 : 0.0) - (k[u]*k[v])/(2m)
        end
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        svec = map(x -> x >= 0 ? 1.0 : -1.0, v1)
        modu[i] = (svec' * (B * svec)) / (4m)

        # nestedness
        nested_sum = 0.0
        nested_count = 0
        for u in 1:S, v in u+1:S
            du, dv = sum(Adj[u,:]), sum(Adj[v,:])
            denom = max(du, dv)
            if denom > 0
                nested_sum  += sum(Adj[u,:] .& Adj[v,:]) / denom
                nested_count += 1
            end
        end
        nest[i] = nested_count > 0 ? nested_sum / nested_count : NaN
    end

    df_local.connectance = conn
    df_local.mean_degree = mdeg
    df_local.modularity  = modu
    df_local.nestedness  = nest
    df_local.degree_cv   = degree_cv

    # 3) Define properties and system metrics
    props   = [:connectance, :modularity, :degree_cv]
    props_names = ["Connectance", "Modularity", "Degree CV"]
    metrics = [:resilience, :reactivity, :rt_pulse, :after_press]
    titles  = Dict(
        :resilience  => "Resilience",
        :reactivity  => "Reactivity",
        :rt_pulse    => "Return Time",
        :after_press => "Persistence"
    )

    fig = Figure(; size = resolution)
    colors = [:red, :blue, :green, :orange, :purple]
    step_names = ["Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS", "Changing groups"]

    for (pi, p) in enumerate(props), (mi, m) in enumerate(metrics)
        ax = Axis(fig[mi, pi];
            title  = "$(titles[m]) Vs $(props_names[pi])",
            xlabel = "$(props_names[pi])",
            ylabel = "Relative error",
            xgridvisible = false,
            ygridvisible = false
        )

        for (si, s) in enumerate(steps)
            scol   = Symbol("$(m)_S$(s)")
            fullco = Symbol("$(m)_full")

            # compute relative error
            e   = abs.(df_local[!, scol] .- df_local[!, fullco])
            e   = e ./ (abs.(df_local[!, fullco]) .+ 1e-6)
            errs = (e .+ 1e-6) ./ (1 + 2e-6)
            xs = df_local[!, p]

            if outlier_quantile !== nothing
                keep = errs .<= quantile(errs, outlier_quantile)
                xs   = xs[keep]
                errs = errs[keep]
            end

            # Scatter plot
            scatter!(ax, xs, errs; color=colors[si], alpha=0.3, markersize=3)

            # Fit linear regression
            df_tmp = DataFrame(y=errs, x=xs)
            model = lm(@formula(y ~ x), df_tmp)
            xfit = range(minimum(xs), maximum(xs), length=100)
            yfit = coef(model)[1] .+ coef(model)[2] .* xfit

            # Plot regression line
            lines!(ax, xfit, yfit; color=colors[si], linewidth=2, label=step_names[si])
        end

        if pi == 3 && mi == 1
            axislegend(ax; position=:lt)
        end
    end

    display(fig)
    if save_plot
        save("lm_error_vs_metrics_rel.png", fig; px_per_unit=6.0)
    end

    return df_local
end

plot_lm_error_vs_structure_metrics(
    G;
    steps=[1, 2, 3, 5],
    remove_unstable=true,
    save_plot=false,
    outlier_quantile=0.9
)