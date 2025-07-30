function plot_error_vs_structure_binned_lines(
    df::DataFrame,
    var::Symbol;
    steps::AbstractVector{Int}=1:5,
    remove_unstable::Bool=false,
    n_bins::Int=50,
    separate::Bool=false,
    save_plot::Bool=false,
    relative_error::Bool=false,
    error_bars::Bool=true
)
    # Optionally filter unstable runs
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # --- Compute network structural properties ---
    N = nrow(df)
    conn = Vector{Float64}(undef, N)
    mdeg = Vector{Float64}(undef, N)
    modu = Vector{Float64}(undef, N)
    nest = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]
        Adj = A .!= 0.0
        S = size(A, 1)

        conn[i] = sum(Adj) / (S*(S-1))
        mdeg[i] = mean(sum(Adj, dims=2))
        degree_cv[i] = std(sum(Adj, dims=2)) / mdeg[i]

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

    # Attach to DataFrame
    df.connectance = conn
    df.mean_degree = mdeg
    df.modularity  = modu
    df.nestedness  = nest
    df.degree_cv   = degree_cv

    # Properties and full-error column
    props = [:connectance, :modularity, :degree_cv]
    fullcol = Symbol("$(var)_full")

    # Setup figure
    if separate
        # grid: rows = length(props), cols = length(steps)
        n_rows = length(props)
        n_cols = length(steps)
        fig = Figure(; size=(300 * n_cols, 250 * n_rows))

        for (ri, p) in enumerate(props)
            for (ci, s) in enumerate(steps)
                ax = Axis(fig[ri, ci];
                    title = "$(var) S$s - $(string(p)) vs |error|",
                    xlabel = string(p),
                    ylabel = "|error|"
                )

                # compute errors
                scol = Symbol("$(var)_S$(s)")
                errs = abs.(df[!, scol] .- df[!, fullcol])
                if relative_error
                    rel_errs = errs ./ (abs.(df[!, fullcol]) .+ 1e-6)
                    errs = (rel_errs .+ 1e-6) ./ (1 + 2e-6)
                end
                xs = df[!, p]

                # binning
                xmin, xmax = minimum(xs), maximum(xs)
                edges = range(xmin, xmax, length=n_bins+1)
                bins = searchsortedlast.(Ref(edges), xs)

                bmx, bmy, bsy = Float64[], Float64[], Float64[]
                for b in 1:n_bins
                    idxs = findall(bins .== b)
                    if !isempty(idxs)
                        push!(bmx, mean(xs[idxs]))
                        push!(bmy, mean(errs[idxs]))
                        push!(bsy, std(errs[idxs]))
                    end
                end

                # plot
                lines!(ax, bmx, bmy; color=:blue, linewidth=2)
                if error_bars
                    errorbars!(ax, bmx, bmy, bsy; color=:blue)
                end
            end
        end
    else
        # one column of plots: one per property, all steps overlaid
        fig = Figure(; size=(1100, 700))
        colors = [:red, :blue, :green, :orange, :purple]

        for (idx, p) in enumerate(props)
            for (s, title) in enumerate(["Absolute error", "Relative error"])
                ax = Axis(fig[idx, s];
                    title = "$(var): $(string(p)) vs $(title)",
                    xlabel = string(p),
                    ylabel = title == "Relative error" ? "% error" : "|error|"
                )

                for (si, s) in enumerate(steps)
                    scol = Symbol("$(var)_S$(s)")
                    errs = abs.(df[!, scol] .- df[!, fullcol])
                    if title == "Relative error"
                        rel_errs = errs ./ (abs.(df[!, fullcol]) .+ 1e-6)
                        errs = (rel_errs .+ 1e-6) ./ (1 + 2e-6)
                    end
                    xs = df[!, p]

                    xmin, xmax = minimum(xs), maximum(xs)
                    edges = range(xmin, xmax, length=n_bins+1)
                    bins = searchsortedlast.(Ref(edges), xs)

                    bmx, bmy, bsy = Float64[], Float64[], Float64[]
                    for b in 1:n_bins
                        idxs = findall(bins .== b)
                        if !isempty(idxs)
                            push!(bmx, mean(xs[idxs]))
                            push!(bmy, mean(errs[idxs]))
                            push!(bsy, std(errs[idxs]))
                        end
                    end

                    lines!(ax, bmx, bmy;
                        color=colors[si], linewidth=2, label="Step $s"
                    )
                    if error_bars
                        errorbars!(ax, bmx, bmy, bsy;
                            color=colors[si]
                        )
                    end
                end
                # axislegend(ax; position=:lt)
            end
        end
    end

    display(fig)
    if save_plot
        save("errorVsStructure_$(var).png", fig; px_per_unit=6.0)
    end

    return df
end

for i in [:resilience, :reactivity, :rt_pulse, :after_press]
    for j in [false]
        fig = plot_error_vs_structure_binned_lines(
            G, 
            i;
            steps=[1, 2, 3, 5],
            remove_unstable=true,
            save_plot=true,
            separate=false,
            n_bins=100,
            relative_error=j,
            error_bars=false
        )
    end 
end

"""
    plot_error_vs_structure_metrics_binned_lines(
        df::DataFrame;
        steps::AbstractVector{Int}=1:5,
        remove_unstable::Bool=false,
        n_bins::Int=50,
        save_plot::Bool=false,
        error_bars::Bool=true
    )

For each of the four system‐level metrics (Resilience, Reactivity,
Return Time, Persistence), plots the *relative* error at each
step against each of the three structural properties
(connectance, modularity, degree_cv).  The result is a 3×4 grid:
rows = structural properties, columns = metrics.
"""
function plot_error_vs_structure_metrics_binned_lines(
    df::DataFrame;
    steps=1:5,
    remove_unstable::Bool=false,
    n_bins::Int=50,
    save_plot::Bool=false,
    error_bars::Bool=true,
    outlier_quantile  = nothing,
    outlier_quantile_x = 0.9
)
    # 1) Filter out runs whose resilience at any step is unstable
    if remove_unstable
        rescols = Symbol.(string(:resilience) .* "_S" .* string.(steps))
        df = filter(row -> all(row[c] < 0 for c in rescols), df)
    end

    # 2) Compute the three structural properties
    N = nrow(df)
    conn      = Vector{Float64}(undef, N)
    mdeg      = Vector{Float64}(undef, N)
    modu      = Vector{Float64}(undef, N)
    degree_cv = Vector{Float64}(undef, N)
    nest      = Vector{Float64}(undef, N)  # if you ever want nestedness

    for (i, row) in enumerate(eachrow(df))
        A   = row.p_final[2]
        Adj = A .!= 0.0
        S   = size(A, 1)

        g = SimpleGraph(A .!= 0)         # unweighted graph
        degs = degree(g)
        degree_cv[i] = std(degs) / mean(degs)

        conn[i]      = sum(Adj) / (S*(S-1))
        mdeg[i]      = mean(sum(Adj, dims=2))

        # modularity (Newman’s leading eigenvector)
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

        # (optional) nestedness
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

    df.connectance = conn
    df.mean_degree = mdeg
    df.modularity  = modu
    df.nestedness  = nest
    df.degree_cv   = degree_cv

    # 3) Which structural props & which system metrics?
    props   = [:connectance, :modularity, :degree_cv]
    props_names = ["Connectance", "Modularity", "Degree CV"]
    metrics = [:resilience, :reactivity, :rt_pulse, :after_press]
    titles  = Dict(
        :resilience  => "Resilience",
        :reactivity  => "Reactivity",
        :rt_pulse    => "Return Time",
        :after_press => "Persistence"
    )

    # 4) Prepare figure: rows = props, cols = metrics
    nP = length(props)
    nM = length(metrics)
    fig = Figure(; size = (1100, 1100))
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

        # overlay one line per step
        for (si, s) in enumerate(steps)
            scol   = Symbol("$(m)_S$(s)")
            fullco = Symbol("$(m)_full")

            # compute squeezed relative error
            e   = abs.(df[!, scol] .- df[!, fullco])
            e   = e ./ (abs.(df[!, fullco]) .+ 1e-6)
            errs = (e .+ 1e-6) ./ (1 + 2e-6)

            xs = df[!, p]
            xs[xs .<= 0.0] .= 0.0

            if outlier_quantile !== nothing
                thresh = quantile(errs, outlier_quantile)
                keep   = errs .<= thresh
                xs     = xs[keep]
                errs   = errs[keep]
                thresh = quantile(xs, outlier_quantile_x)
                keep   = xs .<= thresh
                xs     = xs[keep]
                errs   = errs[keep]
            end
            println("max metric value is: for metric $(p) and step $(s): $(maximum(xs))")
            # bin into n_bins
            xmin, xmax = minimum(xs), maximum(xs)
            edges = range(xmin, xmax, length = n_bins+1)
            bix   = searchsortedlast.(Ref(edges), xs)

            mx, my, sy = Float64[], Float64[], Float64[]
            for b in 1:n_bins
                idxs = findall(bix .== b)
                if !isempty(idxs)
                    push!(mx, mean(xs[idxs]))
                    push!(my, mean(errs[idxs]))
                    push!(sy, std(errs[idxs]))
                end
            end

            lines!(ax, mx, my;
                color     = colors[si],
                linewidth = 2,
                label     = step_names[si]
            )
            if error_bars
                errorbars!(ax, mx, my, sy; color = colors[si])
            end
        end

        # show legend only in top-left panel
        # Show legend only in top-left panel, with smaller size
        if pi == 1 && mi == 1
            axislegend(
                ax;
                position = :rt,
                # labelsize = 11,     # font size of legend labels
                # markersize = 9,     # size of legend markers
                # patchsize = (10, 10) # size of legend box patches
            )
end

    end

    display(fig)
    if save_plot
        save("error_vs_metrics_rel.png", fig; px_per_unit=6.0)
    end

    return df
end

for i in [30]
    plot_error_vs_structure_metrics_binned_lines(
        G;
        steps=[1, 2, 3, 5],
        remove_unstable=false,
        n_bins=i,
        save_plot=true,
        error_bars=false,
        outlier_quantile=0.9,
        outlier_quantile_x=1.0
    )
end