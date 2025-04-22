df = new_results10
# define the steps in order
step_keys = ["Full"; ["S$(i)" for i in 1:15]...]  # adjust the 15 if you have more/less steps

# the metrics you want to plot
metrics = [
    :return_time, 
    :overshoot, 
    :ire, 
    :compound_error,
    :after_persistence
]

# for each metric, build a figure
for metric in metrics
    n = length(step_keys[1:end])
    fig = Figure(; size = (1000, 1000), fontsize=12)
    for i in 1:div(n, 4)
        
        for j in 1:4
            s = step_keys[j + (i-1)*4]
            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(df.connectance)))
            color_indices = [findfirst(==(val), sort(unique(df.connectance))) for val in df.connectance]
            color_vals = cmap[color_indices]
            ax = Axis(
                fig[i, j]; title = s,
                xlabel = "Degree CV",
                ylabel = string(metric),
                # title  = "Step $s"
            )
            # scatter degree_cv vs metric_s
            scatter!(
                ax, 
                df.degree_cv, 
                df[!, Symbol("$(metric)_$s")],
                markersize = 6,
                color = df.connectance,
                colormap = :viridis,
                colorrange = (minimum(df.connectance), maximum(df.connectance))
            )
            # add a little 1:1 reference line?
            # you could also compute & annotate a correlation here
        end
    end
    # fig.layoutgap = (10,10)
    display(fig)
end


begin
    using Statistics, DataFrames, CairoMakie

    df = new_results10   # your big DataFrame
    step_keys = ["Full"; ["S$(i)" for i in 1:15]...]  # adjust to your number of steps
    METS = [
        "return_time", 
        "overshoot", 
        "ire", 
        "compound_error",
        "after_persistence"
    ]

    # 1) Prepare storage
    steps_col   = String[]
    metrics_col = String[]
    r_degree_cv = Float64[]

    # 2) Loop and compute correlations
    for s in step_keys, met in METS
        vec_x = df.degree_cv
        vec_y = df[!, Symbol("$(met)_$s")]
        ok = .!ismissing.(vec_y) .& .!isnan.(vec_y)
        if count(ok) > 5  # require at least a handful of points
            r = cor(vec_x[ok], vec_y[ok])
        else
            r = NaN
        end
        push!(steps_col, s)
        push!(metrics_col, met)
        push!(r_degree_cv, r)
    end

    # 3) Build summary DataFrame
    summary = DataFrame(
        step        = steps_col,
        metric      = metrics_col,
        r_degree_cv = r_degree_cv
    )

    # 4) Plot one line per metric
    fig = Figure(resolution = (800, 200 * length(METS)))
    for (i, met) in enumerate(METS)
        ax = Axis(fig[i, 1];
            title = met,
            xlabel = "Simplification step",
            ylabel = "cor(degree_CV, $met)"
        )
        dat = filter(row -> row.metric == met, summary)
        xs = 1:length(step_keys)
        ys = dat.r_degree_cv
        lines!(ax, xs, ys)
        # axisxticks!(ax, xs, step_keys)
        # rotate_xticks!(ax, 45)
    end
    fig
end

begin
    using DataFrames, GLM, Statistics, CairoMakie

    df = new_results10
    step_keys = ["Full"; ["S$(i)" for i in 1:15]...]  # adjust to your number of steps
    METS = [
        "return_time",
        "overshoot",
        "ire",
        "compound_error",
        "after_persistence"
    ]

    # For each metric, build a separate forest plot of β₁,s ± 95% CI
    for met in METS
        forest_steps = String[]
        betas        = Float64[]
        lo95         = Float64[]
        hi95         = Float64[]

        # fit linear model y_met,s ~ 1 + degree_cv for each step s
        for s in step_keys
            y = df[!, Symbol("$(met)_$s")]
            x = df.degree_cv
            ok = .!ismissing.(y) .& .!isnan.(y)
            if count(ok) < 5
                continue
            end
            sub = DataFrame(y = y[ok], x = x[ok])
            lm_m = lm(@formula(y ~ x), sub)

            # extract slope and its standard error correctly from the fitted model
            β1   = coef(lm_m)[2]
            se1  = stderror(lm_m)[2]
            ci_lo = β1 - 1.96 * se1
            ci_hi = β1 + 1.96 * se1

            push!(forest_steps, s)
            push!(betas, β1)
            push!(lo95, ci_lo)
            push!(hi95, ci_hi)
        end

        forest_df = DataFrame(
            step  = forest_steps,
            β1    = betas,
            lo95  = lo95,
            hi95  = hi95,
        )

        # plot
        n = nrow(forest_df)
        fig = Figure(resolution = (600, 20 * n))
        ax  = Axis(fig[1,1];
            title     = "Effect of degree_CV on $met",
            xlabel    = "Slope β₁ (degree_CV → $met)",
            ylabel    = "Step",
            yreversed = true,
            yticks    = (1:n, forest_df.step)
        )
        for i in 1:n
            y = i
            x = forest_df.β1[i]
            l = forest_df.lo95[i]
            h = forest_df.hi95[i]
            lines!(ax, [l, h], [y, y], linewidth=2)
            scatter!(ax, [x], [y], markersize=6)
        end
        hlines!(ax, 1:n, color=(:gray, 0.3))
        vlines!(ax, [0], color=(:black, 0.5), linestyle=:dash)
        display(fig)
    end
end

begin
    #
    # 3) Heatmap of R² for each (step, metric)
    #
    r2_mat = fill(NaN, length(step_keys), length(METS))
    for (j, s) in enumerate(step_keys), (i, met) in enumerate(METS)
        y = df[!, Symbol("$(met)_$s")]
        x = df.degree_cv
        ok = .!ismissing.(y) .& .!isnan.(y)
        if count(ok) < 5
            continue
        end
        sub = DataFrame(y = y[ok], x = x[ok])
        lm_m = lm(@formula(y ~ x), sub)
        r2_mat[j, i] = r2(lm_m)
    end

    fig2 = Figure(resolution = (800, 400))
    ax2 = Axis(fig2[1,1];
        xlabel     = "Metric",
        ylabel     = "Step",
        xticks     = (1:length(METS), METS),
        yticks     = (1:length(step_keys), step_keys),
        yreversed  = true,
        title      = "Heatmap of R² (degree_CV → metric)"
    )
    heatmap!(ax2, r2_mat; colormap = :viridis, colorrange = (0,1))
    # Colorbar(fig2[1,2], ax2)
    display(fig2)
end