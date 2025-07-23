function plot_scalar_correlations(
    G::DataFrame;
    scenarios = [:ER, :PL, :MOD],
    metrics = [
        (:resilience, "Resilience"), (:reactivity, "Reactivity"),
        # (:mean_tau, "Mean SL"), 
        # (:mean_inverse_tau, "inverse_SL"),
        (:analytical_rmed, "Rmed"),
        (:rt_pulse, "Return Time"),
        (:after_persistence, "Persistence"),
        (:collectivity, "Collectivity"),
        # (:sigma_over_min_d, "σ/min(d)")
    ],
    fit_to_1_1_line::Bool = true,
    save_plot::Bool = false,
    resolution = (900, 650)
)
    step_names = ["sub_grouping", "Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS", "Rewiring + ↻C + ↻IS", "Not recalculating", "Changing groups"]
    for scen in scenarios
        df = G[G.scen .== scen, :]
        fig = Figure(; size=resolution)
        for (i, (sym, label)) in enumerate(metrics)
            for (j, step) in enumerate((2, 3, 4, 5, 7))
                x = df[!, Symbol(string(sym, "_full"))]
                y = df[!, Symbol(string(sym, "_S", step))]
                if step == 7
                    col_step = Symbol(string(sym, "_S7"))
                    col_full = Symbol(string(sym, "_full"))

                    # Ensure the column exists
                    if hasproperty(df, col_step)
                        # Filter only rows where the step-7 resilience is negative
                        new_df = filter(row -> row[:resilience_S7] < 0, df)

                        # Extract x and y vectors
                        x = new_df[!, col_full]
                        y = new_df[!, col_step]
                    else
                        @warn "Column $(col_step) not found in DataFrame."
                        x = Float64[]
                        y = Float64[]
                    end
                end
                mn = min(minimum(x), minimum(y))
                mx = max(maximum(x), maximum(y))

                ax = Axis(
                    fig[j, i];
                    # title="$(scen) $label: Full vs $(step_names[step])",
                    title="$label: Full vs $(step_names[step])",
                    # xlabel="Full", ylabel="step$(step)",
                    titlesize=9,
                    xlabelsize=10, ylabelsize=10,
                    xticklabelsize=10, yticklabelsize=10,
                    limits = ((mn, mx), (mn, mx)),
                    xgridvisible = false,
                    ygridvisible = false
                )

                scatter!(ax, x, y; alpha=0.3)

                # Plot 1:1 line
                lines!(
                    ax, [mn, mx], [mn, mx];
                    color=:black,
                    linestyle=:dash
                )

                if fit_to_1_1_line
                    # R² to the 1:1 line (y_hat = x)
                    y_hat = x
                    ss_tot = sum((y .- mean(y)).^2)
                    ss_res = sum((y .- y_hat).^2)
                    r2_1to1 = 1 - ss_res / ss_tot

                    text!(ax, "R²=$(round(r2_1to1, digits=3))";
                        position=(mx, mn),
                        align=(:right, :bottom),
                        fontsize=10,
                        color=:black
                    )
                else
                    # Pearson correlation coefficient
                    r_val = cor(x, y)
                    text!(ax, "r=$(round(r_val, digits=3))";
                        position = (mx, mn),
                        align    = (:right, :bottom),
                        fontsize = 10,
                        color = :black
                    )
                end
            end
        end
        if save_plot
            filename = "scalar_correlation_$(scen).png"
            save(filename, fig; px_per_unit=6.0)
        end
        display(fig)
    end
end

function plot_vector_correlations(
    df::DataFrame;
    scenarios=[:ER], # [:ER,:PL,:MOD],
    color_by=:conn,
    variable::Symbol = :tau,
    fit_to_1_1_line::Bool = true,
    save_plot::Bool = false,
    resolution = (900, 650),
    pixels_per_unit = 2.0
)
    # mapping from variable name to the full‐and‐step column suffix
    suffix = Dict(
        :tau     => "_full" => "_S",     # tau_full, tau_S1, tau_S2
        :inverse_tau => "_full" => "_S",
        :SL      => "_full" => "_S",
        :ssp_analytical_rmed => "_full" => "_S"
    )[variable]

    step_names = ["sub_grouping", "Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS", "Rewiring + ↻C + ↻IS", "Not recalculating", "Changing groups"]

    # color scale
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    for scen in scenarios
        sub = df[df.scen .== scen, :]
        fig = Figure(; size=resolution)
        # Label(fig[0, 1:4], uppercase(string(variable, " correlations")); fontsize=12)
        for (i, step) in enumerate((2, 3, 4, 5, 6, 7))
            col_full = Symbol(string(variable, suffix[1]))
            col_step = Symbol(string(variable, suffix[2], step))
            ax = Axis(fig[(i-1)÷3+1, (i-1)%3+1];
                # title  = "$(scen) $(variable): full vs S$step",
                title  = "$(step_names[step])",
                xlabel = "SL in full model", ylabel = "SL in step $step",
                titlesize=16,
                xlabelsize=14, ylabelsize=14,
                xticklabelsize=13, yticklabelsize=13,
                xgridvisible = false,
                ygridvisible = false
            )
            
            xs = Float64[]
            ys = Float64[]
            cs = Float64[]

            for row in eachrow(sub)
                v_full = row[col_full]      # e.g. Vector{Float64} or Vector{Point2}
                v_step = row[col_step]      # e.g. Vector{Float64}
                color_val = row[color_by]   # probably scalar
                
                n_points = min(length(v_full), length(v_step)) # just in case
                append!(xs, v_full[1:n_points])
                append!(ys, v_step[1:n_points])
                append!(cs, fill(color_val, n_points))
            end

            scatter!(ax, xs, ys;
                colormap   = :viridis,
                color      = cs,
                colorrange = (cmin, cmax),
                markersize = 4,
                alpha      = 0.6
            )

            mn = min(minimum(xs), minimum(ys))
            mx = max(maximum(xs), maximum(ys))
            lines!(ax, [mn,mx], [mn,mx]; color=:black, linestyle=:dash)

            if fit_to_1_1_line
                # R² to the 1:1 line (y_hat = x)
                y_hat = xs
                ss_tot = sum((ys .- mean(ys)).^2)
                ss_res = sum((ys .- y_hat).^2)
                r2_1to1 = 1 - ss_res / ss_tot

                text!(
                    ax, "R²=$(round(r2_1to1, digits=3))";
                    position=(mx, mn),
                    align=(:right, :bottom),
                    fontsize=10,
                    color=:black
                )
            else
                # Pearson correlation coefficient
                r_val = cor(xs, ys)
                text!(
                    ax, "r=$(round(r_val, digits=3))";
                    position = (mx, mn),
                    align    = (:right, :bottom),
                    fontsize = 10,
                    color = :black
                )
            end
        end

        # Colorbar(fig[1,3], colormap=:viridis,
        # label = string(color_by),
        # limits=(minimum(cs), maximum(cs))
        # )
        if save_plot
            filename = string("vector_correlations_", variable, "_", scen, ".png")
            save(filename, fig; px_per_unit=pixels_per_unit)
        end

        display(fig)
    end
end

function plot_vector_correlations(
    df::DataFrame;
    scenarios=[:ER],
    variable::Symbol = :tau,
    fit_to_1_1_line::Bool = true,
    save_plot::Bool = false,
    resolution = (900, 650),
    pixels_per_unit = 2.0
)
    suffix = Dict(
        :tau => ("_full", "_S"),
        :inverse_tau => ("_full", "_S"),
        :SL => ("_full", "_S"),
        :ssp_analytical_rmed => ("_full", "_S")
    )[variable]

    step_names = [
        "sub_grouping", "Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS",
        "Rewiring + ↻C + ↻IS", "Not recalculating", "Changing groups"
    ]

    for scen in scenarios
        sub = df[df.scen .== scen, :]
        fig = Figure(; size=resolution)

        for (i, step) in enumerate((2, 3, 4, 5, 6, 7))
            col_full = Symbol(string(variable, suffix[1]))
            col_step = Symbol(string(variable, suffix[2], step))

            # Special case for step 7 with additional filtering
            this_df = sub
            if step == 7
                col_step = Symbol(string(variable, "_S7"))
                col_full = Symbol(string(variable, "_full"))
                if hasproperty(df, col_step)
                    this_df = filter(row -> row[:resilience_S7] < 0 &&
                                             all(row[:tau_S7] .< 1000.0), sub)
                else
                    @warn "Column $(col_step) not found in DataFrame."
                    continue
                end
            end

            ax = Axis(fig[(i-1)÷3+1, (i-1)%3+1];
                title  = "$(step_names[step])",
                xlabel = "SL in full model", ylabel = "SL in step $step",
                titlesize=16,
                xlabelsize=14, ylabelsize=14,
                xticklabelsize=13, yticklabelsize=13,
                xgridvisible = false,
                ygridvisible = false
            )

            xs = Float64[]
            ys = Float64[]
            colors = Symbol[]

            for row in eachrow(this_df)
                v_full = row[col_full]
                v_step = row[col_step]

                if length(v_full) < 31 || length(v_step) < 31
                    continue
                end

                v_full_cons = v_full[31:end]
                v_step_cons = v_step[31:end]
                n_points = min(length(v_full_cons), length(v_step_cons))

                append!(xs, v_full_cons[1:n_points])
                append!(ys, v_step_cons[1:n_points])
                append!(colors, fill(:red, n_points))
            end

            scatter!(ax, xs, ys;
                color = colors,
                markersize = 4,
                alpha = 1.0
            )

            if !isempty(xs) && !isempty(ys)
                mn = min(minimum(xs), minimum(ys))
                mx = max(maximum(xs), maximum(ys))
                lines!(ax, [mn, mx], [mn, mx]; color=:black, linestyle=:dash)

                if fit_to_1_1_line
                    ss_tot = sum((ys .- mean(ys)).^2)
                    ss_res = sum((ys .- xs).^2)
                    r2_1to1 = 1 - ss_res / ss_tot

                    text!(ax, "R²=$(round(r2_1to1, digits=3))";
                        position=(mx, mn),
                        align=(:right, :bottom),
                        fontsize=10,
                        color=:black
                    )
                else
                    r_val = cor(xs, ys)
                    text!(ax, "r=$(round(r_val, digits=3))";
                        position = (mx, mn),
                        align=(:right, :bottom),
                        fontsize=10,
                        color=:black
                    )
                end
            end
        end

        if save_plot
            filename = string("vector_correlations_", variable, "_", scen, ".png")
            save(filename, fig; px_per_unit=pixels_per_unit)
        end

        display(fig)
    end
end

