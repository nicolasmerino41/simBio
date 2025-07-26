function plot_scalar_correlations_glv(
    G::DataFrame;
    scenarios = [:ER, :PL, :MOD],
    metrics = [
        (:resilience, "Resilience"), (:reactivity, "Reactivity"),
        # (:mean_SL, "Mean SL"), 
        (:rt_press, "RT_press"),
        (:rt_pulse, "RT_pulse"),
        (:after_press, "after_press"),
        (:after_pulse, "after_pulse"),
        # (:rmed, "rmed"),
        # (:after_persistence, "Persistence"),
        # (:collectivity, "Collectivity"),
        # (:sigma_over_min_d, "σ/min(d)")
    ],
    fit_to_1_1_line::Bool = true,
    save_plot::Bool = false,
    resolution = (900, 650)
)
    step_names = ["Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS", "Rewiring + ↻C + ↻IS", "Changing groups"]
    steps = collect(1:length(step_names))

    df = G
    fig = Figure(; size=resolution)

    for (i, (sym, label)) in enumerate(metrics)
        for (j, step) in enumerate(steps)
            col_full = Symbol(string(sym, "_full"))
            col_step = Symbol(string(sym, "_S", step))

            x_raw = df[!, col_full]
            y_raw = df[!, col_step]

            x_finite = Float64[]
            y_finite = Float64[]

            for k in 1:min(length(x_raw), length(y_raw))
                xi = x_raw[k]
                yi = y_raw[k]
                if !(ismissing(xi) || ismissing(yi)) && isfinite(xi) && isfinite(yi)
                    push!(x_finite, xi)
                    push!(y_finite, yi)
                end
            end

            # Skip plot if invalid data
            if isempty(x_finite) || isempty(y_finite) || any(!isfinite(v) for v in (minimum(x_finite), maximum(x_finite), minimum(y_finite), maximum(y_finite)))
                @warn "Skipping $sym at step $step: invalid or empty data."
                continue
            end


            mn = min(minimum(x_finite), minimum(y_finite))
            mx = max(maximum(x_finite), maximum(y_finite))

            ax = Axis(
                fig[j, i];
                title = "$label: Full vs $(step_names[step])",
                titlesize = 9,
                xlabelsize = 10,
                ylabelsize = 10,
                xticklabelsize = 10,
                yticklabelsize = 10,
                limits = ((mn, mx), (mn, mx))
            )

            scatter!(ax, x_finite, y_finite; alpha=0.3)

            # 1:1 line
            lines!(ax, [mn, mx], [mn, mx]; color = :black, linestyle = :dash)

            if fit_to_1_1_line
                y_hat = x_finite
                ss_tot = sum((y_finite .- mean(y_finite)).^2)
                ss_res = sum((y_finite .- y_hat).^2)
                r2_1to1 = ss_tot == 0 ? NaN : 1 - ss_res / ss_tot

                if isfinite(r2_1to1) && isfinite(mx) && isfinite(mn)
                    text!(ax, "R²=$(round(r2_1to1, digits=3))";
                        position = (mx, mn),
                        align = (:right, :bottom),
                        fontsize = 10,
                        color = :black
                    )
                end
            else
                r_val = cor(x_finite, y_finite)
                if isfinite(r_val) && isfinite(mx) && isfinite(mn)
                    text!(ax, "r=$(round(r_val, digits=3))";
                        position = (mx, mn),
                        align = (:right, :bottom),
                        fontsize = 10,
                        color = :black
                    )
                end
            end
        end
    end

    if save_plot
        filename = "scalar_correlation_all_scenarios.png"
        save(filename, fig; px_per_unit = 6.0)
    end

    display(fig)
end

function plot_vector_correlations_glv(
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

    step_names = ["Rewiring", "Rewiring + ↻C", "Rewiring + ↻IS", "Rewiring + ↻C + ↻IS", "Changing groups"]

    # color scale
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    for scen in scenarios
        sub = df[df.scen .== scen, :]
        fig = Figure(; size=resolution)
        Label(fig[0, 1:4], uppercase(string(variable, " correlations")); fontsize=12)
        for (i, step) in enumerate((1,2, 3, 4, 5))
            col_full = Symbol(string(variable, suffix[1]))
            col_step = Symbol(string(variable, suffix[2], step))
            ax = Axis(fig[1,i];
                # title  = "$(scen) $(variable): full vs S$step",
                title  = "Full vs $(step_names[step])",
                xlabel = "full", ylabel = "S$step",
                titlesize=10,
                xlabelsize=10, ylabelsize=10,
                xticklabelsize=10, yticklabelsize=10
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