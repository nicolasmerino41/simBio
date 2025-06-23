function short_step_correlations_S1S2(
    df::DataFrame;
    color_by::Symbol = :conn_bi,
    remove_unstable::Bool = false
)
    # which metrics to do
    metrics    = (:resilience, :reactivity, :persistence)
    steps      = ("S1","S2")
    step_titles= ("Within-group reshuffle","Between-group reshuffle")
    panel_cols = Symbol.([string(m,"_",s) for m in metrics for s in steps])
    full_cols  = Symbol.(metrics)   # original reality
    color_vals = df[!, color_by]
    cmin, cmax = minimum(color_vals), maximum(color_vals)

    fig = Figure(; size = (900,600))
    for (i, met) in enumerate(metrics), (j, step) in enumerate(steps)
        ax = Axis(fig[i, j],
            title    = step_titles[j],
            xlabel   = string(met),
            ylabel   = string(met, "_", step),
            titlesize= 12)

        xs = df[!, full_cols[i]]
        ys = df[!, Symbol(string(full_cols[i],"_",step))]

        scatter!(ax, xs, ys;
            color       = color_vals,
            colormap    = :viridis,
            colorrange  = (cmin, cmax),
            markersize  = 6,
            alpha       = 0.8)

        # 45° line
        mn, mx = min(minimum(xs), minimum(ys)), max(maximum(xs), maximum(ys))
        lines!(
            ax, [mn,mx], [mn,mx];
            color=:black, linestyle=:dash
        )

        # Pearson r
        # r_val = cor(xs, ys)
        # text!(ax, "r=$(round(r_val, digits=5))";
        #     position = (mx, mn),
        #     align    = (:right, :bottom),
        #     fontsize = 10
        # )
        # R² to the 1:1 line (y_hat = x)
        y_hat = xs
        ss_tot = sum((ys .- mean(ys)).^2)
        ss_res = sum((ys .- y_hat).^2)
        r2_1to1 = 1 - ss_res / ss_tot

        text!(ax, "R²=$(round(r2_1to1, digits=5))";
            position=(mx, mn),
            align=(:right, :bottom),
            fontsize=10,
            color=:black
        )
    end

    # Colorbar(fig[1,3], fig, label = string(color_by))  # put colorbar in the top-right
    display(fig)
end

short_step_correlations_S1S2(df3)