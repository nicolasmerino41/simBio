
function stepwise_scatter(df::DataFrame; color_by=:conn, remove_unstable=false)
    step_keys = ["S1", "S3", "S4", "S5", "S6", "S7", "S8"]
    step_names = [
        "Full Model", "Global A (Global ϵ)", " Global AE",
        "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
        "Global A (Global ϵ) Mean B", "Global AE Mean B"
    ]

    # Optionally remove unstable
    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        println("subset size: ", nrow(df))
    end

    color_vals = df[!, color_by]
    cmin, cmax = minimum(color_vals), maximum(color_vals)

    cols = 3
    ns = length(step_keys)
    rows = ceil(Int, ns / cols)
    fig = Figure(; size=(1000, 570))
    Label(fig[0, 1:cols], "rt_press vs tau for each step"; fontsize=18)
    Label(fig[0, 2:cols], "Remove unstable: $(remove_unstable)"; fontsize=10)

    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1

        rt_press_col = Symbol("rt_press_" * step_keys[idx])
        tau_col = Symbol("tau_" * step_keys[idx])

        ax = Axis(fig[r, c];
            title     = step_names[idx],
            xlabel    = string("tau_", step_keys[idx]),
            ylabel    = string("rt_press_", step_keys[idx]),
            titlesize = 16, xlabelsize = 10, ylabelsize = 10,
            xticksize = 15, yticksize = 15,
        )

        xs = df[!, tau_col]
        ys = df[!, rt_press_col]

        scatter!(ax, xs, ys;
            color       = color_vals,
            colormap    = :viridis,
            colorrange  = (cmin, cmax),
            markersize  = 6,
            alpha       = 0.8
        )

        mnx, mxx = minimum(xs), maximum(xs)
        mny, mxy = minimum(ys), maximum(ys)
        # Optional: draw y=x reference if comparable
        # lines!(ax, [mn, mx], [mn, mx]; color=:black, linestyle=:dash)

        # Pearson r-value
        valid = (!isnan).(xs) .& (!isnan).(ys)
        if sum(valid) > 2
            r_val = cor(xs[valid], ys[valid])
            text!(ax, "r=$(round(r_val, digits=5))";
                position = (mxx, mny),
                align    = (:right, :bottom),
                fontsize = 10
            )
        end
    end

    display(fig)
    return fig
end

stepwise_scatter(df; color_by=:conn, remove_unstable=false)
