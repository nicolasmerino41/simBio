# ───────────────────────────────────────────────────────────────────────────────
# This must be one contiguous function definition; nothing above or below it
# should be at the same indentation level.
# ───────────────────────────────────────────────────────────────────────────────
function step_correlations_END(df::DataFrame, var::Symbol; color_by::Symbol = :conn, remove_unstable = false)
    # We only have steps 2 and 3 (step1 is the full model itself)
    step_keys   = ["S1", "S2", "S3"]
    step_titles = ["Full Model", "Global A (Global ϵ)", "Global AE"]

    # ──────────────────────────────────────────────────────────────────────────
    # 1b) optionally filter out any run that went unstable
    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        println("subset size: ", nrow(df))
    end

    # 1) determine which columns to plot
    full_col   = Symbol(string(var, "_full"))           # e.g. :resilience_full
    panel_cols = Symbol.(string.(var, "_", step_keys))

    # 2) pull out the arrays from the DataFrame
    full_vals  = df[!, full_col]
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    # 3) build a 1×2 grid
    fig = Figure(size=(900, 400))
    Label(fig[0, 1:2], uppercase(string(var, " correlations")), fontsize=18)

    for idx in 1:3
        ax = Axis(fig[1, idx];
            title     = step_titles[idx],
            # xlabel    = full_col,
            # ylabel    = panel_cols[idx],
            # titlesize = 16, labelsize = 12, ticklabelsize = 10
        )

        ys = df[!, panel_cols[idx]]
        scatter!(ax, full_vals, ys;
            colormap   = :viridis,
            color      = color_vals,
            colorrange = (cmin, cmax),
            markersize = 6,
            alpha      = 0.8
        )

        # unity line
        mn = min(minimum(full_vals), minimum(ys))
        mx = max(maximum(full_vals), maximum(ys))
        lines!(ax, [mn, mx], [mn, mx]; color = :black, linestyle = :dash)

        # Pearson r
        r_val = cor(full_vals, ys)
        text!(ax, "r=$(round(r_val, digits=5))";
            position = (mx, mn),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end

    display(fig)
end
# ───────────────────────────────────────────────────────────────────────────────

A = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder.jls")
# And when you call it, do so *after* this file has been included:
begin 
    color_by = :conn
    remove_it = false
    step_correlations_END(A, :resilience; color_by = color_by, remove_unstable = remove_it)
    step_correlations_END(A, :reactivity; color_by = color_by, remove_unstable = remove_it)
    step_correlations_END(A, :Rmed; color_by = color_by, remove_unstable = remove_it)
    step_correlations_END(A, :rt; color_by = color_by, remove_unstable = remove_it)
    step_correlations_END(A, :before; color_by = color_by, remove_unstable = remove_it)
    step_correlations_END(A, :after; color_by = color_by, remove_unstable = remove_it)
end
