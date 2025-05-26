function short_step_correlations(
    df::DataFrame,
    var::Symbol;
    color_by::Symbol = :conn,
    remove_unstable::Bool = false
)
    # ──────────────────────────────────────────────────────────────────────────
    # 1) define your 19 keys & titles
    step_keys = [
        "S1","S2","S3","S4","S5","S6", "S7", "S8"]
        # "S7"
    # ]
    step_names = [
        "Full Model", "Global A (Global ϵ)", " Global AE",
        "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
        "Global A (Global ϵ) Mean B", "Global AE Mean B"
    ]

    # ──────────────────────────────────────────────────────────────────────────
    # 1b) optionally filter out any run that went unstable
    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        println("subset size: ", nrow(df))
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 2) set up full and panel‐columns
    full_col   = Symbol(string(var) * "_full")
    panel_cols = Symbol.(string(var) .* "_" .* step_keys)

    full_vals  = df[!, full_col]
    color_vals = df[!, color_by]
    cmin, cmax = minimum(color_vals), maximum(color_vals)

    # ──────────────────────────────────────────────────────────────────────────
    # 3) grid: 19 panels → 5 rows × 4 cols
    ns   = length(panel_cols)
    cols = 3
    rows = ceil(Int, ns/cols)

    fig = Figure(; size=(1000, 570))
    Label(fig[0, 1:cols], uppercase(string(var, " correlations")); fontsize=18)
    Label(fig[0, 2:cols], "Remove unstable: $(remove_unstable)"; fontsize=10)

    # ──────────────────────────────────────────────────────────────────────────
    # 4) loop panels
    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1

        ax = Axis(fig[r, c];
            title     = step_names[idx],
            xlabel    = string(full_col),
            ylabel    = string(panel_cols[idx]),
            titlesize = 20, xlabelsize = 10, ylabelsize = 10,
            xticksize = 15, yticksize = 15,
        )

        ys = df[!, panel_cols[idx]]

        scatter!(ax, full_vals, ys;
            color       = color_vals,
            colormap    = :viridis,
            colorrange  = (cmin, cmax),
            markersize  = 6,
            alpha       = 0.8
        )

        mn = min(minimum(full_vals), minimum(ys))
        mx = max(maximum(full_vals), maximum(ys))
        lines!(ax, [mn, mx], [mn, mx];
            color     = :black,
            linestyle = :dash
        )

        r_val = cor(full_vals, ys)
        text!(ax, "r=$(round(r_val, digits=5))";
            position = (mx, mn),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end

    display(fig)
    return fig
end

df = T
df.avg_xi = mean.(df.K_Xi_full)

vect = Vector{Vector}()
vect_consumers = Vector{Vector}()
vect_resources = Vector{Vector}()
for i in 1:nrow(df)
    g_vecs = fill(df.g_val[i], 30)
    m_vecs = fill(df.m_val[i], 20)
    SL = vcat(g_vecs, m_vecs) .* vcat(df.R_eq[i], df.C_eq[i]) ./ df.K_Xi_full[i]
    SL_cons = m_vecs .* df.C_eq[i] ./ df.K_Xi_full[i][31:50]
    SL_resources = g_vecs .* df.R_eq[i] ./ df.K_Xi_full[i][1:30]

    push!(vect, SL)
    push!(vect_consumers, SL_cons)
    push!(vect_resources, SL_resources)
end
df.SL = mean.(vect)
df.SL_consumers = mean.(vect_consumers)
df.SL_resources = mean.(vect_resources)

begin
    save_plot = false
    color_by = :conn
    remove_it = true
    rt_press = short_step_correlations(df, :rt_press;  color_by = color_by, remove_unstable=remove_it)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/rt_press.png", rt_press)
    end
    
    persistence = short_step_correlations(df, :after_persistence; remove_unstable=remove_it, color_by = color_by)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/persistence.png", persistence)
    end

    rt_pulse = short_step_correlations(df, :rt_pulse; remove_unstable=remove_it, color_by = color_by)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/rt_pulse.png", rt_pulse)
    end

    collectivity = short_step_correlations(df, :collectivity; remove_unstable=remove_it, color_by = color_by)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/collectivity.png", collectivity)
    end

    resilience = short_step_correlations(df, :resilience; remove_unstable=remove_it, color_by = color_by)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/resilience.png", resilience)
    end

    reactivity = short_step_correlations(df, :reactivity; remove_unstable=remove_it, color_by = color_by)
    if save_plot
        save("ThePaper/Ladder/Scenario2/figures/reactivity.png", reactivity)
    end

    rt_med = short_step_correlations(df, :Rmed; color_by = color_by, remove_unstable=remove_it)
end
