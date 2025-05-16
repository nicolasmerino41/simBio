function short_step_correlations(
    df::DataFrame,
    var::Symbol;
    color_by::Symbol = :conn,
    remove_unstable::Bool = false
)
    # ──────────────────────────────────────────────────────────────────────────
    # 1) define your 19 keys & titles
    step_keys = [
        "S1","S2","S3","S4","S5","S6", "S7"
    ]
    step_names = [
        "Full Model", "Global A (Global ϵ)", " Global AE",
        "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
        "Average Biomass"
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

    fig = Figure(; size=(1050, 700))
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
            titlesize = 12, xlabelsize = 10, ylabelsize = 10,
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
        text!(ax, "r=$(round(r_val, digits=2))";
            position = (mx, mn),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end

    display(fig)
end

df = A
df.avg_xi = mean.(A.K_Xi_full)

vect = Vector{Vector}()
vect_consumers = Vector{Vector}()
vect_resources = Vector{Vector}()
for i in 1:nrow(A)
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
    color_by = :epsi
    remove_it = true
    short_step_correlations(df, :rt_press;  color_by = color_by, remove_unstable=remove_it)

    short_step_correlations(df, :after_persistence;  color_by = color_by, remove_unstable=remove_it)

    short_step_correlations(df, :rt_pulse; remove_unstable=remove_it, color_by = color_by)

    short_step_correlations(df, :collectivity; remove_unstable=remove_it, color_by = color_by)

    short_step_correlations(df, :resilience; remove_unstable=remove_it, color_by = color_by)

    short_step_correlations(df, :reactivity; remove_unstable=remove_it, color_by = color_by)
end

res_df = hcat(A.resilience_S1, A.resilience_S2, A.resilience_S3, A.resilience_S4, A.resilience_S5, A.resilience_S6)

A_matrix = A.p_final[1][8]
epsilon_matrix = A.p_final[1][7]
epsilon_matrix = fill(1.0, 50, 50)

psi_full = compute_collectivity(A_matrix, epsilon_matrix)
a1_A, a1_epsilon  = short_transform_for_ladder_step(1, copy(A_matrix), copy(epsilon_matrix))
a2_A, a2_epsilon = short_transform_for_ladder_step(2, copy(A_matrix), copy(epsilon_matrix))
a3_A, a3_epsilon = short_transform_for_ladder_step(3, copy(A_matrix), copy(epsilon_matrix))

psi_a1 = compute_collectivity(a1_A, a1_epsilon)
psi_a2 = compute_collectivity(a2_A, a2_epsilon) 
psi_a3 = compute_collectivity(a3_A, a3_epsilon)
