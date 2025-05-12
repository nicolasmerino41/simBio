using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")

# ────────────────────────────────────────────────────────────────────────────────
function plot_step_correlations(
    df::DataFrame,
    var::Symbol;            
    color_by::Symbol = :conn, 
    remove_unstable::Bool = false
)
    # 1) your keys & labels
    step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8",
                  "S9","S10","S11","S12","S13","S14","S15", "S16"]
    step_names = [
        "Full Model",
        "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
        "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species-specific ϵ)",
        "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re-randomised ϵ)",
        "Global A (ϵ_full)", "Global A (Species-specific ϵ)",
        "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
        "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)",
        "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
    ]

    # Optionally drop any runs that went unstable at any step
    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
    end

    # The x-axis is always the original “_full” column;
    # the y-axes are the S1–S16 columns.
    full_col   = Symbol(string(var), "_full")
    panel_cols = Symbol.(string(var) .* "_" .* step_keys)

    full_vals  = df[!, full_col]
    color_vals = df[!, color_by]
    cmin, cmax = minimum(color_vals), maximum(color_vals)

    # Layout: 16 panels → 4×4
    ns   = length(panel_cols)
    rows = 4
    cols = 4

    fig = Figure(; size=(900, 580))
    Label(fig[0, 1:cols], uppercase(string(var, " correlations")); fontsize=24)

    # 5) loop panels
    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1

        ax = Axis(fig[r, c];
            xlabel    = string(full_col),
            ylabel    = string(panel_cols[idx]),
            title     = step_names[idx],
            titlesize = 12,
            xlabelsize= 10,
            ylabelsize= 10
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
        lines!(ax, [mn,mx], [mn,mx];
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

A = deserialize("ThePaper/Ladder/Outputs/Final_results.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_bTermOn.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_fixedM.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_param_10000.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_param.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_param_short.jls")
subs = filter(row -> row.stable_S16 == 1, A)

# e.g. colour by connectance
plot_step_correlations(A, :rt_press;  color_by = :conn, remove_unstable=false)

# or colour by number of surviving species after full press
plot_step_correlations(A, :before_persistence; color_by = :conn, remove_unstable=false)

# or by scenario
plot_step_correlations(A, :after_persistence;  color_by = :conn, remove_unstable=false)

# or by delta
plot_step_correlations(A, :rt_pulse; remove_unstable=false)

# or by equilibrium abundances
plot_step_correlations(A, :collectivity; remove_unstable=false, color_by = :IS)

# or by resilience
plot_step_correlations(A, :resilience; remove_unstable=false, color_by = :conn)

# or by reactivity
plot_step_correlations(A, :reactivity; remove_unstable=false, color_by = :conn)
