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
    # ──────────────────────────────────────────────────────────────────────────
    # 1) define your 19 keys & titles
    step_keys = [
        "S1","S2","S3","S4","S5","S6","S7","S8",
        "S9","S10","S11","S12","S13","S14","S15","S16",
        "S17","S18","S19", "S20"
    ]
    step_names = [
        "Full Model",
        "Full A (Species‐specific ϵ)",    "Full A (Global ϵ)",       "Full A (Re‐randomised ϵ)",
        "Global ±A (ϵ_full)",            "Global ±A (Species‐specific ϵ)", "Global ±A (Global ϵ)",  "Global ±A (Re‐randomised ϵ)",
        "Global A (ϵ_full)",            "Global A (Species‐specific ϵ)", "Global A (Global ϵ)",    "Global A (Re‐randomised ϵ)",
        "Random A (ϵ_full)",             "Random A (Species‐specific ϵ)",  "Random A (Global ϵ)",     "Random A (Re‐randomised ϵ)",
        " Global AE", "Randomize m_cons ↻",            "Randomize ξ̂ ↻",               "Randomize K_res ↻"
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
    cols = 4
    rows = ceil(Int, ns/cols)

    fig = Figure(; size=(900, 700))
    Label(fig[0, 1:cols], uppercase(string(var, " correlations")); fontsize=24)

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

A = deserialize("ThePaper/Ladder/Outputs/Final_results.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_bTermOn.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_fixedM.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_param_10000.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_param.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_100__incredibly_short.jls")
A = deserialize("ThePaper/Ladder/Outputs/Final_results_new_1000_very_short.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_1000_th5.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_10000_th10.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_100000_th15.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_1000_th5_fixedRes_and_rerandom.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_100000_short.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_20000_short.jls")
A = deserialize("ThePaper/Ladder/Outputs/Truth_50000_short_ony_eps1.jls")


df = A
# e.g. colour by connectance
plot_step_correlations(df, :rt_press;  color_by = :rt_press_full, remove_unstable=false)

# or colour by number of surviving species after full press
# plot_step_correlations(df, :before_persistence; color_by = :conn, remove_unstable=false)

# or by scenario
plot_step_correlations(df, :after_persistence;  color_by = :conn, remove_unstable=false)

# or by delta
plot_step_correlations(df, :rt_pulse; remove_unstable=false)

# or by equilibrium abundances
plot_step_correlations(df, :collectivity; remove_unstable=false, color_by = :conn)

# or by resilience
plot_step_correlations(df, :resilience; remove_unstable=false, color_by = :conn)

# or by reactivity
plot_step_correlations(df, :reactivity; remove_unstable=false, color_by = :conn)
