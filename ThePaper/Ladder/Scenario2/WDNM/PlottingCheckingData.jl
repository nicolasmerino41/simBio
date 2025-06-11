R = deserialize("ThePaper/Ladder/Outputs/checking/checking20000.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography60000.jls")

desired = [
  :conn, :IS, :scen, :delta, :epsi, :m_val, :g_val, :ite,
  :resilience_full, :reactivity_full, :tau_full, :mean_tau_full, :sigma_over_min_d_full, :SL_full, :mean_SL_full, :inverse_tau_full, :mean_inverse_tau_full,
  :resilience_S1, :reactivity_S1, :tau_S1, :mean_tau_S1, :sigma_over_min_d_S1, :SL_S1, :mean_SL_S1, :inverse_tau_S1, :mean_inverse_tau_S1,
  :resilience_S2, :reactivity_S2, :tau_S2, :mean_tau_S2, :sigma_over_min_d_S2, :SL_S2, :mean_SL_S2, :inverse_tau_S2, :mean_inverse_tau_S2
]

G = R[!, desired]

step_keys = ["_full","_S1","_S2"]
res_cols = Symbol.("resilience" .* step_keys)
G = filter(row -> all(row[c] < 0 for c in res_cols), G)
println("subset size: ", nrow(G))
min_d_cols = Symbol.("sigma_over_min_d" .* step_keys)
G = filter(row -> all(row[c] < 500.0 for c in min_d_cols), G)

# plotting
scenarios = [:ER]
metrics = [(:resilience, "Resilience"), (:reactivity,"Reactivity"), (:mean_tau,"Mean SL"), (:sigma_over_min_d,"σ/min(d)")]
for scen in scenarios
    df = G[G.scen .== scen, :]
    fig = Figure(; size=(1000,450))
    for (i,(sym,label)) in enumerate(metrics)
        for (j,step) in enumerate((1,2))  # comparing step0 to step1 and step2
            ax = Axis(fig[j, i];
                title = "$(scen) $label: full vs rewiring$(step)",
                xlabel="orig", ylabel="step$(step)")
            x = df[!, Symbol(string(sym,"_full"))]
            y = df[!, Symbol(string(sym,"_S",step))]
            scatter!(ax, x, y; alpha=0.3)
            
            mn = min(minimum(x), minimum(y))
            mx = max(maximum(x), maximum(y))

            lines!(
                ax, [mn, mx], [mn, mx];
                color     = :black,
                linestyle = :dash
            )

            r_val = cor(x, y)
            text!(ax, "r=$(round(r_val, digits=5))";
                position = (mx, mn),
                align    = (:right, :bottom),
                fontsize = 10
            )
        end
    end
    display(fig)
end

################### for tau vector correlations ###################
function plot_vector_correlations(df::DataFrame;
                                  scenarios=[:ER], # [:ER,:PL,:MOD],
                                  color_by=:conn,
                                  variable::Symbol = :tau)
    # mapping from variable name to the full‐and‐step column suffix
    suffix = Dict(
        :tau     => "_full" => "_S",     # tau_full, tau_S1, tau_S2
        :inverse_tau => "_full" => "_S",
        :SL      => "_full" => "_S"
    )[variable]

    for scen in scenarios
        sub = df[df.scen .== scen, :]
        fig = Figure(; size=(900,450))

        for (i, step) in enumerate((1,2))
            col_full = Symbol(string(variable, suffix[1]))
            col_step = Symbol(string(variable, suffix[2], step))
            ax = Axis(fig[1,i];
                title  = "$(scen) $(variable): full vs S$step",
                xlabel = "full", ylabel = "S$step"
            )
            xs=[]; ys=[]; cs=[]
            for row in eachrow(sub)
                v_full = row[col_full]
                v_step = row[col_step]
                append!(xs, v_full)
                append!(ys, v_step)
                append!(cs, fill(row[color_by], length(v_full)))
            end

            scatter!(ax, xs, ys;
                colormap   = :viridis,
                # color      = cs,
                colorrange = (minimum(cs), maximum(cs)),
                markersize = 4,
                alpha      = 0.6
            )

            mn = min(minimum(xs), minimum(ys))
            mx = max(maximum(xs), maximum(ys))
            lines!(ax, [mn,mx], [mn,mx]; color=:black, linestyle=:dash)

            r = cor(xs, ys)
            text!(ax, "r=$(round(r,digits=3))",
                position = (mx, mn),
                align    = (:right, :bottom),
                fontsize = 10
            )
        end

        # Colorbar(fig[1,3], colormap=:viridis,
        # label = string(color_by),
        # limits=(minimum(cs), maximum(cs))
        # )

        display(fig)
    end
end

# Example calls:
plot_vector_correlations(G; variable=:tau)
plot_vector_correlations(G; variable=:inverse_tau)
plot_vector_correlations(G; variable=:SL)
