new_results11 = CSV.File("ThePaper/Ladder/Outputs/rppp_guided_results.csv") |> DataFrame
new_results12 = CSV.File("ThePaper/Ladder/Outputs/rppp_guided_results_larger_sample.csv") |> DataFrame
df = df_results

S_vals = [30]

# 1) Grab the column‐names (they are Strings)
rt_cols = filter(n -> startswith(n, "rt_"), names(df))
re_cols = filter(n -> startswith(n, "res_"),   names(df))
@assert length(rt_cols) == length(re_cols)

# 2) Build the “long” DataFrame
n    = nrow(df)
long = DataFrame(step=String[], rt=Float64[], res=Float64[])

for rt_col in rt_cols
    # derive the matching res column
    s    = replace(rt_col, "rt_" => "")          # e.g. "Full", "S1", …
    re_col = "res_$s"
    @assert re_col in names(df)

    df_temp = DataFrame(
      step        = fill(s, n),
      rt = df[!, rt_col],
      res  = df[!, re_col],
    )
    append!(long, df_temp)
end

function get_grid_position(i::Int)
    if i <= 3
        return 1, i
    elseif i <= 7
        return 2, i - 3
    elseif i <= 11
        return 3, i - 7
    else
        return 4, i - 11
    end
end

#############################################################################
#############################################################################
############################# RETURN TIMES ##################################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.C == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.S == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species‑specific ϵ)", "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_rt = df_subset.rt_Full
        Cs = df_subset.C

        fig = Figure(; size = (1400, 900))
        Label(fig[0, :], "RETURN TIMES (S = $S_val)", fontsize=24, tellwidth=false)

        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model RT",
                ylabel = "Step Return Time",
                title  = step_names[i]
            )

            simp_rt = df_subset[!, Symbol("rt_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_rt, simp_rt;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            allvals = vcat(full_rt, simp_rt)
            minv, maxv = minimum(allvals), maximum(allvals)
            lines!(ax, [minv, maxv], [minv, maxv]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_rt, simp_rt), digits=2)
            text!(ax, "r = $r";
                position = (maxv * 0.95, minv + 0.05 * (maxv - minv)),
                align = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)

        display(fig)
    end
end
#############################################################################
#############################################################################
############################ res #####################################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.C == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.S == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species‑specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_res = df_subset.res_Full
        Cs = df_subset.C

        fig = Figure(; size = (1400, 900))
        # Label(fig[0, :], "res (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "res (S = $S_val)", fontsize=24, tellwidth=false)
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model res",
                ylabel = "Step res",
                title  = step_names[i]
            )

            simp_res = df_subset[!, Symbol("res_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_res, simp_res;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            vmin, vmax = minimum(vcat(full_res, simp_res)), maximum(vcat(full_res, simp_res))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_res, simp_res), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end
#############################################################################
#############################################################################
########################### REACTIVITY ######################################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.C == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.S == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species‑specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_rct = df_subset.rea_Full
        Cs = df_subset.C

        fig = Figure(; size = (1400, 900))
        # Label(fig[0, :], "REACTIVITY (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "REACTIVITY (S = $S_val)", fontsize=24, tellwidth=false)

        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Reactivity",
                ylabel = "Step Reactivity",
                title  = step_names[i]
            )

            simp_rct = df_subset[!, Symbol("rea_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_rct, simp_rct;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )
            vmin, vmax = minimum(vcat(full_rct, simp_rct)), maximum(vcat(full_rct, simp_rct))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_rct, simp_rct), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end
#############################################################################
#############################################################################
############################ PERSISTENCE ####################################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.C == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.S == S_val, df)
    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)",
            "Full A (Global ϵ)",
            "Full A (Re‑randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)",
            "Global -Aij&Aij (Species‑specific ϵ)",
            "Global -Aij&Aij (Global ϵ)",
            "Global -Aij&Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)",
            "Global A (Species‑specific ϵ)",
            "Global A (Global ϵ)",
            "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)",
            "Re-randomised A (Species‑specific ϵ)",
            "Re-randomised A (Global ϵ)",
            "Re-randomised A (Re-randomised ϵ)"
        ]

        full_p = df_subset.aft_Full
        Cs = df_subset.C

        fig = Figure(; size = (1400, 900))
        # Label(fig[0, :], "PERSISTENCE (after perturbation) (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "PERSISTENCE (after perturbation) (S = $S_val)", fontsize=24, tellwidth=false)

        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Persistence",
                ylabel = "Step Persistence",
                title  = step_names[i]
            )

            simp_p = df_subset[!, Symbol("aft_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_p, simp_p;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            vmin, vmax = minimum(vcat(full_p, simp_p)), maximum(vcat(full_p, simp_p))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_p, simp_p), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align    = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end
#############################################################################
#############################################################################
######################## SENSITIVITY CORRELATION ############################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.C == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of R and C
    df_subset = filter(row -> row.S == S_val, df)

    begin
        step_keys  = ["Full", "S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full Model", "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        Cs = df_subset.C

        # Label(fig[0, :], "SENSITIVITY CORRELATION (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "SENSITIVITY CORRELATION (S = $S_val)", fontsize=24, tellwidth=false)
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Sensitivity",
                ylabel = "Step Sensitivity",
                title  = step_names[i]
            )

            full_s = df_subset.scorr_Full
            simp_s = df_subset[!, Symbol("scorr_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_s, simp_s;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            vmin, vmax = minimum(vcat(full_s, simp_s)), maximum(vcat(full_s, simp_s))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_s, simp_s), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align    = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end
#############################################################################
#############################################################################
########################### MAXIMUM OVERSHOOT ###############################
#############################################################################   
#############################################################################
for S_val in S_vals
    # Filter the DataFrame for the specific combination of R and C
    df_subset = filter(row -> row.S == S_val, df)

    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        Cs = df_subset.C

        # Label(fig[0, :], "MAXIMUM OVERSHOOT (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "MAXIMUM OVERSHOOT (S = $S_val)", fontsize=24, tellwidth=false)
        
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Maximum Overshoot",
                ylabel = "Step Maximum Overshoot",
                title  = step_names[i]
            )
            full_mo = df_subset.os_Full
            simp_mo = df_subset[!, Symbol("os_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_mo, simp_mo;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            vmin, vmax = minimum(vcat(full_mo, simp_mo)), maximum(vcat(full_mo, simp_mo))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)

            r = round(cor(full_mo, simp_mo), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align    = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end

#############################################################################
#############################################################################
######################## INTEGRATED RECOVERY ERROR ##########################
#############################################################################
#############################################################################
for S_val in S_vals
    # Filter the DataFrame for the specific combination of R and C
    df_subset = filter(row -> row.S == S_val, df)

    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        Cs = df_subset.C
        
        Label(fig[0, :], "INTEGRATED RECOVERY ERROR (S = $S_val)", fontsize=24, tellwidth=false)
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Integrated Recovery Error",
                ylabel = "Step Integrated Recovery Error",
                title  = step_names[i]
            )
            full_ire = df_subset.ire_Full
            simp_ire = df_subset[!, Symbol("ire_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(Cs)))
            color_indices = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_ire, simp_ire;
                markersize = 6,
                color = Cs,
                colormap = :viridis,
                colorrange = (minimum(Cs), maximum(Cs))
            )

            vmin, vmax = minimum(vcat(full_ire, simp_ire)), maximum(vcat(full_ire, simp_ire))
            lines!(ax, [vmin, vmax], [vmin, vmax]; color=:black, linestyle=:dash, linewidth=1)
            r = round(cor(full_ire, simp_ire), digits=2)
            text!(ax, "r = $r";
                position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                align    = (:right, :bottom),
                fontsize = 12
            )
        end
        Colorbar(fig[1, 5], limits = (minimum(Cs), maximum(Cs)), colormap = :viridis)
        display(fig)
    end
end