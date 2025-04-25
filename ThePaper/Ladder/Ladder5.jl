# new_results4 is 10,20,30,40 species, Normal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean*0.1
new_results4 = deserialize("ThePaper/Ladder/Outputs/new_results4_total_species_fixed.jls")

# new_results5 is 10,20,30,40 species, LogNormal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean
new_results5 = deserialize("ThePaper/Ladder/Outputs/new_results5.jls")

# new_results6 is 10,20,30,40 species, LogNormal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean*0.1
# and includes overshoot and ire
# new_results6 = deserialize("ThePaper/Ladder/Outputs/new_results6.jls")

# new_results6.2 is AGAIN 10,20,30,40 species, LogNormal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean
new_results6 = deserialize("ThePaper/Ladder/Outputs/new_results6.2.jls")

# new_results7 is 20,40,60,80 species, LogNormal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean
new_results7 = deserialize("ThePaper/Ladder/Outputs/new_results7.jls")

# new_results8 is 25,50,75,100 species, LogNormal Abundance Distribution & epsilon = 0.2 with sd = epsilon_mean
new_results8 = deserialize("ThePaper/Ladder/Outputs/new_results8.jls")

# new_results9 is 75 species, LogNormal Abundance Distribution & epsilon = 0.1 with sd = epsilon_mean, d_value = 2.0, pred_mortality = 0.05
new_results9 = deserialize("ThePaper/Ladder/Outputs/new_results9.jls")
new_results10 = deserialize("ThePaper/Ladder/Outputs/new_results10.jls")
new_results11 = deserialize("ThePaper/Ladder/Outputs/new_results11.jls")
df = new_results4
if df == new_results4 || df == new_results5 || df == new_results6 
    S_vals = [10, 20, 30, 40]
elseif df == new_results7
    S_vals = [20, 40, 60, 80]
elseif df == new_results8
    S_vals = [25, 50, 75, 100]
end

# 1) Grab the column‐names (they are Strings)
rt_cols = filter(n -> startswith(n, "return_time_"), names(df))
re_cols = filter(n -> startswith(n, "resilience_"),   names(df))
@assert length(rt_cols) == length(re_cols)

# 2) Build the “long” DataFrame
n    = nrow(df)
long = DataFrame(step=String[], return_time=Float64[], resilience=Float64[])

for rt_col in rt_cols
    # derive the matching resilience column
    s    = replace(rt_col, "return_time_" => "")          # e.g. "Full", "S1", …
    re_col = "resilience_$s"
    @assert re_col in names(df)

    df_temp = DataFrame(
      step        = fill(s, n),
      return_time = df[!, rt_col],
      resilience  = df[!, re_col],
    )
    append!(long, df_temp)
end

println("long has $(nrow(long)) rows, cols: ", names(long))

if false
    # Prepare figure and axis
    fig = Figure(; size = (1000, 600))
    ax  = Axis(fig[1, 1];
        xlabel = "Resilience",
        ylabel = "Return Time",
        title  = "Return Time vs Resilience (Full + Simplified)",
        # xscale = log10,
        # yscale = log10
    )

    # Generate one distinct color per step
    steps  = unique(long.step)
    nsteps = length(steps)
    palette = get(ColorSchemes.tab10, range(0, 1, length = nsteps))

    # Plot each step
    for (i, s) in enumerate(steps)
        sel = long.step .== s
        scatter!(
            ax,
            long.resilience[sel],     # invert if you really want negative on the x-axis
            long.return_time[sel],
            color      = palette[i],
            markersize = 6,
            label      = s,
        )
    end

    # Add legend in the bottom‑right
    axislegend(ax; position = :lt, tellwidth = false)

    display(fig)
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
#     df_subset = filter(row -> row.resource_count == r_val && row.consumer_count == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.species_count == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species‑specific ϵ)", "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_rt = df_subset.return_time_Full
        consumer_counts = df_subset.consumer_count

        fig = Figure(; size = (1400, 900))
        Label(fig[0, :], "RETURN TIMES (S = $S_val)", fontsize=24, tellwidth=false)

        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model RT",
                ylabel = "Step Return Time",
                title  = step_names[i]
            )

            simp_rt = df_subset[!, Symbol("return_time_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_rt, simp_rt;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)

        display(fig)
    end
end
#############################################################################
#############################################################################
############################ RESILIENCE #####################################
#############################################################################
#############################################################################
# for r_val in R_vals, c_val in C_vals
#     # Filter the DataFrame for the specific combination of R and C
#     df_subset = filter(row -> row.resource_count == r_val && row.consumer_count == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.species_count == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species‑specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_res = df_subset.resilience_Full
        consumer_counts = df_subset.consumer_count

        fig = Figure(; size = (1400, 900))
        # Label(fig[0, :], "RESILIENCE (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "RESILIENCE (S = $S_val)", fontsize=24, tellwidth=false)
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Resilience",
                ylabel = "Step Resilience",
                title  = step_names[i]
            )

            simp_res = df_subset[!, Symbol("resilience_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_res, simp_res;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
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
#     df_subset = filter(row -> row.resource_count == r_val && row.consumer_count == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.species_count == S_val, df)
    begin
        step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species‑specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re‑randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        full_rct = df_subset.reactivity_Full
        consumer_counts = df_subset.consumer_count

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

            simp_rct = df_subset[!, Symbol("reactivity_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_rct, simp_rct;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
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
#     df_subset = filter(row -> row.resource_count == r_val && row.consumer_count == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of S
    df_subset = filter(row -> row.species_count == S_val, df)
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

        full_p = df_subset.after_persistence_Full
        consumer_counts = df_subset.consumer_count

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

            simp_p = df_subset[!, Symbol("after_persistence_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_p, simp_p;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
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
#     df_subset = filter(row -> row.resource_count == r_val && row.consumer_count == c_val, df)
for S_val in S_vals
    # Filter the DataFrame for the specific combination of R and C
    df_subset = filter(row -> row.species_count == S_val, df)

    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        consumer_counts = df_subset.consumer_count

        # Label(fig[0, :], "SENSITIVITY CORRELATION (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "SENSITIVITY CORRELATION (S = $S_val)", fontsize=24, tellwidth=false)
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Sensitivity",
                ylabel = "Step Sensitivity",
                title  = step_names[i]
            )

            full_s = df_subset.sensitivity_corr_Full
            simp_s = df_subset[!, Symbol("sensitivity_corr_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_s, simp_s;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
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
    df_subset = filter(row -> row.species_count == S_val, df)

    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        consumer_counts = df_subset.consumer_count

        # Label(fig[0, :], "MAXIMUM OVERSHOOT (R = $r_val, C = $c_val)", fontsize=24, tellwidth=false)
        Label(fig[0, :], "MAXIMUM OVERSHOOT (S = $S_val)", fontsize=24, tellwidth=false)
        
        for (i, s) in enumerate(step_keys)
            row, col = get_grid_position(i)

            ax = Axis(fig[row, col];
                xlabel = "Full Model Maximum Overshoot",
                ylabel = "Step Maximum Overshoot",
                title  = step_names[i]
            )
            full_mo = df_subset.overshoot_Full
            simp_mo = df_subset[!, Symbol("overshoot_$s")]

            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_mo, simp_mo;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
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
    df_subset = filter(row -> row.species_count == S_val, df)

    begin
        step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11", "S12", "S13", "S14", "S15"]
        step_names = [
            "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
            "Global -Aij&Aij (ϵ_full)", "Global -Aij&Aij (Species-specific ϵ)", "Global -Aij&Aij (Global ϵ)", "Global -Aij&Aij (Re-randomised ϵ)",
            "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
            "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
        ]

        fig = Figure(; size = (1400, 900))
        consumer_counts = df_subset.consumer_count
        
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
            cmap = cgrad(:viridis, length(unique(consumer_counts)))
            color_indices = [findfirst(==(val), sort(unique(consumer_counts))) for val in consumer_counts]
            color_vals = cmap[color_indices]
            
            scatter!(
                ax, full_ire, simp_ire;
                markersize = 6,
                color = consumer_counts,
                colormap = :viridis,
                colorrange = (minimum(consumer_counts), maximum(consumer_counts))
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
        Colorbar(fig[1, 5], limits = (minimum(consumer_counts), maximum(consumer_counts)), colormap = :viridis)
        display(fig)
    end
end