#############################################################################
#############################################################################
########################### RETURN TIMES ####################################
#############################################################################
#############################################################################
function plot_return_times_by_param(df, S_vals; color_by = :skew)
    step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8",
                  "S9","S10","S11","S12","S13","S14","S15"]
    step_names = [
        "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
        "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species-specific ϵ)",
        "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re-randomised ϵ)",
        "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)",
        "Global A (Re-randomised ϵ)", "Re-randomised A (ϵ_full)",
        "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)",
        "Re-randomised A (Re-randomised ϵ)"
    ]

    for S_val in S_vals
        for α in sort(unique(df.param))
            # subset just PL scenario with this S and α
            df_sub = filter(r -> r.S == S_val && r.scenario == :PL && r.param == α, df)
            full_rt = df_sub.rt_Full
            Cs      = df_sub[!, color_by]

            fig = Figure(; size = (1400, 900))
            Label(fig[1, 1:3], 
                  "RETURN TIME correlations for S=$S_val, α=$α",
                  fontsize = 24)

            # plot each simplification step
            for (i, key) in enumerate(step_keys)
                # compute grid position 4×3
                row = 1 + fld(i-1, 3)  # 1..4
                col = 1 + mod(i-1, 3)  # 1..3

                ax = Axis(fig[row+1, col];
                          xlabel = "Full Model RT",
                          ylabel = "Step RT",
                          title  = step_names[i])

                simp_rt = df_sub[!, Symbol("rt_$key")]

                cmap    = cgrad(:viridis, length(unique(Cs)))
                color_indices    = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
                color_vals       = cmap[color_indices]

                scatter!(
                    ax, full_rt, simp_rt;
                    markersize = 6,
                    color = Cs,
                    colormap = :viridis,
                    colorrange = (minimum(Cs), maximum(Cs))    
                )

                # unity line
                vmin, vmax = minimum(vcat(full_rt, simp_rt)), maximum(vcat(full_rt, simp_rt))
                lines!(ax, [vmin, vmax], [vmin, vmax];
                       color = :black, linestyle = :dash, linewidth = 1)

                # annotation
                r = round(cor(full_rt, simp_rt), digits = 2)
                text!(ax, "r = $r";
                      position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                      align    = (:right, :bottom),
                      fontsize = 12)
            end

            # colorbar alongside the 4×3 grid
            Colorbar(fig[2, 4], 
                     colormap = :viridis,
                     limits   = (minimum(Cs), maximum(Cs)),
                     label    = String(color_by))

            display(fig)
        end
    end
end

# then simply call
plot_return_times_by_param(df, [25, 50, 75, 100]; color_by = :skew)
plot_return_times_by_param(df, [25, 50, 75, 100]; color_by = :C)
plot_return_times_by_param(df, [25, 50, 75, 100]; color_by = :param)

#############################################################################
#############################################################################
########################### REACTIVITY ######################################
#############################################################################
#############################################################################
function plot_reactivity_by_param(df, S_vals; color_by = :skew)
    step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8",
                  "S9","S10","S11","S12","S13","S14","S15"]
    step_names = [
        "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
        "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species-specific ϵ)",
        "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re-randomised ϵ)",
        "Global A (ϵ_full)", "Global A (Species-specific ϵ)", "Global A (Global ϵ)",
        "Global A (Re-randomised ϵ)", "Re-randomised A (ϵ_full)",
        "Re-randomised A (Species-specific ϵ)", "Re-randomised A (Global ϵ)",
        "Re-randomised A (Re-randomised ϵ)"
    ]

    for S_val in S_vals
        for α in sort(unique(df.param))
            # subset just PL scenario with this S and α
            df_sub = filter(r -> r.S == S_val && r.scenario == :PL && r.param == α, df)
            full_rea = df_sub.rea_Full
            Cs      = df_sub[!, color_by]

            fig = Figure(; size = (1400, 900))
            Label(fig[1, 1:3], 
                  "REACTIVITY correlations for S=$S_val, α=$α",
                  fontsize = 24)

            # plot each simplification step
            for (i, key) in enumerate(step_keys)
                # compute grid position 4×3
                row = 1 + fld(i-1, 3)  # 1..4
                col = 1 + mod(i-1, 3)  # 1..3

                ax = Axis(fig[row+1, col];
                          xlabel = "Full Model RT",
                          ylabel = "Step RT",
                          title  = step_names[i])

                simp_rea = df_sub[!, Symbol("rea_$key")]

                cmap    = cgrad(:viridis, length(unique(Cs)))
                color_indices    = [findfirst(==(val), sort(unique(Cs))) for val in Cs]
                color_vals       = cmap[color_indices]

                scatter!(
                    ax, full_rea, simp_rea;
                    markersize = 6,
                    color = Cs,
                    colormap = :viridis,
                    colorrange = (minimum(Cs), maximum(Cs))    
                )

                # unity line
                vmin, vmax = minimum(vcat(full_rea, simp_rea)), maximum(vcat(full_rea, simp_rea))
                lines!(ax, [vmin, vmax], [vmin, vmax];
                       color = :black, linestyle = :dash, linewidth = 1)

                # annotation
                r = round(cor(full_rea, simp_rea), digits = 2)
                text!(ax, "r = $r";
                      position = (vmax * 0.95, vmin + 0.05 * (vmax - vmin)),
                      align    = (:right, :bottom),
                      fontsize = 12)
            end

            # colorbar alongside the 4×3 grid
            Colorbar(fig[2, 4], 
                     colormap = :viridis,
                     limits   = (minimum(Cs), maximum(Cs)),
                     label    = String(color_by))

            display(fig)
        end
    end
end

plot_reactivity_by_param(df, [25, 50, 75, 100]; color_by = :skew)
plot_reactivity_by_param(df, [25, 50, 75, 100]; color_by = :C)
plot_reactivity_by_param(df, [25, 50, 75, 100]; color_by = :param)