################ INFORMATION #######################
# a) Ladder1.jl is for all the functions
# b) Ladder2.jl is for running the simulations
# c) Ladder3.jl is for post-processing and plotting
#############################################################################
#############################################################################
################ Post-processing Analysis and Plotting ######################
#############################################################################
#############################################################################
using GLM
# (A) Scatter plot of resilience (predictor) vs. full-model mean return time.
begin
    fig1 = Figure(; size=(1200,600))
    ax1 = Axis(fig1[1,1], xlabel="Resilience", ylabel="Mean Return Time (Full Model)",
            title="Return Time vs. Resilience")
    # df_filtered = filter(row -> row.resilience > -2.0 && row.resilience < 0.0, df_results)
    df_filtered = filter(row -> !iszero(row.return_time_full), df_results2)
    df_filtered = filter(row -> row.fully_stable, df_filtered)
    MK.scatter!(ax1, df_filtered.resilience, df_filtered.return_time_full, markersize=8, color=:blue)
    
    # Fit a linear model for return_time_full ~ resilience.
    model1 = lm(@formula(return_time_full ~ resilience), df_filtered)
    res_lin = sort(df_filtered.resilience)
    predicted = predict(model1)
    # For plotting a regression line, we sort the data.
    sorted_idx = sortperm(df_filtered.resilience)
    MK.lines!(ax1, df_filtered.resilience[sorted_idx], predicted[sorted_idx], color=:red, linewidth=2)
    # Label(fig1, "Pearson cor = $(round(cor(df_filtered.resilience, df_filtered.return_time_full), digits=2))", position=(20,20))
    display(fig1)
end

begin
    # (B) Scatter plot comparing full model vs. simplified model return times.
    fig2 = Figure(; size=(1200,600))
    ax2 = Axis(fig2[1,1], xlabel="Return Time (Full Model)", ylabel="Return Time (Simplified Model)",
            title="Full vs. Simplified Return Times")
    MK.scatter!(ax2, df_filtered.return_time_full, df_filtered.return_time_simplified, markersize=8, color=:green)
    # Plot a 1:1 line.
    min_val = minimum(vcat(df_filtered.return_time_full, df_filtered.return_time_simplified))
    max_val = maximum(vcat(df_filtered.return_time_full, df_filtered.return_time_simplified))
    lines!(ax2, [min_val, max_val], [min_val, max_val], color=:black, linestyle=:dash, linewidth=2)
    display(fig2)
end

begin
    # (C) Analyze persistence and relative variance vs. connectance.
    fig3 = Figure(; size=(1500,1000))
    ax3 = Axis(fig3[1,1], xlabel="Connectance", ylabel="Persistence (Initial)",
            title="Initial Persistence vs. Connectance")
    ax4 = Axis(fig3[1,2], xlabel="Connectance", ylabel="Persistence Full(Post-Perturbation)",
            title="Post-Perturbation Persistence Full vs. Connectance")
    ax5 = Axis(fig3[2, 1], xlabel="Connectance", ylabel="Persistence Simplified(Post-Perturbation)",
            title="Post-Perturbation Persistence Simplified vs. Connectance")
    ax6 = Axis(fig3[2, 2], xlabel="Persistence Full(Post-Perturbation)", ylabel="Persistence Simplified(Post-Perturbation)")
    
    MK.scatter!(ax3, df_filtered.connectance, df_filtered.persistence, markersize=8, color=:purple)
    MK.scatter!(ax4, df_filtered.connectance, df_filtered.after_persistence_full, markersize=8, color=:orange)
    MK.scatter!(ax5, df_filtered.connectance, df_filtered.after_persistence_simpl, markersize=8, color=:brown)
    MK.scatter!(ax6, df_filtered.after_persistence_full, df_filtered.after_persistence_simpl, markersize=8, color=:red)

    display(fig3)
end

# Print simple correlation statistics:
println("\nCorrelation between resilience and full-model return time: ",
        round(cor(df_filtered.resilience, df_filtered.return_time_full), digits=3))
println("Correlation between full and simplified return times: ",
        round(cor(df_filtered.return_time_full, df_filtered.return_time_simplified), digits=3))
println("Correlation between relative variance and persistence: ",
        round(cor(df_results.relVar, df_results.persistence), digits=3))