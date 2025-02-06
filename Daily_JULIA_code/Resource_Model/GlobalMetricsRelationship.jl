function measure_cell_sensitivity(all_results_list::Vector{DataFrame}; capped = true, cap_val = 5.0)
    # Create an empty DataFrame to store the cell-level sensitivity metrics.
    cell_sensitivity_df = DataFrame(
         cell_id = Int[],
         avg_sensitivity = Float64[],
         sensitivity_sd = Float64[],
         density = Float64[],
         avg_degree = Float64[],
         avg_clustering = Float64[],
         global_betweenness = Float64[],
         global_closeness = Float64[]
    )
    
    # Loop over each cell (each element of all_results_list)
    for cell_index in 1:length(all_results_list)
         df = all_results_list[cell_index]
         real_cell_index = df[1, :cell]
         # Identify the baseline row: the one where no species is removed.
         baseline_rows = df[df.sp_removed .== "none", :]
         if nrow(baseline_rows) == 0
             @warn "Cell $cell_index has no baseline row; skipping."
             continue
         end
         baseline_effect = baseline_rows[1, :delta_total_biomass]
         
         # Extract the rows corresponding to species removal.
         removal_df = df[df.sp_removed .!= "none", :]
         if nrow(removal_df) == 0
             @warn "Cell $cell_index has no species removal rows; skipping."
             continue
         end
         
         # Define sensitivity as the absolute difference between the removal effect
         # and the baseline effect.
         effects = abs.(removal_df.delta_total_biomass)
         effects = capped ? min.(effects, cap_val) : effects
         avg_sens = mean(effects)
         sens_sd = std(effects)
         
         # Retrieve global network metrics for the cell using your function.
         metrics = compute_food_web_metrics(cell_index; round=false)
         gm = metrics.global_metrics
         
         # Append a new row with this cell's id, sensitivity, and global metrics.
         push!(cell_sensitivity_df, (
            cell_id = real_cell_index,
            avg_sensitivity = avg_sens,
            sensitivity_sd = sens_sd,
            density = gm.density,
            avg_degree = gm.avg_degree,
            avg_clustering = gm.avg_clustering,
            global_betweenness = gm.global_betweenness,
            global_closeness = gm.global_closeness
         ))
    end
    return cell_sensitivity_df
end

# cell_sensitivity_df = measure_cell_sensitivity(all_results_list_even_pi)
# corrected_cell_sensitivity_df = cell_sensitivity_df[cell_sensitivity_df.avg_sensitivity .< 5.0, :]

function map_cell_sensitivity(cell_sensitivity_df::DataFrame; disp = true)
   grid = deepcopy(float(DA_sum))
   for i in axes(grid, 1), j in axes(grid, 2)
        if grid[i, j] == 0
            grid[i, j] = NaN
        elseif grid[i, j] == 1.0
            grid[i, j] = 0.0
        end
    end

   for i in 1:nrow(cell_sensitivity_df)
      
      coord = idx[cell_sensitivity_df[i, :cell_id]]
      local_i, local_j = coord[1], coord[2]

      grid[local_i, local_j] = cell_sensitivity_df[i, :avg_sensitivity] 
    
    end
    if disp
        display(map_plot(grid; palette = custom_palette))
    end
    return grid
end

# grid = map_cell_sensitivity(cell_sensitivity_df; capped = true)

######## PLOTTING THE CORRELATION BETWEEN GLOBAL METRICS AND SENSITIVITY ########
function plot_global_metrics_vs_sensitivity(
    cell_sensitivity_df::DataFrame;
    save = false
)

    # Define the list of global metrics and corresponding labels.
    metrics = [:density, :avg_degree, :avg_clustering, :global_betweenness, :global_closeness]
    metric_labels = ["Density", "Average Degree", "Average Clustering", "Global Betweenness", "Global Closeness"]

    # Create a figure with one subplot per metric.
    n = length(metrics)
    num_cols = 3  # For example, arrange subplots in 3 columns
    num_rows = ceil(Int, n / num_cols)
    fig = Figure(resolution = (1200, 300 * num_rows))
    
    for (i, metric) in enumerate(metrics)
        row = div(i - 1, num_cols) + 1
        col = mod(i - 1, num_cols) + 1
        ax = Axis(fig[row, col],
            title = "$(metric_labels[i]) vs. Average Sensitivity",
            xlabel = metric_labels[i],
            ylabel = "Average Sensitivity",
            limits = (0, maximum(cell_sensitivity_df[!, metric]+cell_sensitivity_df[!, metric]/10), 0, maximum(cell_sensitivity_df[!, :avg_sensitivity]+cell_sensitivity_df[!, :avg_sensitivity]/10))
        )
        scatter!(ax, cell_sensitivity_df[!, metric], cell_sensitivity_df[!, :avg_sensitivity],
                 markersize = 8, color = :blue)
                 
        # Optionally, compute and display the Pearson correlation coefficient.
        # corr_val = cor(corrected_cell_sensitivity_df[!, metric], corrected_cell_sensitivity_df[!, :avg_sensitivity])
        # text!(
        #     ax, "r = $(round(corr_val, digits=2))", position = (0.5, 03),
        #     align = (:left, :top), color = :black, fontsize = 12
        # )
    end

    display(fig)
    for metric in metrics
        max_val = maximum(cell_sensitivity_df[!, metric])
        println("Max $(metric_labels[findfirst(isequal(metric), metrics)]) = $max_val")
    end
    if save 
        save("Plots/global_metrics_vs_sensitivity.png", fig)
    end
end

# save("Plots/global_metrics_vs_sensitivity.png", fig)

##### HERE I WANT TO PLOT THE HERP/PRED RATIO OF OPTIM CONFIGURATION IN THE IBERIAN GRID #####
if false
    begin
    
    # Choose the DataFrame that holds your optimal parametrisation.
    # It should contain a baseline row for each cell with sp_removed == false (i.e. baseline)
    # and have columns "i", "j", and "herb_pred_ratio".
    # Replace Big_P_results with Big_P_even_pi if desired.
    df = deepcopy(Big_P_even_pi)
    # df.sp_removed = parse.(Bool, df.sp_removed)
    baseline_df = filter(row -> !ismissing(row.sp_removed), df)
    baseline_df.sp_removed = map(x -> x == "false" ? false : true, baseline_df.sp_removed)
    baseline_df = filter(row -> row.sp_removed == false, baseline_df)
    
    # Determine the grid dimensions.
    # Here we assume that "i" and "j" indicate the cell coordinates.
    baseline_df.i = parse.(Int, baseline_df.i)
    baseline_df.j = parse.(Int, baseline_df.j)
    # println(baseline_df)
    i_max = maximum(baseline_df.i)
    j_max = maximum(baseline_df.j)
    
    # Create a matrix to hold the herbivore/predator ratio.
    ratio_matrix = fill(NaN, i_max, j_max)
    
    # Fill the matrix using the i, j coordinates.
    for row in eachrow(baseline_df)
        ratio_matrix[row.i, row.j] = row.herb_pred_ratio
    end

    # Now plot the matrix as a heatmap.
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        title = "Herbivore/Predator Ratio (Baseline)",
        xlabel = "Column (j)",
        ylabel = "Row (i)",
        yreversed = true  # So that row 1 appears at the top.
    )
    
    heatmap!(ax, ratio_matrix; colormap = :viridis)
    # Colorbar(fig[1, 2], ax; label = "Herb/Pred Ratio")
    display(fig)
    end

end