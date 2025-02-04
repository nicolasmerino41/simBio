# ----------------------------
# 1. Compute CVI and NRI for each cell
# ----------------------------
function new_measure_cell_stability(all_results_list::Vector{DataFrame})
    # Create a DataFrame to store cell-level stability metrics.
    cell_stability_df = DataFrame(
         cell_id = Int[],
         baseline = Float64[],
         ARI = Float64[],
         RV = Float64[],
         CVI = Float64[],
         density = Float64[],
         avg_clustering = Float64[],
         global_betweenness = Float64[],
         global_closeness = Float64[],
         NRI = Float64[]
    )

    max_density = 0.0
    max_avg_clust = 0.0
    max_betweenness = 0.0
    # Find maximum NRI
    for i in 1:length(idx)
        metrics = compute_food_web_metrics(i; round=false)
        gm = metrics.global_metrics
        max_density = max(gm.density, max_density)
        max_avg_clust = max(gm.avg_clustering, max_avg_clust)
        max_betweenness = max(gm.global_betweenness, max_betweenness)
    end
    
    for i in 1:length(all_results_list)
        df = all_results_list[i]
        cell_index = df[1, :cell]
        # Identify the baseline row: sp_removed == "none"
        baseline_rows = df[df.sp_removed .== "none", :]
        if nrow(baseline_rows) == 0
            @warn "Cell $cell_index has no baseline row; skipping."
            continue
        end
        baseline_effect = baseline_rows[1, :total_biomass]
        
        # Extract the rows corresponding to species removals.
        removal_df = df[df.sp_removed .!= "none", :]
        if nrow(removal_df) == 0
            @warn "Cell $cell_index has no species removal rows; skipping."
            continue
        end
        
        # Compute effects as the absolute difference from baseline.
        effects = abs.(removal_df.delta_total_biomass .- baseline_effect)
        ARI = mean(effects)
        RV = std(effects)
         
        # Compute Cell Vulnerability Index (CVI)
        # CVI = (ARI / baseline) * (1 + (RV / ARI))
        CVI = (ARI / baseline_effect) * (1 + (RV / ARI))
        
        # Retrieve global network metrics for the cell.
        metrics = compute_food_web_metrics(cell_index; round=false)
        gm = metrics.global_metrics
        density = gm.density
        avg_clust = gm.avg_clustering
        betweenness = gm.global_betweenness
        closeness = gm.global_closeness
        
        # Compute Network Robustness Index (NRI)
        # For simplicity, we use weights w1 = 1, w2 = 1, w3 = 1:
        NRI = density/max_density + avg_clust/max_avg_clust - betweenness/max_betweenness
         
        push!(cell_stability_df, (
           cell_id = cell_index,
           baseline = baseline_effect,
           ARI = ARI,
           RV = RV,
           CVI = CVI,
           density = density,
           avg_clustering = avg_clust,
           global_betweenness = betweenness,
           global_closeness = closeness,
           NRI = NRI
        ))
    end
    println("The 5 highest CVI values are: \n")
    for i in 1:5
        println(sort(cell_stability_df.CVI, rev=true)[i])
    end
    println("The 5 lowest CVI values are: \n")
    for i in 1:5
        println(sort(cell_stability_df.CVI)[i])
    end
    return cell_stability_df
end

# Compute the stability DataFrame (using, for example, all_results_list_even_pi)
cell_stability_df_even_pi = new_measure_cell_stability(all_results_list_even_pi)
cell_stability_df_not_even_pi = new_measure_cell_stability(all_results_list)

# ----------------------------
# 2. Mapping functions to visualize cell-level metrics on the grid
# ----------------------------
function map_cell_metric(
    cell_stability_df::DataFrame, col::Symbol; 
    capped::Bool=false, cap_val::Union{Nothing,Real}=nothing,
    disp::Bool=true, title::Union{Nothing,String}=nothing,
    standardize_by_NPP::Bool=false
)
    # Copy the DA_sum grid as a starting point (convert to float)
    grid = deepcopy(float(DA_sum))
    
    # Replace 0 or 1 entries with NaN or 0 as needed for visualization.
    for i in axes(grid, 1), j in axes(grid, 2)
        if grid[i, j] == 0
            grid[i, j] = NaN
        elseif grid[i, j] == 1.0
            grid[i, j] = NaN
        end
    end
    max_NPP = maximum(npp_DA[.!isnan.(npp_DA)])
    
    # Loop over each cell in the stability DataFrame.
    for i in 1:nrow(cell_stability_df)
         coord = idx[cell_stability_df[i, :cell_id]]
         local_i, local_j = coord[1], coord[2]
         value = cell_stability_df[i, col]
         if capped && cap_val !== nothing && value > cap_val
             grid[local_i, local_j] = cap_val
             if standardize_by_NPP
                npp_value = npp_DA[local_i, local_j]/max_NPP
                grid[local_i, local_j] = value/npp_value
             end
         else
             grid[local_i, local_j] = value
             if standardize_by_NPP
                npp_value = npp_DA[local_i, local_j]/max_NPP
                grid[local_i, local_j] = value/npp_value
             end
         end
    end
    if disp
        display(map_plot(grid; palette = custom_palette, title = isnothing(title) ? identity : title))
    end
    return grid
end

# Map CVI and NRI on the grid.
grid_CVI_even_pi = map_cell_metric(
    cell_stability_df_even_pi, :CVI;
    title = "Cell Vulnerability Index (CVI) with even pi",
    standardize_by_NPP = true
    )#, capped = true, cap_val = 1.18)
grid_CVI_not_even_pi = map_cell_metric(cell_stability_df_not_even_pi, :CVI; title = "Cell Vulnerability Index (CVI) with not even pi")#, capped = true, cap_val = 1.18)
grid_NRI = map_cell_metric(cell_stability_df, :NRI; title = "Network Robustness Index (NRI)")


# ----------------------------
# 3. Scatter plot correlating CVI and NRI across cells
# ----------------------------
begin
    fig3 = Figure(resolution = (600,600))
    ax3 = Axis(fig3[1,1],
        title = "CVI vs. NRI",
        xlabel = "NRI",
        ylabel = "CVI"
    )
    subset = cell_stability_df_even_pi
    scatter!(ax3, subset.NRI, subset.CVI,
             markersize = 10, color = :blue)
    # Optionally, compute and display the Pearson correlation coefficient.
    corr_val = cor(cell_stability_df.CVI, cell_stability_df.NRI)
    text!(ax3, "r = $(round(corr_val, digits=2))", position = (0.05, 0.95),
          align = (:left, :top), color = :black, fontsize = 12)
    display(fig3)
end
