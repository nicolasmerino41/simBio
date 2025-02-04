function compute_CVI(all_results_list::Vector{DataFrame})
    # Create a DataFrame to hold the results
    cell_stability = DataFrame(
         cell_id = Int[],
         baseline = Float64[],
         ARI = Float64[],
         RV = Float64[],
         CVI = Float64[],
         density = Float64[],
         avg_clustering = Float64[],
         global_betweenness = Float64[],
         NRI = Float64[],
         CSM = Float64[]
    )
    
    for i in 1:length(all_results_list)
        cell_index = all_results_list[i][1, :cell]
        df = all_results_list[i]
        # Get baseline performance: use sp_removed=="none" (or false)
        baseline_row = df[df.sp_removed .== "none", :]
        if nrow(baseline_row) == 0
            continue
        end
        baseline_value = baseline_row[1, :total_biomass]  # or another indicator

        # Get removal effects:
        removal_df = df[df.sp_removed .!= "none", :]
        if nrow(removal_df) == 0
            continue
        end
        # Let Δ = change in biomass due to removal. You might want absolute value.
        Δ = abs.(removal_df.delta_total_biomass)
        ARI = mean(Δ)
        RV = std(Δ)
        
        # Compute cell vulnerability index (CVI)
        CVI = (ARI / baseline_value) * (1 + (RV / ARI))
        
        # Get global network metrics for the cell
        metrics = compute_food_web_metrics(cell_index; round=false)
        gm = metrics.global_metrics
        density = gm.density
        avg_clust = gm.avg_clustering
        betweenness = gm.global_betweenness
        
        # Compute a Network Robustness Index (NRI)
        # You can adjust the weights (w1, w2, w3)
        w1, w2, w3 = 1.0, 1.0, 1.0
        NRI = w1*density + w2*avg_clust - w3*betweenness
        # println(
        #     "Density: ", round(density, digits=2), 
        #     " Avg Clustering: ", round(avg_clust, digits=2),
        #     " Betweenness: ", round(betweenness, digits=2),
        # )
        # Composite Stability Metric (CSM)
        α, β = 0.5, 0.5  # For instance, equal weighting
        CSM = α*(1 - CVI) + β*(NRI)
        
        # Append row to results
        push!(cell_stability, (
            cell_id = cell_index,
            baseline = baseline_value,
            ARI = ARI,
            RV = RV,
            CVI = CVI,
            density = density,
            avg_clustering = avg_clust,
            global_betweenness = betweenness,
            NRI = NRI,
            CSM = CSM
        ))
    end
    return cell_stability
end

# Then run and analyze:
stability_df = compute_CVI(all_results_list_even_pi)

###### PLOTTING THE DATA ############
function map_CSM(stability_df; plot = true, palette = custom_palette, resolution = (600, 600), title = "CSM")
    
    grid = deepcopy(float(DA_sum))
    for i in axes(grid, 1), j in axes(grid, 2)
        if grid[i, j] == 0
            grid[i, j] = NaN
        elseif grid[i, j] == 1.0
            grid[i, j] = 0.0
        end
    end

    for i in 1:nrow(stability_df)
        coord = idx[stability_df[i, :cell_id]]
        local_i, local_j = coord[1], coord[2]
        grid[local_i, local_j] = stability_df[i, :CSM]
        # println(stability_df[i, :CSM])
    end

    if plot
        fig = Figure(resolution = resolution)
        ax = Axis(fig[1, 1], title = title)
        MK.heatmap!(
            ax, grid, interpolate = false, colormap = palette,
            colorrange = (0.2, 0.8)
        )
        ax.yreversed[] = true
        display(fig)
    end

    return grid
end

map_CSM(stability_df)