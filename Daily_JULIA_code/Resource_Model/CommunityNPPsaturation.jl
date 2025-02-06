using GLM

function CommunityNPPsaturation(df;
        scatter::Bool = true,
        map::Bool = true,
        palette = custom_palette,
        resolution_scatter = (600,600),
        resolution_map = (1000,600),
        scatter_title = "NPP vs. Total Biomass",
        map_title = "Residuals (Observed - Predicted Biomass)",
        NPP_aside::Bool = true,
        richness_aside::Bool = false,
        evaluate_richness::Bool = false
)
    
    # --- Fit a linear regression model ---
    lm_model = lm(@formula(biomass_at_the_end ~ NPP), df)
    predicted_biomass = float.(predict(lm_model, df))
    
    # Add predicted biomass and residuals to the DataFrame.
    df[!, :predicted_biomass] = predicted_biomass
    df[!, :residual] = df.biomass_at_the_end .- predicted_biomass

    if evaluate_richness
        # Calculate richness for each cell
        vect = []
        for i in 1:nrow(df)
            cell = df[i, :cell_id]
            local_i, local_j = idx[cell][1], idx[cell][2]
            richness = DA_richness_birmmals[local_i, local_j]
            push!(vect, richness)
        end
        df[!, :richness] = vect
    end 

    # --- Scatter Plot: NPP vs. Total Biomass ---
    if scatter
        fig_scatter = Figure(resolution = resolution_scatter)
        ax_scatter = Axis(fig_scatter[1,1],
            title = scatter_title,
            xlabel = "NPP",
            ylabel = "Total Biomass"
        )
        
        # Normalize richness values for color mapping
        min_richness = minimum(df.richness)
        max_richness = maximum(df.richness)
        norm_richness = (df.richness .- min_richness) ./ (max_richness - min_richness)

        # Choose a colormap (e.g., :viridis, :plasma, :inferno)
        colors = cgrad(:viridis)[norm_richness]

        # Scatter plot with color mapped to richness
        scatter!(
            ax_scatter, df.NPP, df.biomass_at_the_end, 
            color = colors, markersize = 10
        )

        # Plot the regression line (sort by NPP for a smooth line).
        sorted_idx = sortperm(df.NPP)
        lines!(ax_scatter, df.NPP[sorted_idx], predicted_biomass[sorted_idx],
               color = :red, linewidth = 2)
        
        # Plot the regression line (sort by NPP for a smooth line).
        sorted_idx = sortperm(df.NPP)
        lines!(ax_scatter, df.NPP[sorted_idx], predicted_biomass[sorted_idx],
               color = :red, linewidth = 2)
        
        # Compute and display the Pearson correlation coefficient.
        r = cor(df.NPP, df.biomass_at_the_end)
        text!(ax_scatter, "r = $(round(r, digits=2))", position = (0.05, 0.95),
              align = (:left, :top), color = :black, fontsize = 12)
        
        display(fig_scatter)
    end

    # --- Map Plot: Residuals Across the Iberian Grid ---
    if map
        # Create a grid based on the DA_sum template.
        grid = deepcopy(float(DA_sum))
        # Initialize grid: set cells that equal 0 or 1.0 to NaN
        for i in axes(grid, 1), j in axes(grid, 2)
            if grid[i, j] == 0 || grid[i, j] == 1.0
                grid[i, j] = NaN
            end
        end

        # For each cell in the DataFrame, place the residual into the grid using idx.
        for i in 1:nrow(df)
            coord = idx[df[i, :cell_id]]  # Assumes idx[cell_id] returns [local_i, local_j]
            local_i, local_j = coord[1], coord[2]
            grid[local_i, local_j] = df[i, :residual]
        end

        # Plot the residual grid as a heatmap.
        fig_map = Figure(resolution = resolution_map)
        ax_map = Axis(fig_map[1,1],
            title = map_title
        )
        heatmap!(ax_map, grid; interpolate = false, colormap = palette)
        ax_map.yreversed[] = true
        if richness_aside && NPP_aside
            error("richness_aside and NPP_aside cannot be true at the same time")
        end
        # Optionally, show the NPP map alongside.
        if NPP_aside
            ax_map2 = Axis(fig_map[1,2], title = "NPP")
            heatmap!(ax_map2, npp_DA; interpolate = false, colormap = palette)
            # You might not want to reverse the y-axis on the NPP plot.
            ax_map2.yreversed[] = false
        end
        if richness_aside
            ax_map2 = Axis(fig_map[1,2], title = "Richness")
            heatmap!(ax_map2, DA_richness_birmmals; interpolate = false, colormap = palette)
            # You might not want to reverse the y-axis on the NPP plot.
            ax_map2.yreversed[] = true
        end


        display(fig_map)
    end

    return df
end

if false
A = CommunityNPPsaturation(
    Big_P_even_pi_maximised; 
    scatter = true, map = true, NPP_aside = true,
    resolution_scatter = (600,600),
    resolution_map = (1000,400)
    )

#### STUDYING HOW THE SIMULATION CHANGES WITH A CHANGE IN NPP ####
A = single_run(
    1, 0.36666666666664, 0.031007751937984496, 0.5736842105263158, true;
    artificial_pi = true,
    NPP = 840.0718383789062
    )
B = single_run(
    1, 0.36666666666664/10, 0.031007751937984496/10, 0.5736842105263158/10, true;
    artificial_pi = true,
    NPP = 840.0718383789062*10
    )
end