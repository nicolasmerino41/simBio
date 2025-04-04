using Graphs
"""
compute_food_web_metrics(cell_index::Int) -> (global_metrics, species_metrics)

Given the index of a cell, this function:
    1. Extracts the species present using extract_species_names_from_a_cell.
    2. Constructs a sub–network (by sub–setting iberian_interact_NA) corresponding to these species.
    3. Builds a directed graph from the sub–matrix.
    4. Computes global network metrics:
        • density,
        • average degree,
        • average clustering coefficient,
        • (modularity and eigenvector centrality are not computed here, and are set to NaN),
        • average betweenness centrality,
        • average closeness centrality.
    5. Computes species–specific metrics for each node (species) in the sub–network:
        • in–degree, out–degree, total degree,
        • betweenness centrality,
        • closeness centrality,
        • local clustering coefficient.
    
For detailed definitions and the rationale for each metric, see the Supplementary Information.
"""

function compute_food_web_metrics(cell_index::Int; round = false)
    # --- Extract species names from the chosen cell ---
    local_i, local_j = idx[cell_index][1], idx[cell_index][2]
    cell = DA_birmmals_with_pi[local_i, local_j]
    species_names = extract_species_names_from_a_cell(cell)
    if isempty(species_names)
        error("No species found in cell $cell_index.")
    end

    # --- Build the sub-network adjacency matrix ---
    indices = [species_dict[sp] for sp in species_names if haskey(species_dict, sp)]
    submatrix = iberian_interact_NA[indices, indices]
    # Convert to a standard Array and then to a Boolean matrix:
    A_bool = map(x -> x > 0.0, Array(submatrix))
    n = size(A_bool, 1)

    # --- Construct a directed graph from the Boolean matrix ---
    g = SimpleDiGraph(n)
    for i in 1:n, j in 1:n
        if A_bool[i, j]
            add_edge!(g, i, j)
        end
    end
    # Create an undirected version for clustering computation.
    ug = Graph(g)

    # --- Compute Global Metrics ---
    density = Graphs.density(g)
    avg_degree = sum(degree(g, v) for v in vertices(g)) / n
    avg_clustering = mean(local_clustering_coefficient(ug, v) for v in vertices(ug))
    modularity = NaN         # Not computed in this implementation.
    # Compute global betweenness and closeness using available functions.
    global_betweenness = try
        mean(betweenness_centrality(g))
    catch e
        @warn("Betweenness centrality failed; setting to NaN. $e")
        NaN
    end
    global_closeness = try
        mean(closeness_centrality(g))
    catch e
        @warn("Closeness centrality failed; setting to NaN. $e")
        NaN
    end
    global_eigenvector = NaN # Eigenvector centrality not computed.
    global_metrics = (
        density = density,
        avg_degree = avg_degree,
        avg_clustering = avg_clustering,
        modularity = modularity,
        global_betweenness = global_betweenness,
        global_closeness = global_closeness,
        global_eigenvector = global_eigenvector
    )

    # --- Compute Species-Specific Metrics ---
    # We'll collect per-species metrics in a vector of named tuples.
    species_metrics_list = Vector{NamedTuple{(:species_name, :indegree, :outdegree, :total_degree, :betweenness, :closeness, :clustering),NTuple{7, Any}}}(undef, n)
    # Attempt to compute centrality measures; if they fail, assign NaN.
    local_betweenness = try
        betweenness_centrality(g)
    catch
        fill(NaN, n)
    end
    local_closeness = try
        closeness_centrality(g)
    catch
        fill(NaN, n)
    end

    # --- Compute Species-Specific Metrics ---
    for v in vertices(g)
        indeg = indegree(g, v)
        outdeg = outdegree(g, v)
        tot_deg = degree(g, v)
        local_clust = local_clustering_coefficient(ug, v)  # Fixed here
        species_metrics_list[v] = (
            species_name = species_names[v],
            indegree = indeg,
            outdegree = outdeg,
            total_degree = tot_deg,
            betweenness = local_betweenness[v],
            closeness = local_closeness[v],
            clustering = local_clust
        )
    end
    species_metrics = DataFrame(species_metrics_list)
    
    return (global_metrics = global_metrics, species_metrics = species_metrics)
end

    # # Example usage:
    # cell = 459
    # metrics_results = compute_food_web_metrics(4094)  # Replace 1 with the desired cell index.
    # @info "Global metrics for cell $cell:" metrics_results.global_metrics
    # @info "Species-specific metrics:" metrics_results.species_metrics
    # global_metrics = metrics_results.global_metrics
    # species_metrics = metrics_results.species_metrics

function extract_metrics_map(metric)
    DA_metric = deepcopy(float.(DA_sum))
    val = 0
    Threads.@threads for cell in 1:length(idx)
        metrics_results = compute_food_web_metrics(cell)
        value = metrics_results.global_metrics[metric]
        DA_metric[idx[cell][1], idx[cell][2]] = value
    end

    return DA_metric
end

# This function aggregates species-specific metrics across the given cell indices.
function compute_average_species_metrics(cell_indices::AbstractVector{Int})
    all_species_metrics = DataFrame()
    
    # Iterate over the cell indices.
    for cell in cell_indices
        try
            # Call your function to compute the metrics for the cell.
            metrics = compute_food_web_metrics(cell; round=false)
            species_metrics = metrics.species_metrics
            
            # (Optional) Add a column indicating the originating cell.
            species_metrics[!, :cell_id] .= cell
            
            # Vertically concatenate the species metrics from this cell.
            all_species_metrics = vcat(all_species_metrics, species_metrics)
        catch e
            @warn "Error computing metrics for cell $cell: $e"
            continue
        end
    end

    # Group by species (the column "species" holds the species name).
    grouped = groupby(all_species_metrics, :species_name)

    # For each species, compute the mean of each metric across cells.
    species_avg_df = DataFrames.combine(grouped,
        :indegree    => mean  => :mean_indegree,
        :outdegree   => mean  => :mean_outdegree,
        :total_degree=> mean  => :mean_total_degree,
        :betweenness => mean  => :mean_betweenness,
        :closeness   => mean  => :mean_closeness,
        :clustering  => mean  => :mean_clustering,
        nrow         => :cell_count  # How many cells the species appears in.
    )

    return species_avg_df
end

# Example usage:
# Assume that idx is defined and contains the cell indices (or use 1:number_of_cells)
cell_range = 1:5950
@time avg_species_metrics = compute_average_species_metrics(cell_range)

# Display the aggregated DataFrame.
display(avg_species_metrics)

##### COMPUTE THE NRI OF EACH CELL AND PLOT IT #####
function compute_and_map_NRI(; plot = true, title = "NRI", standardise_by_NPP = false, resolution = (600, 600))
    
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

    max_betweenness = 0.0
    max_density = 0.0
    max_avg_clust = 0.0
    # Loop over each cell in the stability DataFrame.
    for i in 1:length(idx)
        local_i, local_j = idx[i][1], idx[i][2]
        metrics = compute_food_web_metrics(i; round=false)
        max_betweenness = max(max_betweenness, metrics.global_metrics.global_betweenness)
        max_density = max(max_density, metrics.global_metrics.density)
        max_avg_clust = max(max_avg_clust, metrics.global_metrics.avg_clustering)
    end
    max_NPP = maximum(npp_DA[.!isnan.(npp_DA)])
    new_npp_da = npp_DA ./ max_NPP
    # Loop over each cell in the stability DataFrame.
    Threads.@threads for i in 1:length(idx)
        local_i, local_j = idx[i][1], idx[i][2]
        metrics = compute_food_web_metrics(i; round=false)
        betweenness = metrics.global_metrics.global_betweenness
        density = metrics.global_metrics.density
        avg_clust = metrics.global_metrics.avg_clustering
        NRI = density/max_density + avg_clust/max_avg_clust - betweenness/max_betweenness
        grid[local_i, local_j] = NRI
        if standardise_by_NPP
            grid[local_i, local_j] = NRI * new_npp_da[local_i, local_j]
        end
    end

    if plot
        fig = Figure(resolution = resolution)
        ax = Axis(fig[1, 1], title = title)
        Makie.heatmap!(ax, grid; interpolate=false, colormap=custom_palette)
        Colorbar(fig[1, 2], colormap=custom_palette)
        ax.yreversed[] = true
        display(fig)
    end
    return grid
end
# compute_and_map_NRI(; plot = true, title = "NRI", standardise_by_NPP = false)
##### TRYING THE FUNCTION #####
if true
    
    if isempty(va)
        DA_density, DA_avg_degree, DA_avg_clustering, DA_global_betweenness, DA_global_closeness =
            extract_metrics_map(1), extract_metrics_map(2), extract_metrics_map(3), extract_metrics_map(5), extract_metrics_map(6)
    end
    va = [
        DA_density,
        DA_avg_degree,
        DA_avg_clustering,
        DA_global_betweenness,
        DA_global_closeness
    ]
    va_names = ["Connectance", "Average Degree", "Average Clustering", "Global Betweenness", "Global Closeness"]
    together = true
    number_of_the_metric = 3
    # max_value = 1
    begin
        if together
            fig = Figure(resolution = (2300, 900))
            ax1 = Axis(fig[1, 1], title = "Connectance")
            ax2 = Axis(fig[1, 2], title = "Average Degree")
            ax3 = Axis(fig[1, 3], title = "Average Clustering")
            ax4 = Axis(fig[1, 4], title = "Global Betweenness")
            ax5 = Axis(fig[1, 5], title = "Global Closeness")

            Makie.heatmap!(ax1, DA_density; interpolate=false, colormap=custom_palette)
            Makie.heatmap!(ax2, DA_avg_degree; interpolate=false, colormap=custom_palette)
            Makie.heatmap!(ax3, DA_avg_clustering; interpolate=false, colormap=custom_palette)
            Makie.heatmap!(ax4, DA_global_betweenness; interpolate=false, colormap=custom_palette)
            Makie.heatmap!(ax5, DA_global_closeness; interpolate=false, colormap=custom_palette)

            ax1.yreversed[] = true
            ax2.yreversed[] = true
            ax3.yreversed[] = true
            ax4.yreversed[] = true
            ax5.yreversed[] = true

            display(fig)
        else
            name = va_names[number_of_the_metric]
            fig = Figure(resolution = (600, 600))
            ax = Axis(fig[1, 1], title = "$name")
            Makie.heatmap!(
                ax, va[number_of_the_metric];
                interpolate=false, colormap=custom_palette,
                # colorrange = (0, 0.15)
                )
            ax.yreversed[] = true
            display(fig)
        end
    end
end

# COMPARING AVG CLUSTERING WITH NRI
if false
begin
    a = map_cell_metric(cell_stability_df_even_pi, :CVI; disp = false, title = "Cell Vulnerability Index (CVI) with even pi")#, capped = true, cap_val = 1.18) 
    b = compute_and_map_NRI(; plot = false, title = "NRI")
    c = compute_and_map_NRI(; plot = false, title = "NRI", standardise_by_NPP = true)

    fig = Figure(resolution = (1200, 600))
    ax = Axis(fig[1, 1], title = "Cell Vulnerability Index (CVI) with even pi")
    Makie.heatmap!(ax, a; interpolate=false, colormap=custom_palette)
    ax.yreversed[] = true
    
    ax2 = Axis(fig[1, 2], title = "NRI")
    Makie.heatmap!(ax2, b; interpolate=false, colormap=custom_palette)
    ax2.yreversed[] = true

    ax3 = Axis(fig[1, 3], title = "NRI * NPP")
    Makie.heatmap!(ax3, c; interpolate=false, colormap=custom_palette)
    ax3.yreversed[] = true

    display(fig)
end
end