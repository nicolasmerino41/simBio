begin
    using Graphs
    using DataFrames
    using Statistics
    using Logging

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

    function compute_food_web_metrics(cell_index::Int)
        # --- Extract species names from the chosen cell ---
        local_i, local_j = idx[cell_index][1], idx[cell_index][2]
        cell = DA_birmmals_with_pi[local_i, local_j]
        species_names = extract_species_names_from_a_cell(cell)
        if isempty(species_names)
            error("No species found in cell $cell_index.")
        end

        # --- Build the sub-network adjacency matrix ---
        indices = [species_dict[sp] for sp in species_names if haskey(species_dict, sp)]
        global submatrix = iberian_interact_NA[indices, indices]
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
        global_metrics = (density = density,
                          avg_degree = avg_degree,
                          avg_clustering = avg_clustering,
                          modularity = modularity,
                          global_betweenness = global_betweenness,
                          global_closeness = global_closeness,
                          global_eigenvector = global_eigenvector)

        # --- Compute Species-Specific Metrics ---
        # We'll collect per-species metrics in a vector of named tuples.
        species_metrics_list = Vector{NamedTuple{(:species, :indegree, :outdegree, :total_degree, :betweenness, :closeness, :clustering),NTuple{7, Any}}}(undef, n)
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
            species_metrics_list[v] = (species = species_names[v],
                                    indegree = indeg,
                                    outdegree = outdeg,
                                    total_degree = tot_deg,
                                    betweenness = local_betweenness[v],
                                    closeness = local_closeness[v],
                                    clustering = local_clust)
        end
        species_metrics = DataFrame(species_metrics_list)

        return (global_metrics = global_metrics, species_metrics = species_metrics)
    end

    # Example usage:
    cell = 459
    metrics_results = compute_food_web_metrics(cell)  # Replace 1 with the desired cell index.
    @info "Global metrics for cell $cell:" metrics_results.global_metrics
    @info "Species-specific metrics:" metrics_results.species_metrics
    global_metrics = metrics_results.global_metrics
    species_metrics = metrics_results.species_metrics
end

function extract_metrics_map(metric = 1)
    DA_metric = deepcopy(DA_richness)
    val = 0
    Threads.@threads for cell in idx
        val += 1
        metrics_results  = compute_food_web_metrics(val)
        value = metrics_results.global_metrics[metric]
        DA_metric[cell] = value
    end

    return DA_metric
end
DA_density, DA_avg_degree, DA_avg_clustering, DA_global_betweenness, DA_global_closeness =
 extract_metrics_map(1), extract_metrics_map(2), extract_metrics_map(3), extract_metrics_map(5), extract_metrics_map(6)
begin
fig = Figure(resolution = (1900, 600))
ax1 = Axis(fig[1, 1], title = "Density")
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
end