using Graphs
using CommunityDetection

# Example: Suppose adj_matrix is your adjacency matrix (NxN matrix)
# For an undirected graph:
g = Graph(Matrix(full_symmetric_IM))

# If your graph is directed, use DiGraph:
g = DiGraph(Matrix(iberian_interact_df))

global_clustering_coefficient(g)

############### MAPPING MEAN DEGREE ################
# Assuming lambda_DA.multiplicative is a grid and idx is defined
# Deep copy of the grid
clustering_grid = deepcopy(lambda_DA.multiplicative)
value = 0
# Loop through the sampled indices and calculate mean degrees
for idx_val in idx
    # Extract the corresponding submatrix based on idx
    abundances = DA_with_abundances[idx_val].a

    # Find non-zero indices based on abundances
    non_zero_indices = findall(x -> x > 0, abundances)

    if isempty(non_zero_indices)
        println("No non-zero indices for $(idx_val), skipping...")
        continue  # Skip if no non-zero indices
    end

    # Create a local adjacency matrix based on non-zero indices
    submatrix = Matrix(iberian_interact_df[non_zero_indices, non_zero_indices])

    # Calculate degrees for the local graph
    g = DiGraph(submatrix)
    if global_clustering_coefficient(g) > value
        value = global_clustering_coefficient(g)
    end
    # if global_clustering_coefficient(g) == 0
    #     println("Clustering coefficient is 0 for $(idx_val), skipping...")
    #     println(length(non_zero_indices))
    #     error("Clustering coefficient is 0")
    # end
    # Assign the mean degree to the corresponding cell in the grid
    clustering_grid[idx_val[1], idx_val[2]] = global_clustering_coefficient(g)
end
value