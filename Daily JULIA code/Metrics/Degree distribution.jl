using Graphs
using Makie
using Random

# Sample data (replace this with your actual data)
idx_sampled = sample(idx, 16, replace = false)
# Create a grid layout for the plots
fig = Figure(resolution=(800, 800))
axess = [Axis(fig[i, j]; title="Sample $(i + (j - 1) * 4)") for i in 1:4, j in 1:4]

# Loop through sampled indices and plot
for (ax, idx_val) in zip(axess, idx_sampled)
    # Extract abundances for the current sampled index
    abundances = DA_with_abundances[idx_val].a

    # Find non-zero indices based on abundances
    non_zero_indices = findall(x -> x > 0, abundances)

    if isempty(non_zero_indices)
        println("No non-zero indices for $(idx_val), skipping...")
        continue  # Skip if no non-zero indices
    end

    # Create a local adjacency matrix based on non-zero indices
    submatrix = Matrix(iberian_interact_df[non_zero_indices, non_zero_indices])

    # Ensure the submatrix is symmetric
    submatrix = (submatrix .+ submatrix') ./ 2

    # Calculate degrees and prepare for plotting
    g = SimpleGraph(submatrix)
    degrees = degree(g)
    
    degree_counts = Dict{Int, Int}()
    for d in degrees
        degree_counts[d] = get(degree_counts, d, 0) + 1
    end

    # Prepare data for bins
    x = collect(keys(degree_counts))  # Degrees
    y = collect(values(degree_counts)) # Frequencies

    # Plot bar graph
    MK.barplot!(ax, x, y, width=0.8, color=:blue)
    ax.xlabel = "Degree"
    ax.ylabel = "Frequency"
    MK.xlims!(ax, 0, 80)
end

############### MAPPING MEAN DEGREE ################
# Assuming lambda_DA.multiplicative is a grid and idx is defined
# Deep copy of the grid
mean_degree_grid = deepcopy(lambda_DA.multiplicative)

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

    # Ensure the submatrix is symmetric
    submatrix = (submatrix .+ submatrix') ./ 2

    # Calculate degrees for the local graph
    g = SimpleGraph(submatrix)
    degrees = degree(g)

    # Calculate mean degree
    mean_degree = mean(degrees)

    # Assign the mean degree to the corresponding cell in the grid
    mean_degree_grid[idx_val[1], idx_val[2]] = mean_degree
end

begin
    
f = Figure(resolution = (800,500))
ax1 = Axis(f[1,1])
ax2 = Axis(f[1,2])
ax3 = Axis(f[2,1])
ax4 = Axis(f[2,2])
image!(ax1, connectance_map; interpolate=false, colormap=custom_palette, colorrange = (0, 0.15))
image!(ax2, DA_richness; interpolate=false, colormap=custom_palette)
image!(ax3, mean_degree_grid; interpolate=false, colormap=custom_palette)
image!(ax4, clustering_grid; interpolate=false, colormap=custom_palette, colorrange = (0, 0.2))
ax1.yreversed = true
ax2.yreversed = true
ax3.yreversed = true
ax4.yreversed = true
hidexdecorations!(ax1; grid=false)
hideydecorations!(ax1; grid=false)
hidexdecorations!(ax2; grid=false)
hideydecorations!(ax2; grid=false)
hidexdecorations!(ax3; grid=false)
hideydecorations!(ax3; grid=false)
hidexdecorations!(ax4; grid=false)
hideydecorations!(ax4; grid=false)
ax1.title = "Connectance"
ax2.title = "Richness"
ax3.title = "Mean Degree"
ax4.title = "Clustering"

display(f)
end
