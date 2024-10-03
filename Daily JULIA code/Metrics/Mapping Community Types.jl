# Convert species presence dictionary to a matrix
species_presence_list = collect(values(species_presence))
species_presence_matrix = hcat([BitArray(presence) for presence in species_presence_list]...)'

using Clustering, Distances

# Compute Jaccard distance matrix
jaccard_distance_matrix = pairwise(Jaccard(), species_presence_matrix, dims=1)

# Perform hierarchical clustering
hc = hclust(jaccard_distance_matrix, linkage=:average)

# Decide on the number of clusters (k)
k = 2  # Adjust based on your data and analyses
cluster_assignments = cutree(hc; k = k)

# Create a grid to store community types
community_types_grid = fill(NaN, 125, 76)  # Adjust dimensions based on your grid size

# Map cluster assignments to cells
for (i, cell) in enumerate(idx)
    community_types_grid[cell[1], cell[2]] = cluster_assignments[i]
end

using DimensionalData

community_types_grid = DimArray(reshape([NaN for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

for (i, cell) in enumerate(idx)
    community_types_grid[cell] = cluster_assignments[i]
end

# Define a categorical color palette
using ColorSchemes

palette = ColorSchemes.Set1_5.colors  # For 5 clusters; adjust if k changes

# Modify map_plot to handle categorical data
function map_plot_categorical(plot_data; palette = palette, legend = false, flip = true, title = nothing)
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Plot the data using categorical mapping
    plt_obj = image!(ax, plot_data; colormap = palette, colorrange=(1, length(palette)))

    # Reverse y-axis if needed
    ax.yreversed = flip

    # Add legend if requested
    if legend
        # Create a custom colorbar for categorical data
        labels = ["Community $i" for i in 1:length(palette)]
        cb = Colorbar(fig[1, 2], plt_obj, label = "Community Type", ticks = 1:length(palette), ticklabels = labels)
    end

    # Add title if provided
    if !isnothing(title)
        ax.title = title
    end

    hidexdecorations!(ax; grid = false)
    hideydecorations!(ax; grid = false)

    fig
end

# Plot the community types
fig = map_plot_categorical(community_types_grid; palette = palette, title = "Community Types")
display(fig)
