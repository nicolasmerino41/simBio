using Clustering
function merge_community_characteristics(
    all_results_list_even_pi, raw_sensitivity,
    CVI, NRI, CSM, NPP,
    richness, over_under
)
    # Create an empty DataFrame with proper column names
    df = DataFrame(
        i = Int[],
        j = Int[],
        raw_sensitivity = Float64[],
        CVI = Float64[],
        NRI = Float64[],
        CSM = Float64[],
        NPP = Float64[],
        richness = Float64[],
        over_under = Float64[]
    )

    # Iterate through the global list (here we use all_results_list_even_pi)
    for dff in all_results_list_even_pi
        
        cell = dff[1, :cell]
        
        i = idx[cell][1]
        j = idx[cell][2]
        # Collect values for each metric from the corresponding matrices.
        row = [
            i, j,
            raw_sensitivity[i, j], CVI[i, j], NRI[i, j], CSM[i, j],
            NPP[i, j], richness[i, j], over_under[i, j]
        ]

        # (Optional) Convert any NaN to missing if desired.
        row = map(x -> isnan(x) ? NaN : x, row)
        push!(df, row)
    end

    # Remove rows with missing values.
    dropmissing!(df)
    return df
end

# Create the DataFrame of cell features.
df_cells = merge_community_characteristics(
    all_results_list_even_pi,
    Matrix(sensitivity_grid), Matrix(grid_CVI_even_pi), Matrix(grid_NRI), 
    Matrix(csm_grid), Float64.(Matrix(npp_DA)), Matrix(DA_richness_birmmals), Matrix(npp_saturation_grid)
)

# Normalize the features for clustering.
feature_cols = names(df_cells, Not([:i, :j]))
features = Matrix(select(df_cells, feature_cols))
features_norm = (features .- mean(features, dims=1)) ./ std(features, dims=1)

# Run k-means clustering with transposed data.
k = number_clusters  # Choose the number of clusters. You define it in MakingTheCase.jl
clust_res = kmeans(features_norm', k)  # Transpose so that each cell is a column.

# Remove any preexisting :cluster column if necessary.
if :cluster in names(df_cells)
    select!(df_cells, Not(:cluster))
end

df_cells[!, :cluster] = clust_res.assignments

println("features_norm size: ", size(features_norm))  # Should be (3701, 7)
println("df_cells rows: ", nrow(df_cells))

# === Build a grid (matrix) of cluster labels ===
dims = size(sensitivity_grid)  # same as any metric grid
cluster_grid = fill(NaN, dims...)
for row in eachrow(df_cells)
    # Here we assume that idx or the coordinate system is simply (i,j)
    # If your coordinates are stored in :i and :j, then:
    cluster_grid[row.i, row.j] = row.cluster
end

"""
    plot_cluster_map(cluster_grid; resolution=(800,600), title="Community Clusters", xlabel="Column (j)", ylabel="Row (i)", yreversed=true, colormap=:Set1)

Plots the cluster grid as a heatmap. 
- `cluster_grid`: a matrix with cluster labels (one per cell).
- `resolution`: the figure resolution.
- `title`, `xlabel`, `ylabel`: plot labels.
- `yreversed`: if true, reverses the y-axis (common in spatial maps).
- `colormap`: the discrete colormap to distinguish clusters.
Returns the created figure.
"""
function plot_cluster_map(
    cluster_grid; 
    resolution=(800,600), 
    title="Community Clusters",
    colormap=:Set1
)
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1,1],
        title = title,
        yreversed = true
    )
    heatmap!(ax, cluster_grid; interpolate = false, colormap = colormap)
    display(fig)
    return cluster_grid
end

"""
    plot_cluster_scatter(df::DataFrame; xvar::Symbol, yvar::Symbol, cluster_col::Symbol = :cluster, resolution=(600,600), title="", xlabel="", ylabel="", colormap=:Set1, markersize=8)

Creates a scatter plot of two variables (specified by `xvar` and `yvar`) from `df`, 
with points colored according to the values in `cluster_col`. If `title`, `xlabel`, or `ylabel` are not provided, defaults are created.
Returns the created figure.
"""
function plot_cluster_scatter(
    df::DataFrame;
    xvar::Symbol = :NPP, yvar::Symbol = :CSM,
    cluster_col::Symbol = :cluster,
    resolution=(600,600),
    colormap=:custom_palette, markersize=8
)
    # Provide default labels if not given.
    title = "$(xvar) vs. $(yvar) by Cluster"
    xlabel = string(xvar)
    ylabel = string(yvar)
        
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1,1],
        title = title,
        xlabel = xlabel,
        ylabel = ylabel
    )
    scatter!(ax, df[!, xvar], df[!, yvar], 
        color = df[!, cluster_col], 
        colormap = colormap, 
        markersize = markersize
    )
    display(fig)
    return fig
end

using DataFrames, MultivariateStats, Statistics, CairoMakie

"""
    pca_explain_over_under(df::DataFrame; vars = [:raw_sensitivity, :CVI, :NRI, :CSM, :NPP, :richness, :over_under])

Performs a PCA on the columns given by `vars` of the DataFrame `df`.
Then, it computes the Pearson correlation between each PC’s scores and the
`over_under` variable. The function prints the correlation coefficients for
each PC, identifies the PC with the highest absolute correlation with `over_under`,
and prints its loadings (i.e. the contribution of each variable).
It also produces a biplot of PC1 vs. PC2.

Returns a tuple with the PCA model, the PC scores, and the vector of correlations.
"""

using MultivariateStats

function pca_explain_over_under(df::DataFrame; vars = [:raw_sensitivity, :CVI, :NRI, :CSM, :NPP, :richness, :over_under])
    # Select the data for the given variables.
    X = Matrix(select(df, vars))
    
    # Create a Boolean mask for rows that are complete (no missing values).
    complete_inds = map(r -> all(!ismissing, r), eachrow(X))
    X = X[complete_inds, :]
    
    # Manually center X (subtract the mean of each column).
    X_centered = X .- mean(X, dims=1)
    
    # Fit PCA on the centered data.
    pca_model = fit(PCA, X_centered; maxoutdim = size(X,2))
    
    # Compute the scores (transformed data). Each row is an observation.
    T = MultivariateStats.transform(pca_model, X_centered)
    
    # Identify the column index for "over_under" within our variable list.
    over_under_index = findfirst(x -> x == :over_under, vars)
    over_under_vals = X[:, over_under_index]
    
    # Compute correlation between each PC (each column of T) and over_under.
    numPC = size(T,2)
    pc_corr = [cor(T[:,i], over_under_vals) for i in 1:numPC]
    
    println("Correlation of each PC with over_under:")
    for i in 1:numPC
        println("PC$(i): $(round(pc_corr[i], digits=3))")
    end
    
    # Identify the PC with the highest absolute correlation with over_under.
    idx_max = argmax(abs.(pc_corr))
    println("PC$(idx_max) has the highest absolute correlation with over_under: $(round(pc_corr[idx_max], digits=3))")
    
    # Print loadings for that PC.
    loadings = pca_model.projection[:, idx_max]
    println("Loadings for PC$(idx_max):")
    for (var, val) in zip(vars, loadings)
        println("  $(var): $(round(val, digits=3))")
    end

    # === Biplot: Plot PC1 vs. PC2 ===
    fig = Figure(resolution = (800,600))
    ax = Axis(fig[1,1],
        title = "PCA Biplot (PC1 vs. PC2)",
        xlabel = "PC1",
        ylabel = "PC2"
    )
    scatter!(ax, T[:,1], T[:,2], markersize = 8, color = :blue)
    # Scale loadings for display: scale by the maximum absolute value of scores.
    scale1 = maximum(abs.(T[:,1]))
    scale2 = maximum(abs.(T[:,2]))
    for (j, var) in enumerate(vars)
        arrow!(ax, 0, 0, pca_model.projection[j,1]*scale1, pca_model.projection[j,2]*scale2,
            linewidth = 2, color = :red)
        text!(ax, string(var),
            position = (pca_model.projection[j,1]*scale1*1.1, pca_model.projection[j,2]*scale2*1.1),
            align = (:left, :center), color = :red, textsize = 10)
    end
    display(fig)
    
    return (pca_model, T, pc_corr)
end
