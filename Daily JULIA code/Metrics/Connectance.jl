function count_non_zero(matrix::AbstractMatrix)
    # Count the number of non-zero elements in the matrix
    return sum(x != 0 for x in matrix)
end

function calculate_connectance(matrix::AbstractMatrix)
    # Get the number of rows and columns in the matrix
    num_rows, num_cols = size(matrix)
    
    # Calculate the total number of possible links
    total_links = num_rows * num_cols
    
    # Count the number of actual links (non-zero entries)
    actual_links = count_non_zero(matrix)  # Use the custom count function
    
    # Calculate connectance
    connectance = actual_links / total_links
    
    return connectance
end

# Example usage
example_matrix = [0 1 0; 1 0 1; 0 0 0]
connectance_value = calculate_connectance(full_IM)
println("Connectance: ", connectance_value)
value = 0
spe = Float64
connectance_map = deepcopy(lambda_DA.multiplicative)
above = 0
below_20 = 0
for i in idx
    vector = DA_with_abundances[i[1], i[2]].a
    # Extract indices of non-zero abundances
    non_zero_indices = findall(x -> x > 0, vector)

    # Create a submatrix of iberian_interact_df
    # Using the indices from non_zero_indices to filter the rows and columns
    submatrix = Matrix(iberian_interact_df[non_zero_indices, non_zero_indices])

    # Now you can calculate connectance for the submatrix
    L = count_non_zero(submatrix) # Count non-zero entries in the submatrix
    S = length(non_zero_indices)        # Number of species (rows/columns in submatrix)
    
    if S > 0
        connectance = L / (S^2)  # Calculate connectance
    else
        connectance = 0  # Handle the case where there are no species
    end
    d = value
    if S > 30
        value = max(value, connectance)
    end
    if connectance > 0.15
        above += 1
    end
    if S < 20
        below_20 += 1
    end
    if d != value
        spe = S
    end
    connectance_map[i[1], i[2]] = connectance
    # Do something with the connectance value (e.g., store it, print it, etc.)
    # println("Connectance at $(i[1]), $(i[2]): ", round(connectance, digits=2), " and S: ", S)
end
spe
value
above
below_20

###################################################
###################################################
f = Figure(resolution = (800,500))
ax1 = Axis(f[1,1])
ax2 = Axis(f[1,2])
ax3 = Axis(f[2,1])
image!(ax1, connectance_map; interpolate=false, colormap=custom_palette, colorrange = (0, 0.15))
image!(ax2, DA_richness; interpolate=false, colormap=custom_palette)
image!(ax3, modularity_map; interpolate=false, colormap=custom_palette)
ax1.yreversed = true
ax2.yreversed = true
ax3.yreversed = true
hidexdecorations!(ax1; grid=false)
hideydecorations!(ax1; grid=false)
hidexdecorations!(ax2; grid=false)
hideydecorations!(ax2; grid=false)
hidexdecorations!(ax3; grid=false)
hideydecorations!(ax3; grid=false)
ax1.title = "Connectance"
ax2.title = "Richness"
ax3.title = "Modularity"
display(f)