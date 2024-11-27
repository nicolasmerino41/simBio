using Statistics

function compute_local_morans_i_total_biomass(array_output, position; modified = false, caca = false)
    # Compute combined abundances based on your existing structure
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    else
        error("Invalid combination of parameters")
    end

    # Initialize a grid to store total biomass per cell
    total_biomass_grid = DimArray(reshape([NaN for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

    # Populate the total biomass grid
    for cell in idx
        if !any(isnan, combined_abundances[cell].a)
            total_biomass_grid[cell] = combined_abundances[cell].b
        else
            total_biomass_grid[cell] = NaN
        end
    end

    # Compute the global mean and variance of total biomass
    valid_values = [total_biomass_grid[cell] for cell in idx if !isnan(total_biomass_grid[cell])]
    mean_biomass = mean(valid_values)
    var_biomass = var(valid_values)

    # Initialize a grid to store local Moran's I values
    local_morans_i_grid = DimArray(reshape([NaN for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

    # For each cell, compute local Moran's I
    for cell in idx
        if isnan(total_biomass_grid[cell])
            continue
        end

        # Get the neighbors of the cell
        neighbors = [
            CartesianIndex(cell[1]-1, cell[2]),  # Up
            CartesianIndex(cell[1]+1, cell[2]),  # Down
            CartesianIndex(cell[1], cell[2]-1),  # Left
            CartesianIndex(cell[1], cell[2]+1)   # Right
        ]

        # Filter valid neighbors
        valid_neighbors = [n for n in neighbors if n in idx && !isnan(total_biomass_grid[n])]

        # If there are no valid neighbors, skip the cell
        if isempty(valid_neighbors)
            continue
        end

        # Compute the local Moran's I numerator
        xi = total_biomass_grid[cell] - mean_biomass
        sum_wij_xj = 0.0
        for neighbor in valid_neighbors
            xj = total_biomass_grid[neighbor] - mean_biomass
            sum_wij_xj += xj  # Assuming weights wij = 1 for neighbors
        end

        local_I = (xi / var_biomass) * sum_wij_xj

        # Store the local Moran's I value in the grid
        local_morans_i_grid[cell] = local_I
    end

    return local_morans_i_grid
end

celo = compute_local_morans_i_total_biomass(a, 1; modified = true)
# Your plot data (e.g., local Moran's I grid)
plot_data = celo  # Replace with your actual data

# Call the map_plot function with the desired parameters
map_plot(plot_data; legend = true, title = "Local Moran's I for Total Biomass")

MK.plot(npp_DA)