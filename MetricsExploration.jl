# Let's explore how we can start evaluating model outputs
function alpha_richness(model_output::Union{DimArray{MyStructs256{Float64}, 2}, AbstractArray{MyStructs256{Float64}, 2}})
    # Filter and compute richness where b is non-zero
    rich_values = [count(!iszero, cell.a) for cell in model_output if cell.b != 0.0]

    # Print number of non-zero cells
    println("non_zero_cells: ", length(rich_values))

    # Return the mean of the richness values, handle empty array case
    return isempty(rich_values) ? 0.0 : mean(rich_values)
end

alpha_richness(DA_with_abundances_p)
fig, ax, pl = myname(Matrix(DA_with_abundances_p), axis = (aspect = DataAspect(),))

function gamma_diversity(array_output)
    # Assuming each element in array_output is an instance of MyStructs256 and array_output is 3D: [time, x, y]
    num_species = 256
    time_steps = length(array_output)
    total_species_presence = zeros(Int, num_species)

    # Loop over each time step
    for t in 1:time_steps
        
        species_present_this_step = zeros(Bool, num_species)
        # println("time step: ", array_output)

        # Loop over each cell in the 2D matrix for this time step
        for row in axes(array_output[t], 1), col in axes(array_output[t], 2)
                # Safely access the MyStructs256 instance and its 'a' field
                cell = array_output[t][row, col]
                # println("cell type: ", typeof(cell))
                species_present_this_step .|= (cell.a .!= 0)
        end
        
        # Update the total species presence across all time steps
        total_species_presence .+= species_present_this_step
    end

    # Compute the average presence across time steps
    average_species_presence = total_species_presence / time_steps 
    return sum(average_species_presence)
end
gamma_diversity(r)

function total_biomass_abundance(array_output)
    num_species = 256
    time_steps = length(array_output)
    total_species_presence = 0.0

    # Loop over each time step
    for t in 1:time_steps
        
        # Loop over each cell in the 2D matrix for this time step
        for row in axes(array_output[t], 1), col in axes(array_output[t], 2)
                # Safely access the MyStructs256 instance and its 'a' field
                cell = array_output[t][row, col]
                # println("cell type: ", typeof(cell))
                total_species_presence += cell.b
        end
    end
    return total_species_presence/time_steps
end
total_biomass_abundance(r)

function average_connectance(array_output, random_, type::String; time::Int = -1, coords = (0,0))
    time_steps = length(array_output)
    total_connectance = 0.0
    count_non_zero_cells = 0

    if type == "time"
        # Assuming you need `time` parameter here for some computations.
        if time == -1
            throw(ArgumentError("Time required for computations when type is 'time'."))
        end
        for t in time
            for row in axes(array_output[t], 1), col in axes(array_output[t], 2)
                cell = array_output[t][row, col]
                if cell.b != 0.0
                    species_indices = findall(!iszero, cell.a)
                    if !isempty(species_indices)
                        interaction_matrix = iberian_interact_NA[species_indices, species_indices]
                        total_connectance += connectance(interaction_matrix)
                        count_non_zero_cells += 1
                    end
                end
            end
        end
        return count_non_zero_cells > 0 ? total_connectance / count_non_zero_cells : 0.0

    elseif type == "cell"
        row, col = coords
        valid_time_steps = 0
        time_steps = length(array_output)
        for t in 1:time_steps
            if array_output[t][row, col].b != 0.0
                species_indices = findall(!iszero, array_output[t][row, col].a)
                if !isempty(species_indices)
                    interaction_matrix = iberian_interact_NA[species_indices, species_indices]
                    total_connectance += connectance(interaction_matrix)
                    valid_time_steps += 1
                end
            end
        end
        return valid_time_steps > 0 ? total_connectance / valid_time_steps : 0.0
    else
        throw(ArgumentError("Invalid type specified. Use 'time' or 'cell'."))
    end
end

average_connectance(r, random_iberian_interact, "cell", time = 1, coords = (10,10))

idx = findall(x -> x == 1.0, DA_sum)
DA_with_presences = DimArray([fill(0.0, 256) for _ in 1:125, _ in 1:76], (Dim{:a}(1:125), Dim{:b}(1:76)))

for row in axes(DA_with_abundances, 1), col in axes(DA_with_abundances, 2)
    if DA_with_abundances[row, col].b != 0.0
        for i in findall(!iszero, DA_with_abundances[row, col].a)
            DA_with_presences[row, col][i] = 1.0
        end
    end
end
function richness_evaluation(array_output, DA_with_presences)
    matches = 0
    for i in idx 
        above_ten = [x > 1 ? 1.0 : 0.0 for x in array_output[i].a]
        matches += sum(above_ten .== DA_with_presences[i])
    end
    return matches/(length(idx)*256)
end

richness_evaluation(r[1000].state, DA_with_presences)

