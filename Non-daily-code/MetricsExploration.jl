function composition_evaluation(array_output, DA_with_presences, thresh)
    matches = 0
    for i in idx 
        if any(isnan, array_output[i].a)
            throw(ArgumentError("NaN found in array_output"))
        end
        above_threshold = [x > thresh ? 1.0 : 0.0 for x in array_output[i].a]
        matches += sum(above_threshold .== DA_with_presences[i])
    end
    return matches/(length(idx)*256)
end
function richness_evaluation(array_output, DA_with_presences, body_mass_vector)
    matches = 0

    for i in idx
        if any(isnan, array_output[i].a)
            println("NaN found in array_output")
            return 0.0
        end
        # Ensure we are within bounds
        if length(array_output[i].a) != length(body_mass_vector)
            println("Error: Mismatch in lengths at index $i")
            println("Length of array_output[$i].a: ", length(array_output[i].a))
            println("Length of body_mass_vector: ", length(body_mass_vector))
            return 0.0
        end

        # Use the species-specific thresholds
        above_threshold = [array_output[i].a[j] > body_mass_vector[j] ? 1.0 : 0.0 for j in 1:length(array_output[i].a)]
        matches += sum((above_threshold .== 1.0) .& (DA_with_presences[i] .== 1.0)) + sum((above_threshold .== 0.0) .& (DA_with_presences[i] .== 0.0))
    end

    total_presences = sum([sum(DA_with_presences[i] .== 1.0) for i in idx])
    total_absences = sum([sum(DA_with_presences[i] .== 0.0) for i in idx])

    return matches / (total_presences + total_absences)
end

# Example call
richness_evaluation_result = richness_evaluation(
    deepcopy(p[end].birmmals)+deepcopy(p[end].herps), DA_with_presences, body_mass_vector
)
println("Richness Evaluation Result: ", richness_evaluation_result)

function presence_absence_prediction_accuracy(array_output, DA_with_presences, threshold)
    true_positives = 0
    true_negatives = 0
    total_presences = 0
    total_absences = 0
    
    for i in idx
        current_cell = array_output[i].a
        # Debug information
        # println("Current cell: ", current_cell)
        # println("Threshold: ", threshold)
        
        # Predict presence based on threshold
        predicted_presence = current_cell .> threshold
        
        # Actual presence
        actual_presence = DA_with_presences[i]
        
        # Count true positives and true negatives
        true_positives += sum((predicted_presence .== 1) .& (actual_presence .== 1))
        true_negatives += sum((predicted_presence .== 0) .& (actual_presence .== 0))
        
        # Count total presences and absences in the actual data
        total_presences += sum(actual_presence .== 1)
        total_absences += sum(actual_presence .== 0)
    end
    
    # Calculate TPR and TNR
    tpr = total_presences > 0 ? true_positives / total_presences : 0
    tnr = total_absences > 0 ? true_negatives / total_absences : 0
    
    # Return NaN if no presences or absences to avoid division by zero
    if total_presences == 0 || total_absences == 0
        return NaN
    else
        return 2 * (tpr * tnr) / (tpr + tnr)  # Harmonic mean of TPR and TNR
    end
end

# Let's explore how we can start evaluating model outputs
function alpha_richness(model_output::Union{DimArray{MyStructs256{Float64}, 2}, AbstractArray{MyStructs256{Float64}, 2}})
    # Filter and compute richness where b is non-zero
    rich_values = [count(!iszero, cell.a) for cell in model_output if cell.b != 0.0]

    # Print number of non-zero cells
    println("non_zero_cells: ", length(rich_values))

    # Return the mean of the richness values, handle empty array case
    return isempty(rich_values) ? 0.0 : mean(rich_values)
end

# alpha_richness(DA_with_abundances_p)
# fig, ax, pl = myname(Matrix(DA_with_abundances_p), axis = (aspect = DataAspect(),))

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
# gamma_diversity(r)

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
# total_biomass_abundance(r)

function average_connectance(array_output, type::String; time::Int = -1, coords = (0,0))
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

# average_connectance(r, random_iberian_interact, "cell", time = 1, coords = (10,10))

# idx = findall(x -> x == 1.0, DA_sum)
# DA_with_presences = DimArray([fill(0.0, 256) for _ in 1:125, _ in 1:76], (Dim{:a}(1:125), Dim{:b}(1:76)))

# for row in axes(DA_with_abundances, 1), col in axes(DA_with_abundances, 2)
#     if DA_with_abundances[row, col].b != 0.0
#         for i in findall(!iszero, DA_with_abundances[row, col].a)
#             DA_with_presences[row, col][i] = 1.0
#         end
#     end
# end