############ SHANNON INDEX ###############
##########################################
# Function to compute Shannon index for a given abundance vector
function shannon_index(abundance_vector)
    total_abundance = sum(abundance_vector)
    if total_abundance == 0
        return 0.0  # Handle case where total abundance is zero
    end
    proportions = abundance_vector / total_abundance
    return -sum(p * log(p) for p in proportions if p > 0)
end

# Function to compute the average Shannon index over specified indices
function average_shannon_index(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    
    # Initialize a list to store Shannon indices
    shannon_indices = Float64[]

    # Iterate over the specified indices
    for index in idx
        cell_abundance_vector = combined_abundances[index].a
        if !any(isnan, cell_abundance_vector)
            push!(shannon_indices, shannon_index(cell_abundance_vector))
        end
    end

    # Compute and return the average Shannon index
    return mean(shannon_indices)
end

# Example usage
# Assuming array_output and idx are defined
average_shannon = average_shannon_index(p)
println("Average Shannon Index: ", average_shannon)

# Function to compute the Shannon index matrix
function map_shannon_index(array_output)
    combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA
    shannon_matrix = DimArray(reshape([NaN32 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

    for cell in idx
        cell_abundance_vector = combined_abundances[cell].a
        if any(isnan, cell_abundance_vector)
            shannon_matrix[cell] = NaN32  # Skip cells with NaNs
        else
            shannon_matrix[cell] = shannon_index(cell_abundance_vector)
        end
    end
    return shannon_matrix
end

shannon_matrix = map_shannon_index(p)

# Plot the shannon_matrix to visualize the Shannon index distribution
map_plot(shannon_matrix; palette = custom_palette, legend = true);
###################### SIMPSON INDEX ######################
###########################################################
# Function to compute Simpson index for a given abundance vector
function simpson_index(abundance_vector)
    total_abundance = sum(abundance_vector)
    if total_abundance == 0
        return 0.0  # Handle case where total abundance is zero
    end
    proportions = abundance_vector / total_abundance
    return sum(p^2 for p in proportions)
end

# Function to compute the average Simpson index over specified indices
function average_simpson_index(array_output)
    # Merge birmmals and herps
    combined_abundances = deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)
    
    # Initialize a list to store Simpson indices
    simpson_indices = []

    # Iterate over the specified indices
    for index in idx
        cell_abundance_vector = combined_abundances[index].a
        if any(isnan, cell_abundance_vector)
            continue  # Skip cells with NaNs
        end
        push!(simpson_indices, simpson_index(cell_abundance_vector))
    end

    # Compute and return the average Simpson index
    return mean(simpson_indices)
end

# Function to compute the average inverse Simpson index over specified indices
function average_inverse_simpson_index(array_output)
    # Merge birmmals and herps
    combined_abundances = deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)
    
    # Initialize a list to store inverse Simpson indices
    inverse_simpson_indices = []

    # Iterate over the specified indices
    for index in idx
        cell_abundance_vector = combined_abundances[index].a
        if any(isnan, cell_abundance_vector)
            continue  # Skip cells with NaNs
        end
        push!(inverse_simpson_indices, 1 - simpson_index(cell_abundance_vector))
    end

    # Compute and return the average inverse Simpson index
    return mean(inverse_simpson_indices)
end

# Assuming array_output and idx are defined
average_simpson = average_simpson_index(p)
average_inverse_simpson = average_inverse_simpson_index(p)
println("Average Simpson Index: ", average_simpson)
println("Average Inverse Simpson Index: ", average_inverse_simpson)

############### Mean Trophic Level #################
####################################################
TrophInd = CSV.File("DFs/TLs.csv") |> DataFrame
TrophInd = TrophInd[1:256, 1:2]
TrophInd[findall(x -> x < 1.05, TrophInd[:, 2]), 2] .= 1.0
# TrophInd[:, 2] = TrophInd[:, 2].-1
# TrophInd[256, 2] = 1.0 # For some reason last line had floating point error
rename!(TrophInd, Symbol("Column1") => :Species, Symbol("TL") => :TL)
TrophInd[findall(x -> 1.98 < x < 2.05, TrophInd[:, 2]), 2] .= 2.0
order_indices = indexin(spain_names, TrophInd[:, :Species])
TrophInd = TrophInd[order_indices, :]
TrophInd_vector = TrophInd[:, :TL]

# Function to calculate mean trophic level
function calculate_mean_tl(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    
    meanTL_matrix = DimArray(reshape([NaN32 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector
        # Calculate mean trophic level for the present species
        if sum(presence) > 0 && !any(isnan, abundances)
            meanTL_matrix[cell] = mean(TrophInd_vector[presence])
        else
            meanTL_matrix[cell] = NaN
        end
    end
    # Exclude NaN values from the mean
    return mean(filter(!isnan, meanTL_matrix))
end
meanTL_matrix = calculate_mean_tl(p)

# Calling the function with the specified palette and legend
map_plot(meanTL_matrix; type = "heatmap", palette = :inferno, legend = true)

############## BIOMASS DISTRIBUTION ################
####################################################
function average_bbp(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    bbp_vector = Float64[]
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector
        
        if sum(presence) > 0 && !any(isnan, abundances)
            present_abundances = abundances[presence]
            present_trophic_levels = TrophInd_vector[presence]
            total_biomass = sum(present_abundances)
            weighted_trophic_sum = sum(present_abundances .* present_trophic_levels)
            bbp_vector = push!(bbp_vector, weighted_trophic_sum / total_biomass)
        else
            bbp_vector = push!(bbp_vector, NaN)
        end
    end
    return mean(filter(!isnan, bbp_vector))
end
function generate_bbp_matrix(array_output)
    combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA
    bbp_matrix = DimArray(reshape([NaN32 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector
        
        if sum(presence) > 0
            present_abundances = abundances[presence]
            present_trophic_levels = TrophInd_vector[presence]
            total_biomass = sum(present_abundances)
            weighted_trophic_sum = sum(present_abundances .* present_trophic_levels)
            bbp_matrix[cell] = weighted_trophic_sum / total_biomass
        else
            bbp_matrix[cell] = NaN32
        end
    end
    return bbp_matrix
end

bbp_matrix = generate_bbp_matrix(p)
bbp_vector = average_bbp(p)
# Plot the meanTL_matrix to visualize the biomass balance point distribution
map_plot(bbp_matrix; palette = :inferno, legend = true)

######### RICHNESS SIMILARITY ################
function richness_similarity(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end

    # Create a matrix to store simulated species richness
    simulated_richness = DimArray(reshape([0.0 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
    indices_to_remove = [
    CartesianIndex(34, 9), CartesianIndex(34, 10),
    CartesianIndex(34, 11), CartesianIndex(34, 12), 
    CartesianIndex(34, 13)
    ]

    # Use filter! to remove unwanted indices from idx
    idx_removed = filter!(x -> !(x in indices_to_remove), idx)

    # Calculate presence/absence and simulated richness
    for cell in idx_removed
        if !any(isnan, combined_abundances[cell].a)
            abundances = combined_abundances[cell].a
            presence = abundances .> body_mass_vector
            simulated_richness[cell] = sum(presence)
            # if simulated_richness[cell] != DA_richness[cell]
            #     print("cell is: ", cell, "\n")
            # end
            # println(simulated_richness[cell])
            # println(DA_richness[cell])
        elseif any(isnan, combined_abundances[cell].a)
            simulated_richness[cell] = 0.0
        end
    end

    # Calculate Mean Absolute Error (MAE) between real and simulated richness
    differences = [abs(DA_richness[cell] - simulated_richness[cell]) for cell in idx_removed]

    return mean(differences)
end

# Example usage
# Assuming `p`, `DA_sum`, `body_mass_vector`, `lambda_DA`, and `idx` are defined
mae = richness_similarity(p)
println("Mean Absolute Error (MAE) between real and simulated richness: ", mae)

carnivores_vector = deepcopy(herb_carv_vector)
carnivores_vector[carnivores_vector .== 1.0] .= 0.0
carnivores_vector[carnivores_vector .== 0.00000001] .= 1.0

function alive_predators(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end

    carn_alive = Float64[]
    # Calculate presence/absence and simulated richness
    for cell in idx
        if !any(isnan, combined_abundances[cell].a)
            abundances = combined_abundances[cell].a
            presence = abundances .> body_mass_vector
            
            simulated_predator_richness = sum(presence.*carnivores_vector)/106
            carn_alive = push!(carn_alive, simulated_predator_richness)
        end
    end

    return mean(carn_alive)
end
alive_predators(p)

# CONVERGENCE OF POPULATION SIZES
# Function to compute actual abundances
function compute_actual_abundances(array_output; type = "separated")
    if type == "separated"
        corrected_array_output = map(i -> (deepcopy(array_output[i].herps) + deepcopy(array_output[i].birmmals)) .* lambda_DA, eachindex(array_output))
    elseif type == "combined"
        corrected_array_output = map(i -> (deepcopy(array_output[i].state)) .* lambda_DA, eachindex(array_output))
    end
    return corrected_array_output
end
cucu = compute_actual_abundances(s, type = "combined")

# Updated check_convergence function
function check_convergence(array_output,threshold=0.01, window=20)
    num_steps = length(array_output)
    if num_steps < window
        println("Not enough steps to check convergence")
        return false, NaN
    end

    # Compute actual abundances for all steps
    actual_abundances = compute_actual_abundances(array_output; type = "combined")

    # Compute the mean and variance of the last 'window' steps
    recent_steps = actual_abundances[(end-window+1):end]
    combined_abundances = [sum([sum(cell.a) for cell in step[idx]]) for step in recent_steps]

    # Filter out NaN values
    combined_abundances = filter(!isnan, combined_abundances)

    if length(combined_abundances) < window
        return false, NaN
    end

    mean_abundance = mean(combined_abundances)
    var_abundance = var(combined_abundances)

    # Check if the variance is below the threshold
    stable = var_abundance < threshold * mean_abundance
    return stable, var_abundance
end

# Example usage
stable, var_abundance = check_convergence(s)
println("System stabilized: ", stable, " with variance: ", var_abundance)

# STABILISATION OF DIVERSITY INDICES
function check_diversity_stabilization(array_output, threshold=0.01, window=20)
    num_steps = length(array_output)
    if num_steps < window
        println("Not enough steps to check convergence")
        return false, NaN
    end

    # Compute actual abundances for all steps
    actual_abundances = compute_actual_abundances(array_output)

    # Compute the mean and variance of the last 'window' steps
    recent_steps = actual_abundances[(end-window+1):end]
    shannon_indices = [average_shannon_index(step; modified = true) for step in recent_steps]
    # simpson_indices = [average_simpson_index(step; modified = true) for step in recent_steps]
    
    # Compute the mean and variance
    mean_shannon = mean(shannon_indices)
    var_shannon = var(shannon_indices)
    # mean_simpson = mean(simpson_indices)
    # var_simpson = var(simpson_indices)
    
    # Check if the variance is below the threshold
    stable_shannon = var_shannon < threshold * mean_shannon
    # stable_simpson = var_simpson < threshold * mean_simpson
    
    return stable_shannon, var_shannon
end

# Example usage
stable_shannon, var_shannon = check_diversity_stabilization(p)
println("System stabilized: ", stable_shannon, " with variance: ", var_shannon)

# BETA DIVERSITY
# Function to calculate Bray-Curtis dissimilarity between two abundance vectors
using Distances
function bray_curtis(abundance_vector1, abundance_vector2)
    return pairwise(BrayCurtis(), [abundance_vector1, abundance_vector2])[1, 2]
end
# Function to compute the average Bray-Curtis dissimilarity over all pairs of land cells
function average_bray_curtis(array_output; type = "separated")
    if type == "separated"
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals))
    elseif type == "combined"
        combined_abundances = deepcopy(array_output[end].state).*lambda_DA
    end

    dissimilarities = []
    for i in 1:length(idx)-1
        for j in i+1:length(idx)
            cell1 = combined_abundances[idx[i]].a
            cell2 = combined_abundances[idx[j]].a
            push!(dissimilarities, bray_curtis(cell1, cell2))
        end
    end

    return mean(dissimilarities)
end

# Example usage
avg_bray_curtis = average_bray_curtis(s; type = "combined")

# JACCARD INDEX
# function jaccard_index(presence_vector1, presence_vector2)
#     intersection = sum((presence_vector1 .> body_mass_vector) .& (presence_vector2 .> body_mass_vector))
#     union = sum((presence_vector1 .> 0) .| (presence_vector2 .> 0))
#     return union == 0 ? 0.0 : intersection / union
# end
function jaccard_index(presence_vector1, presence_vector2)
    presence1 = presence_vector1 .> body_mass_vector
    presence2 = presence_vector2 .> body_mass_vector
    intersection = sum(presence1 .& presence2)
    union = sum(presence1 .| presence2)
    return union == 0 ? 0.0 : intersection / union
end

# Function to compute the average Jaccard index over all pairs of land cells
function average_jaccard(array_output)
    combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA

    jaccard_indices = []
    for i in 1:length(idx)-1
        for j in i+1:length(idx)
            cell1 = combined_abundances[idx[i]].a
            cell2 = combined_abundances[idx[j]].a
            push!(jaccard_indices, jaccard_index(cell1, cell2))
        end
    end

    return mean(jaccard_indices)
end

# Example usage
avg_jaccard = average_jaccard(p)
println("Average Jaccard Index: ", avg_jaccard)

# Function to compute Jaccard similarity matrix
function jaccard_similarity_matrix(array_output)
    num_cells = length(idx)
    combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA

    similarity_matrix = fill(NaN, num_cells, num_cells)
    for i in 1:num_cells
        for j in i:num_cells
            cell1 = combined_abundances[idx[i]].a
            cell2 = combined_abundances[idx[j]].a
            similarity = jaccard_index(cell1, cell2)
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
        end
    end

    return similarity_matrix
end

# Example usage
similarity_matrix = jaccard_similarity_matrix(p)
MK.plot(similarity_matrix)