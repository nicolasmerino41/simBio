using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio")) 
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
###################################FUNCTIONS###########################
#######################################################################
############### GROWTH ########################
function growth(abundance, self_regulation, K)
    # Assuming growth is a simple function of abundance, self-regulation, and carrying capacity (K)
    return self_regulation * K * abundance * (1 - abundance / K)
end
cucu = growth(vabundances[1][1][1], self_regulation[1], K[1])
########################## SIMBIO ####################################
######################################################################
function simbio(K, self_regulation, abundances, A_matrix, num_steps=nothing)
    num_species = size(abundances, 1)
    num_steps = isnothing(num_steps) ? size(abundances, 2) : num_steps

    if num_species != length(self_regulation)
        throw(ArgumentError("The number of species must be equal to the length of the self_regulation vector"))
    end
    if num_species != length(K)
        throw(ArgumentError("The number of species must be equal to the length of the K vector"))
    end

    for t in 1:num_steps-1
        for species in 1:num_species
            growth_rate = growth(abundances[species, t], self_regulation[species], K[species])
            predation_effect = sum(A_matrix[species, :] .* abundances[:, t] .* abundances[species, t])
            abundances[species, t + 1] = max(0, abundances[species, t] + growth_rate + predation_effect)
            if isnan(abundances[species, t + 1])
                throw(ArgumentError("Abundances at t+1 is NA"))
            end
        end
    end
    
    return abundances
end

########################### SIMBIO EFFICIENT ############################
#########################################################################
function simbio_efficient(K, self_regulation, abundances, A_matrix, num_steps)
    num_species = length(abundances)
    steps = num_steps 
    # if num_species != length(self_regulation)
    #     throw(ArgumentError("The number of species must be equal to the length of the self_regulation vector"))
    # end
    # if num_species != length(K)
    #     throw(ArgumentError("The number of species must be equal to the length of the K vector"))
    # end
    # println(typeof(abund))
    # println(typeof(K))
    # println(typeof(self_regulation))	
    for t in 1:steps
        new_abundances = similar(abundances)
        for species in 1:num_species
            # println("species", species)
            growth_rate = growth(abundances[species], self_regulation[species], K[species])
            if false
                println("gr", growth_rate)
            end
            predation_effect = sum(A_matrix[species, :] .* self_regulation .* abundances[] .* abundances[species])
            if false
                println("pred", predation_effect)
            end
            new_abundances[species] = max(0, abundances[species] + growth_rate + predation_effect)
            if false
                println(species, "new_ab", new_abundances[species])    
            end
            if any(isnan.(new_abundances[species]))
                println("new_abundances t+1 is NA in cell", cell, "and time t", t, "and species", species)
                throw(ArgumentError("new_abundances t+1 is NA in cell"))
            end
        end
        abundances = new_abundances
        # println(t)
    end
    
    return abundances
end

######## SIMBIO_EFFICIENT_MAP #################################
###############################################################
function simbio_efficient_map!(K, self_regulation, abundances, A_matrix, num_steps)
    num_species = length(abundances)
    
    # Check if lengths match
    if num_species != length(self_regulation) || num_species != length(K)
        throw(ArgumentError("The number of species must be equal to the length of the self_regulation and K vectors"))
    end

    for t in 1:num_steps
        # Calculate growth rates
        growth_rates = growth.(abundances, self_regulation, K)
        # println(growth_rates)
        # Calculate total change in abundances
        delta_abundances = growth_rates .+ sum(A_matrix .* self_regulation .* adjoint(abundances) .*abundances, dims=2)[:, 1]
        # println(delta_abundances[])
        # Update abundances in-place
        @. abundances += delta_abundances
        # println(abundances)
    end
    
    return abundances
end

# Define your matrix
matrix = [1 1 1; 1 1 1; 1 1 1]

# Define your vector
vector = [0.01, 0.02, 0.03]
vector_ = reshape(vector, 1, :)
pepe = matrix .* vector
# Sum rows of a matrix and return a vector
row_sums = sum(pepe, dims=2)[:, 1]
######################## Size_selection kernel ########################
#######################################################################
function size_selection_kernel(predator_mass, prey_mass, sd, beta)
    intensity = exp(-((log(float(predator_mass)) / (beta * float(prey_mass)))^2.0) / (2.0 * float(sd)^2.0))
    return float(intensity)
end
beta = float(3)

gbif_sizes = CSV.read(joinpath(meta_path, "Lists\\gbif_sizes.csv"), DataFrame)[:, 2:end]
####################### DISTANCE DECAY ###############################
######################################################################
function distance_decay(distance, ro, abundance, K)
    exp(-distance^2 / (2 * ro^2)) * (abundance / K)
end 
distance_decay(0.2, 0.05, 1000, 1200)
####################### CONNECTANCE ###############################
######################################################################
function connectance(matrix)
    num_species = size(matrix, 1)
    num_links = count(!iszero, matrix)  # Count the number of non-zero elements
    total_possible_links = num_species^2  # Total possible links in a square matrix
    
    return num_links / total_possible_links
end
###################### ZERO DIAGONAL ##################################
######################################################################
function zero_out_diagonal!(matrix)
    n = size(matrix, 1)
    for i in 1:n
        matrix[i, i] = 0
    end
    return matrix
end
#################### FILL DIAGONAL ###################################
######################################################################
# Assuming there exists a function with name `fill_diagonal!` to fill the diagonal of a matrix.
# If not, it needs to be defined as follows:
function fill_diagonal!(mat, val)
    for i in 1:min(size(mat)...)
        mat[i, i] = val
    end
end

web = CSV.read(joinpath(meta_path, "Metaweb_data\\TetraEU_pairwise_interactions.csv"), DataFrame)

web = DataFrame(predator = web.sourceTaxonName, prey = web.targetTaxonName)

web.predator = string.(web.predator)
web.prey = string.(web.prey)
web = web[:, [:predator, :prey]]

unique_predators = unique(web.predator)
unique_preys = unique(web.prey)

x = vcat(unique_predators, unique_preys)
unique_species = unique(x)

# Read the CSV file
diets = CSV.File(joinpath(meta_path, "Metaweb_data\\TetraEU_generic_diet.csv")) |> DataFrame
diets = hcat(diets.sourceTaxonName, diets.targetGenericItemName)

# fn = download("C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling\\Rasters\\iberian_temperature.tif")
temp = Rasters.Raster(joinpath(meta_path, "Rasters\\iberian_temperature.tif"))
Plots.plot(temp);
using Plots

Amph = CSV.read(joinpath(meta_path, "Atlas_data/DB_Amphibians_IP.txt"), delim='\t', DataFrame)
Bird = CSV.read(joinpath(meta_path, "Atlas_data/DB_Birds_IP.txt"), delim='\t', DataFrame)
Mamm = CSV.read(joinpath(meta_path, "Atlas_data/DB_Mammals_IP.txt"), delim='\t', DataFrame)
Rept = CSV.read(joinpath(meta_path, "Atlas_data/DB_Reptiles_IP.txt"), delim='\t', DataFrame)
# data = load("Abundance lists ATLAS\\eq_dens_5928cells.RData")
amphibian_names = names(Amph)
reptile_names = names(Rept)
mammal_names = names(Mamm)
bird_names = names(Bird)
spain_fauna = append!(amphibian_names, reptile_names, mammal_names, bird_names)

# Refactored code to merge columns of DataFrames for Spanish fauna, keeping UTMCODE just once
spanish_fauna = hcat(Amph, Rept[:, Not(:UTMCODE)], Mamm[:, Not(:UTMCODE)], Bird[:, Not(:UTMCODE)])

# Merge the `web` DataFrame with `spain_fauna` using inner join on 'predator' from `web` and 'species' from `spain_fauna`
merged_web = innerjoin(web, DataFrame(species=spain_fauna), on=(:predator => :species))
# Filter the merged DataFrame based on the prey column to include only species found in Spain
merged_web = merged_web[in.(merged_web.prey, Ref(spain_fauna)), :]
# Obtaining the species names that are at least predator/prey
unique_species_in_web = unique(vcat(merged_web.predator, merged_web.prey))
println("There are ", length(unique_species_in_web), " unique species in the food web")

# Initializing an empty matrix with zeros for the Iberian species interaction
n = length(unique_species_in_web)
iberian_interact_matrix = zeros(Int, n, n)
iberian_interact_matrix = NamedArray(iberian_interact_matrix, (unique_species_in_web, unique_species_in_web))
## Creating a mapping from species names to matrix indices
species_to_index = Dict(zip(unique_species_in_web, 1:n))
species_names = collect(keys(species_to_index))
#Filling the matrix with 1s where there are predator-prey interactions
for i in 1:nrow(merged_web)
    iberian_interact_matrix[merged_web.predator[i], merged_web.prey[i]] = 1
end
iberian_interact_matrix = float(iberian_interact_matrix)
# Count the amount of 1s in the iberian_interact_matrix
interaction_count = sum(iberian_interact_matrix .== 1)
println("The number of 1s in the iberian_interact_matrix is $interaction_count")

# Turn iberian_interact_matrix into a DataFrame
# Convert the iberian_interact_matrix into a DataFrame with appropriate column and row names
iberian_interact_df = DataFrame(iberian_interact_matrix, species_names)

abundances = rand(Float64, 250, 5)
abundances[:, 2:end] .= 0
K = rand(50:1000, 1, 250)
ro = 0.05
predation_matrix = rand(250,250)
self_regulation = rand(1, 250)*0.001
diagm(predation_matrix) = self_regulation

lista = CSV.File(joinpath(meta_path, "listamatrices.csv")) |> DataFrame
############################## UNNECESSARY ################################
###########################################################################
###########################################################################


####################### Re-doing the list abundances ####################
#########################################################################
# Remove the first column and store its values in species_names
species_names = abundances_df[!, 1]
abundances_df1 = select(abundances_df, Not(1))

# Define the number of columns per sublist
cols_per_sublist = 120

# Get the column names of abundances_df
col_names = names(abundances_df1)

# Initialize an empty list to store the submatrices
submatrices = []

# Iterate over the column names and group them into submatrices of 120 columns each
for i in 1:cols_per_sublist:length(col_names)
    # Extract a sublist of column names
    sublist_cols = col_names[i:min(i+cols_per_sublist-1, end)]
    
    # Extract the corresponding columns from the DataFrame
    sublist_df = select(abundances_df, sublist_cols)
    
    try
        # # Convert the sublist DataFrame to a matrix
        sublist_matrix = float(Matrix(sublist_df))
        
        # Append the sublist to the list of submatrices
        push!(submatrices, sublist_matrix)
    catch err
        println("Error converting sublist to matrix: $err")
    end
end

# Convert submatrices into a list of vectors containing only the first column of each matrix
vector_list = [submatrix[:, 1] for submatrix in submatrices]

############## Creating random K's per cell ##################
# Define the number of species
num_species = length(unique_species_in_web)

# Define the number of cells
num_cells = length(submatrices)

# Initialize an empty list to store the K vectors for each cell
K_list = []
K_list = Vector{Vector{Float64}}()
# Generate a random K vector for each cell (each element of submatrices)
for _ in 1:num_cells
    # Generate random K values for each species in a cell
    k_cell = rand(10000.0:100000.0, num_species)

    # Append the K vector for the current cell to the K list
    push!(K_list, k_cell)
end
# Define other parameters needed for the simbio function

ro = 0.5  # Example value for ro (you can adjust as needed)
self_regulation = 1  # Example self-regulation values

# Ensure all values in the matrices contained in `submatrices` are floats
for i in 1:length(submatrices)
    submatrices[i] = float(submatrices[i])
end

iberian_interact_matrix = float(iberian_interact_matrix)
function turn_adj_into_inter(adjacencyy, sigma, epsilon)
    adjacency = deepcopy(adjacencyy)
    epsilon = float(epsilon)
    u = adjacency
    for i in names(adjacency, 1)
        for j in names(adjacency, 2)
            if adjacency[i, j] != 0.0 && i != j && adjacency[i, j] > 0.0 && adjacency[j, i] == 0.0
                predator_mass = Float64(gbif_sizes[gbif_sizes.species .== i, :bodyMass][1])
                prey_mass = Float64(gbif_sizes[gbif_sizes.species .== j, :bodyMass][1])
                sd = Float64(gbif_sizes[gbif_sizes.species .== i, :sigma][1])  # Use the sigma value for the predator
                # println(length(predator_mass))
                # Calculate interaction strength based on size-selection kernel
                kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                intensity = max(0.001*sd, kernel)
                # We need to create a Normal distribution first
                normal_dist = Normal(0, sigma*intensity)
                x = round(abs(rand(normal_dist)), digits = 20)
                u[i, j] = x / epsilon
                u[j, i] = -x
            elseif adjacency[i, j] != 0.0 && i != j && adjacency[i, j] > 0.0 && adjacency[j, i] > 0.0
                predator_mass = Float64(gbif_sizes[gbif_sizes.species .== i, :bodyMass][1])
                prey_mass = Float64(gbif_sizes[gbif_sizes.species .== j, :bodyMass][1])
                sd = Float64(gbif_sizes[gbif_sizes.species .== i, :sigma][1])  # Use the sigma value for the predator
                # println(length(predator_mass))
                # Calculate interaction strength based on size-selection kernel
                kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                intensity = max(0.001 * sigma, kernel)
                
                # Draw from a semi-Gaussian distribution
                normal_dist = Normal(0, sigma * intensity)
                x = round(abs(rand(normal_dist)), digits = 20)
                u[i, j] = x / 4.0
            elseif i == j
                u[i, j] = -self_regulation
            end
        end
    end
    adjacency = u
    return adjacency
end
@time begin
    iberian_interact_matrix_finished = turn_adj_into_inter(iberian_interact_matrix, 0.01, beta) 
end

num_cells = 20
num_steps = 120
@time begin
    # Refactored code to use a more concise and Julia-idiomatic approach
A_matrix_list = [turn_adj_into_inter(iberian_interact_matrix, 0.01, beta) for _ in 1:num_cells] 
end

P_matrix_list = [turn_adj_into_inter(iberian_interact_matrix, 0.01, beta) for _ in 1:num_cells]

# Preallocate the simbio_results array with the correct size to improve performance
simbio_results_matr = Vector{Matrix{Float64}}(undef, num_cells)
simbio_results_vect = Vector{Vector{Float64}}(undef, num_cells)
# Running the simulation in a preallocated array
# @time for cell in 1:num_cells
#     simbio_results_matr[cell] = simbio(K_list[cell], self_regulation, submatrices[cell], P_matrix_list[cell])
#     println(cell)
# end

# # Running the simulation in a preallocated array
# @time for cell in 1:num_cells
#     K = K_list[cell]
#     abundances = vector_list[cell]
#     A_matrix = P_matrix_list[cell]
#     simbio_results_vect[cell] = simbio_efficient(K_list[cell], self_regulation, abundances, A_matrix, num_steps)
#     println(cell)
# end
###################### Non equilbrium simulations, building A_matrices ##########################
#################################################################################################
# Creating submatrices per each cell
# Function that extracts the species present in each cell
function species_in_cell(df::DataFrame, cell_index::Int, relevant_species::Set{String})
    filter(col -> col in relevant_species && df[cell_index, col] == 1, names(df))
end

# Filter out the columns from `spanish_fauna` that are not in `unique_species_in_web`
filtered_fauna = spanish_fauna[:, intersect(names(spanish_fauna), unique_species_in_web)]

# Assuming `filtered_fauna` DataFrame's first column is the cell index
num_cells = nrow(filtered_fauna)
relevant_species_set = Set(unique_species_in_web)

# Vector holding species list for each cell, using array comprehension
species_presence_list = [species_in_cell(filtered_fauna, i, relevant_species_set) for i in 1:num_cells]
species_presence_list[1]

################# Creating A_matrix list ############################################
#####################################################################################
# Preallocate the list to store the matrices
A_matrix_list = Vector{Matrix{Float64}}(undef, num_cells)

# Building the list of matrices
if false
for (i, species_list) in enumerate(species_presence_list)
    # Initialize a matrix with all zeros
    species_matrix = zeros(Int, length(species_list), length(species_list))
    
    # Create a mapping from species names to matrix indices
    species_to_index = Dict(zip(species_list, 1:length(species_list)))
    
    # Fill the matrix with 1's where there are predator-prey interactions
    for row in eachrow(merged_web)
        predator_index = get(species_to_index, row.predator, nothing)
        prey_index = get(species_to_index, row.prey, nothing)
        if predator_index !== nothing && prey_index !== nothing
            species_matrix[predator_index, prey_index] = 1
        end
    end
    
    # Assign the matrix to the corresponding cell in the list
    A_matrix_list[i] = species_matrix
end
end # Remove if needed
# # Serialize object and save to file
# serialize("A_matrix_list.jls", A_matrix_list)
# Load the serialized object from the file and cast to the specific type
A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}

################# Creating fullA_matrix list #######################################
####################################################################################
# Preallocate the list to store the matrices
fullA_matrix_list = Vector{Matrix{Float64}}(undef, num_cells)
# Create a mapping from all unique species to global indices
global_species_to_index = Dict(zip(unique_species_in_web, 1:length(unique_species_in_web)))

# Building the list of full matrices
if false # Remove if needed
for i in 1:num_cells
    # Initialize a matrix for all species, with zeros
    full_species_matrix = zeros(Int, length(unique_species_in_web), length(unique_species_in_web))
    
    # Get the species list for the current cell
    species_list = unique_species_in_web
    # Create a local index mapping for species present in the cell
    local_species_to_index = Dict(zip(species_list, 1:length(species_list)))
    
    # Fill the matrix with 1's where there are predator-prey interactions
    for row in eachrow(merged_web)
        global_predator_index = get(global_species_to_index, row.predator, nothing)
        global_prey_index = get(global_species_to_index, row.prey, nothing)
        
        # Check if both predator and prey are present in the current cell's species list
        if (global_predator_index !== nothing && global_prey_index !== nothing) &&
           (row.predator in species_list && row.prey in species_list)
            # Use local indices to set the value in the cell's matrix
            local_predator_index = local_species_to_index[row.predator]
            local_prey_index = local_species_to_index[row.prey]
            full_species_matrix[local_predator_index, local_prey_index] = 1
        end
    end
    
    # Assign the full matrix to the corresponding cell in the list
    fullA_matrix_list[i] = full_species_matrix
end
end # Remove if needed

# # Serialize object and save to file
# serialize("fullA_matrix_list.jls", fullA_matrix_list)
# Load the serialized object from the file and cast to the specific type
fullA_matrix_list = deserialize(joinpath(meta_path, "fullA_matrix_list.jls"))::Vector{Matrix{Float64}}
################ Building IM's #################################################
#################################################################################
sigmas = [1.000, 0.100, 0.010, 0.001]
epsilon = 1
num_cells = 10
# Create a list of lists, each sub-list with num_cells elements, for each sigma
IM_list = [Vector{Matrix{Float64}}(undef, num_cells) for _ in 1:length(sigmas)]
# Julia equivalent for the provided R code snippet
@time for sigma_index in 1:length(sigmas)
    sigma_value = sigmas[sigma_index]
    for cell in 1:num_cells
        short_A_matrix = A_matrix_list[cell]
        u = copy(short_A_matrix)
        for i in 1:size(short_A_matrix, 1)
            for j in 1:size(short_A_matrix, 2)
                # Get predator and prey species names
                predator_species = species_presence_list[cell][i]
                prey_species = species_presence_list[cell][j]
                
                # Extract predator and prey masses and standard deviation from gbif_sizes
                predator_mass = gbif_sizes[gbif_sizes.species .== predator_species, :bodyMass][1]
                prey_mass = gbif_sizes[gbif_sizes.species .== prey_species, :bodyMass][1]
                sd = gbif_sizes[gbif_sizes.species .== predator_species, :sigma][1]

                # Calculate kernel and intensity
                kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                intensity = max(0.001 * sigma_value, kernel)
                    
                # Generate random value from normal distribution
                normal_dist = Normal(0, sigma_value * intensity)
                x = round(abs(rand(normal_dist)), digits = 20)
                    
                # Update interaction matrix based on conditions
                if short_A_matrix[i, j] != 0 && i != j && short_A_matrix[i, j] > 0 && short_A_matrix[j, i] == 0
                    u[i, j] = x / epsilon
                    u[j, i] = -x
                elseif short_A_matrix[i, j] != 0 && i != j && short_A_matrix[i, j] > 0 && short_A_matrix[j, i] > 0
                    u[i, j] = x / 4
                end
            end
        end
        IM_list[sigma_index][cell] = u
        fill_diagonal!(A_matrix_list[cell], -1)
    end
    println("Done with sigma = ", sigma_value)
end

# # Serialize object and save to file
# serialize("IM_list.jls", IM_list)
# Load the serialized object from the file and cast to the specific type
# IM_list = deserialize("IM_list.jls")::Vector{Vector{Matrix{Float64}}}

################ Builiding Threaded IM's #########################################
##################################################################################
sigmas = [1.000, 0.100, 0.010, 0.001]
epsilon = 1
num_cells = 100
if false # Remove if needed
# Create a list of lists, each sub-list with num_cells elements, for each sigma
fullIM_list = [Vector{Matrix{Float64}}(undef, num_cells) for _ in 1:length(sigmas)]

function fullIM_builder(sigmas, A_matrix_list)
    # Corrected parallelized operation with specified sigma_value for each thread
    for sigma_index in 1:length(sigmas)
        # Get the sigma value for the current thread
        local sigma_value = sigmas[sigma_index]
        
        # Iterate over each cell in A_matrix_list
        for cell in 1:num_cells
            # Extract the short A_matrix for the current cell
            short_A_matrix = A_matrix_list[cell]
            
            # Create a copy of short_A_matrix (unnecessary if short_A_matrix is not reused)
            u = copy(short_A_matrix)
            
            # Parallel loop over rows of short_A_matrix
            for i in 1:size(short_A_matrix, 1)
                # Parallel loop over columns of short_A_matrix
               for j in 1:size(short_A_matrix, 2)
                    # Get predator and prey species names
                    predator_species = unique_species_in_web[i]
                    prey_species = unique_species_in_web[j]
                     
                    # Extract predator and prey masses and standard deviation from gbif_sizes
                    predator_mass = gbif_sizes[gbif_sizes.species .== predator_species, :bodyMass][1]
                    prey_mass = gbif_sizes[gbif_sizes.species .== prey_species, :bodyMass][1]
                    sd = gbif_sizes[gbif_sizes.species .== predator_species, :sigma][1]

                    # Calculate kernel and intensity
                    kernel = size_selection_kernel(predator_mass, prey_mass, sd, beta)
                    intensity = max(0.001 * sigma_value, kernel)
                            
                    # Generate random value from normal distribution
                    normal_dist = Normal(0, sigma_value * intensity)
                    x = round(abs(rand(normal_dist)), digits = 20)
                            
                    # Update interaction matrix based on conditions
                    if short_A_matrix[i, j] != 0 && i != j && short_A_matrix[i, j] > 0 && short_A_matrix[j, i] == 0
                        u[i, j] = x / epsilon
                        u[j, i] = -x
                    elseif short_A_matrix[i, j] != 0 && i != j && short_A_matrix[i, j] > 0 && short_A_matrix[j, i] > 0
                        u[i, j] = x / 4
                    end
                        
                end
            end
            
            # Store updated interaction matrix in IM_list
            fullIM_list[sigma_index][cell] = u
            # Fill diagonal of A_matrix_list with -1
            # fill_diagonal!(A_matrix_list[cell], -1)
        end
        println("Done with sigma = ", sigma_value)
    end
end
@time fullIM_builder(sigmas, fullA_matrix_list)
end
# serialize("fullIM_list.jls", fullIM_list)
fullIM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
for sigma in 1:length(fullIM_list) 
    for cell in 1:length(fullIM_list[sigma])
        fullIM_list[sigma][cell] = zero_out_diagonal!(fullIM_list[sigma][cell])
    end
end
aaa = fullIM_list[3][1]
############################# EMPTY ABUNDANCES #############################
# Build a list of matrices for each cell
num_species = 256
num_steps = 120
mabundances = [NamedArray(zeros(num_species, num_steps), (unique_species_in_web, 1:num_steps)) for _ in 1:num_cells] 
vabundances = [NamedArray(zeros(num_species), (unique_species_in_web)) for _ in 1:num_cells]
vupdated_abundances = [zeros(num_species) for _ in 1:num_cells]
vupdated_abundances = [copy(vupdated_abundances) for _ in 1:length(sigmas)]
######################## FOR MABUNDANCES ################################
for (index, species_list) in enumerate(species_presence_list[1:num_cells])
    matrix = mabundances[index]

    for species in species_list
        species_index = findfirst(isequal(species), unique_species_in_web)
        matrix[species_index, 1] = rand(10.0:100.0)
    end
    mabundances[index] = matrix
end
mabundances = [copy(mabundances) for _ in 1:length(sigmas)]
################## FOR VABUNDANCES ####################################
# Refactored code block to assign random values between 0 and 10 for species present in each cell
# and assign species names to the vectors in vabundances.

# Assuming species_names is a vector of species names
for (index, species_list) in enumerate(species_presence_list[1:num_cells])
    vector = vabundances[index]

    for species in species_list
        species_index = findfirst(isequal(species), unique_species_in_web)
        vector[species_index] = rand(10:100)
    end
    vabundances[index] = vector # Convert NamedVector to Vector
end
vabundances = [copy(vabundances) for _ in 1:length(sigmas)]
vabundances_notnamed = [zeros(num_species) for _ in 1:num_cells]
for cell in 1:num_cells
    vabundances_notnamed[cell] = vabundances[1][cell][:] 
end
vabundances_notnamed = [copy(vabundances_notnamed) for _ in 1:length(sigmas)]

K1 = rand(100:1000.0, num_species)
self_regulation = fill(0.0001, num_species)

############ SIMBIO #########################
if false
updated_abundances = deepcopy(mabundances)
@time for sigma in 1:length(sigmas)
    for cell in 1:20 
        updated_abundances[sigma][cell] = simbio(K1, self_regulation, mabundances[sigma][cell], fullIM_list[sigma][cell]) 
    end 
end
end
############# SIMBIO EFFICIENT ################
if false
vupdated_abundances = deepcopy(vabundances)
@time for sigma in 1:length(sigmas)
    for cell in 1:10
        vupdated_abundances[sigma][cell] = simbio_efficient(K1, self_regulation, vabundances[sigma][cell], fullIM_list[sigma][cell], num_steps)
    end
end
end
############ SIMBIO EFFICIENT THREADS ################
if true
vupdated_abundances = deepcopy(vabundances)
@time Threads.@threads for sigma in 1:length(sigmas)
    for cell in 1:10
        # Perform the same operation on each thread
        vupdated_abundances[sigma][cell] = simbio_efficient(K1, self_regulation, vabundances[sigma][cell], fullIM_list[sigma][cell], num_steps)
    end
end
end
############### SIMBIO_EFFICIENT_MAP ################
if true
    vupdated_abundances = deepcopy(vabundances_notnamed)
@time for sigma in 1:length(sigmas)
    for cell in 1:100
    # println("eso ", typeof(vabundances_notnamed[sigma][cell]))
    vupdated_abundances[sigma][cell] = simbio_efficient_map!(K1, self_regulation, vabundances_notnamed[sigma][cell], fullIM_list[sigma][cell], num_steps)
    end
end
end
############### SIMBIO_EFFICIENT_MAP THREADS ################
if false
    vupdated_abundances = deepcopy(vabundances)
@time Threads.@threads for sigma in 1:length(sigmas)
    for cell in 1:100
    # println("eso ", typeof(vabundances_notnamed[sigma][cell]))
    vupdated_abundances[sigma][cell] = simbio_efficient_map!(K1, self_regulation, vabundances_notnamed[sigma][cell], fullIM_list[sigma][cell], num_steps)
    end
end
end
############## simbio_spatial_map ################
# Read distance matrix
distance_matrix = CSV.read(joinpath(meta_path, "distance_matrix.csv"), DataFrame)[:,2:end]

function simbio_spatial_map!(K1, self_regulation, abundances, A_matrix, num_steps)
    num_species = length(abundances)
    
    # Check if lengths match
    if num_species != length(self_regulation) || num_species != length(K1)
        throw(ArgumentError("The number of species must be equal to the length of the self_regulation and K vectors"))
    end

    for t in 1:num_steps
        # Calculate growth rates
        growth_rates = growth.(abundances, self_regulation, K1)
        # println(growth_rates)
        # Calculate total change in abundances
        delta_abundances = growth_rates .+ sum(A_matrix .* self_regulation .* adjoint(abundances) .*abundances, dims=2)[:, 1]
        # println(delta_abundances[])
        # Emmigration

        # Update abundances in-place
        @. abundances += delta_abundances
        # println(abundances)
    end
    
    return abundances
end

@time Threads.@threads for sigma in 1:length(sigmas)
    for cell in 1:100
    # println("eso ", typeof(vabundances_notnamed[sigma][cell]))
    vupdated_abundances[sigma][cell] = simbio_spatial_map!(K1, self_regulation, vabundances_notnamed[sigma][cell], fullIM_list[sigma][cell], num_steps)
    end
end

############# CHECKING IF ALL ELEMENTS IN VUPDATED_ABUNDANCES ARE NON-ZERO ################
for i in 1:length(sigmas)
    for j in 1:length(vupdated_abundances[i])
        if any(vupdated_abundances[i][j] .!= 0.0)
            println("At least 1 element in vector vupdated_abundances[$i][$j] is non-zero.")
        end
    end
end

############## CHECKING IF VABUNDANCES AND VUPDATED_ABUNDANCES ARE THE SAME OBJECT ################
if vabundances === vupdated_abundances
    println("vabundances and vupdated_abundances are the same object.")
else
    println("vabundances and vupdated_abundances are not the same object.")
end


