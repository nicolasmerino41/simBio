
using Pkg
PC = "nicol"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")

################################################################################################
####################### REAL SIMULATION ############################
####################################################################
####################################################################
####################################################################
######################## Size_selection kernel ########################
#######################################################################
function size_selection_kernel(predator_mass, prey_mass, sd, beta)
    intensity = exp(-((log(float(predator_mass)) / (beta * float(prey_mass)))^2.0) / (2.0 * float(sd)^2.0))
    return float(intensity)
end
beta = float(3)

gbif_sizes = CSV.read(joinpath(meta_path, "Lists\\gbif_sizes.csv"), DataFrame)[:, 2:end]
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
#################### get_neighbors ###################################
######################################################################
# Helper function to get the neighbors
function get_neighbors(matrix, row, col)
    neighbors = []
    rows, cols = size(matrix)
    
    for r in max(row-1, 1):min(row+1, rows)
        for c in max(col-1, 1):min(col+1, cols)
            if (r != row || c != col) && !isnan(matrix[r, c])
                push!(neighbors, matrix[r, c])
            end
        end
    end
    return neighbors
end
#################### when_NA ###################################
######################################################################
function when_NA(array_output)
    for time in 1:length(array_output)
        value = 0.0
        for index in idx
            if any(isinf, array_output[time].state[index].a) || any(isnan, array_output[time].state[index].a)
              println("Time:", time, " Index:", index)
              value += 1
            end
        end
        if value != 0
            println("In time ", time,", ", value, " NA's where generated for the first time")
            return
        end
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

Amph = CSV.read(joinpath(meta_path, "Atlas_data/DB_Amphibians_IP.txt"), delim='\t', DataFrame)
Bird = CSV.read(joinpath(meta_path, "Atlas_data/DB_Birds_IP.txt"), delim='\t', DataFrame)
Mamm = CSV.read(joinpath(meta_path, "Atlas_data/DB_Mammals_IP.txt"), delim='\t', DataFrame)
Rept = CSV.read(joinpath(meta_path, "Atlas_data/DB_Reptiles_IP.txt"), delim='\t', DataFrame)

amphibian_names = names(Amph)[2:end]
reptile_names = names(Rept)[2:end]
mammal_names = names(Mamm)[2:end]
bird_names = names(Bird)[2:end]
herps_names = append!(deepcopy(amphibian_names), deepcopy( reptile_names))
birmmals_names = append!(deepcopy(mammal_names), deepcopy(bird_names))
spain_fauna = append!(deepcopy(herps_names), deepcopy(birmmals_names)) 

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
# Ordering the matrix by Apmhibians, Reptiles, Mammals, Birds
spain_names = filter(name -> name in unique_species_in_web, names(hcat(Amph[:, Not(:UTMCODE)], Rept[:, Not(:UTMCODE)], Mamm[:, Not(:UTMCODE)], Bird[:, Not(:UTMCODE)])))
iberian_interact_matrix = iberian_interact_matrix[:, spain_names]
iberian_interact_matrix = iberian_interact_matrix[spain_names, :]

## Creating a mapping from species names to matrix indices
species_to_index = Dict(zip(spain_names, 1:n))
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
self_regulation = 0.001
iberian_interact_df = DataFrame(iberian_interact_matrix, species_names)
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

##################### UTMtoRASTER ##################################
####################################################################
####################################################################
####################################################################
species_df = CSV.File("DFs\\Species_spain_df.csv") |> DataFrame

variables = species_df[!, 2:5]
rename!(variables, [:ID, :Value, :sum, :UTMCODE])

species = species_df[!, unique_species_in_web]
species = species[:, spain_names]

species_df = hcat(variables, species, makeunique=true)
species_df_matrix = Matrix(species_df)

utmraster = Raster("Rasters\\updated_utmraster.tif") 
utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)

DA = deserialize("Objects\\DA.jls")
DA_herps = deserialize("Objects\\DA_herps.jls")
DA_birmmals = deserialize("Objects\\DA_birmmals.jls")

DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if DA[row, col] != MyStructs256(Vector{Float64}(fill(0.0, 256)))
        new_a = Vector{Float64}([DA[row, col].a[i] != 0.0 ? 10.0 : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = MyStructs256(new_a)
    end
end
DA_birmmals_with_abundances = deepcopy(DA_birmmals)
# Iterate over rows and columns
for row in axes(DA_birmmals, 1), col in axes(DA_birmmals, 2)
    current = DA_birmmals[row, col]
    empty_birmmals = MyBirmmals(fill(0.0, 207))
    
    if current != empty_birmmals
        new_a = [current.a[i] != 0.0 ? 0.01 : current.a[i] for i in 1:207]
        DA_birmmals_with_abundances[row, col] = MyBirmmals(new_a)
    end
end
DA_herps_with_abundances = deepcopy(DA_herps)
# Iterate over rows and columns
for row in axes(DA, 1), col in axes(DA, 2)
    current = DA[row, col]
    empty_struct = MyHerps(fill(0.0, 49))
    
    if current != empty_struct
        new_a = [current.a[i] != 0.0 ? 0.01 : current.a[i] for i in 1:49]
        DA_herps_with_abundances[row, col] = MyHerps(new_a)
    end
end

DA_sum = deserialize("Objects\\DA_sum.jls")
DA_sum_r = reverse(DA_sum, dims=1)
DA_sum_p = permutedims(DA_sum, (2, 1))

DA_with_abundances_r = reverse(DA_with_abundances, dims=1)
DA_with_abundances_p = permutedims(DA_with_abundances, (2, 1))
DA_with_abundances_p_masked = deepcopy(DA_with_abundances_p)

DA_richness = deserialize("Objects\\DA_richness.jls")::DimArray{Float64,2}

########################## IDX #####################################
idx = findall(x -> x == 1.0, DA_sum)
DA_with_presences = DimArray([fill(0.0, 256) for _ in 1:125, _ in 1:76], (Dim{:a}(1:125), Dim{:b}(1:76)))

for row in axes(DA_with_abundances, 1), col in axes(DA_with_abundances, 2)
    if DA_with_abundances[row, col].b != 0.0
        for i in findall(!iszero, DA_with_abundances[row, col].a)
            DA_with_presences[row, col][i] = 1.0
        end
    end
end

######################### NPP ####################################
npp_absolute = CSV.File("DFs\\npp_absolute_df.csv") |> DataFrame
rename!(npp_absolute, [:ID, :UTMCODE, :npp, :X, :Y]) 
npp_absolute_in_kg = deepcopy(npp_absolute)
npp_absolute_in_kg.npp = npp_absolute.npp .* 1000
npp_absolute_in_kg = npp_absolute_in_kg[:, [2, 3]]

species_df = leftjoin(species_df, npp_absolute_in_kg, on = :UTMCODE, makeunique = true)
species_df_matrix = Matrix(species_df)

npp_DA = deserialize("Objects\\npp_DA.jls")
npp_DA = npp_DA./10000
################### EFFICIENT MATRIX FRAMEWORK #####################
####################################################################
####################################################################
####################################################################
# Load a DataFrame from a serialized file ('.jls' format).
iberian_interact_df = deserialize("Objects\\iberian_interact_df.jls")
# Convert the DataFrame to a matrix for easier manipulation.
iberian_interact_matrix = iberian_interact_df |> Matrix
# Convert the modified matrix back to a DataFrame, preserving the original column names.
iberian_interact_df = DataFrame(iberian_interact_matrix, names(iberian_interact_df))
# Create a NamedArray from the matrix, using the DataFrame's column names for both dimensions.
iberian_interact_NA = NamedArray(
    iberian_interact_matrix, 
    (names(iberian_interact_df), names(iberian_interact_df)),
    ("Species", "Species")
)
iberian_interact_NA = iberian_interact_NA[spain_names, spain_names]


##################### NEW NICHES ###########################
######## bio rasters  ##############
bio5_DA = deserialize("Objects\\bio5.jls")
bio6_DA = deserialize("Objects\\bio6.jls")
bio12_DA = deserialize("Objects\\bio12.jls")
futurebio5_DA = bio5_DA .+ 1.0
futurebio6_DA = bio6_DA .+ 1.0
futurebio12_DA = bio12_DA .+ rand(Normal(0, 100), 125, 76)

######## niches_df  ##############
species_niches = CSV.File("DFs\\iberian_species_niches_withbinned_TH.csv", decimal = ',') |> DataFrame
order_indices = indexin(spain_names, species_niches[:, :Species])
species_niches = species_niches[order_indices, :]

lax_species_niches = CSV.File("DFs\\iberian_species_niches_withLaxNiche.csv", decimal = ',') |> DataFrame
order_indices = indexin(spain_names, lax_species_niches[:, :Species])
lax_species_niches = lax_species_niches[order_indices, :]

strict_species_niches = CSV.File("DFs\\iberian_species_niches_withVeryStrictNiche.csv", decimal = ',') |> DataFrame
order_indices = indexin(spain_names, strict_species_niches[:, :Species])
strict_species_niches = strict_species_niches[order_indices, :]

herbivore_names = CSV.File(joinpath(meta_path, "herbivore_names.csv")) |> DataFrame
herbivore_names = convert(Vector{String}, herbivore_names[:, 2])
binary_vector = [name in herbivore_names ? 1 : 0 for name in names(iberian_interact_df)]
opposite_binary_vector = [name in herbivore_names ? 0 : 1 for name in names(iberian_interact_df)]

function int_Gr_for_biotic(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(self_regulation * (npp .+ 0.1) .* 
    1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5) .* 
    1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6) .* 
    1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12) .* 
    state.a .* (1.0 .- (state.a ./ (npp .* binary_vector .+ npp .* opposite_binary_vector ./ 100.0))))
end

function int_Gr_for_biotic_k(state::MyStructs256, self_regulation::AbstractFloat, k_DA::MyStructs256)
    return MyStructs256(self_regulation .* (k_DA.a .+ 0.00001) .* 
    state.a .* (1.0 .- (state.a ./ (k_DA.a .* binary_vector .+ k_DA.a .* opposite_binary_vector))))
end

function lax_int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(self_regulation * (npp+0.1) .*
    (1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)) .*
    (1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)) .*
    (1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)) .*
    state.a .* (1.0 - (state.b / ((npp+0.1)))))
end

function trophic_optimized(abundances, full_IM)
    # Calculate the weighted interaction directly
    interaction = full_IM * abundances.a
    return MyStructs256(interaction .* abundances.a)
end
###################### GENERATE Ki(z) DimArray ################
###############################################################
k_DA = DimArray(reshape([MyStructs256(fill(0.0, 256)) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

for row in axes(k_DA, 1), col in axes(k_DA, 2)
    if isone(DA_sum[row, col])
    
        k_DA[row, col] = MyStructs256(1 ./ (1 .+ abs.(bio5_DA[row, col] .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5) .*
                                        1 ./ (1 .+ abs.(bio6_DA[row, col] .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6) .*
                                        1 ./ (1 .+ abs.(bio12_DA[row, col] .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12))
    end
end

dimensions = (125, 76)
idx_tupled = [(i, j) for i in 1:dimensions[1], j in 1:dimensions[2]]

lax_climatic_niche_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    return state + lax_int_Gr(state, self_regulation, npp, bio5, bio6, bio12)
end

biotic_rule_k = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    merged_state = state + 
        int_Gr_for_biotic_k(state, self_regulation, k_DA)  +
        trophic_optimized(state, full_IM)
    return MyStructs256(max.(0.0, merged_state.a))
end
biotic_rule_k_herps = Cell{Tuple{:herps, :birmmals, :k_DA}, :herps}() do data, (herps, birmmals, k_DA), I
    if any(isinf, birmmals.a) || any(isnan, birmmals.a)
        @error "state has NA values in birmmals"
        println(I)
    end
    if any(isinf, herps.a) || any(isnan, herps.a)
        @error "state has NA values in herps"
        println(I)
    end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyHerps(max.(0.0, merged_state.a[1:49]))
end
biotic_rule_k_birmmals = Cell{Tuple{:herps, :birmmals, :k_DA}, :birmmals}() do data, (herps, birmmals, k_DA), I
    if typeof(birmmals) != MyBirmmals{Float64}
        @error "birmmals is not MyBirmmals"
    end
    if typeof(herps) != MyHerps{Float64}
        @error "herps is not MyHerps"
    end
    if any(isinf, birmmals.a) || any(isnan, birmmals.a)
        @error "state has NA values in birmmals"
        println(I)
    end
    if any(isinf, herps.a) || any(isnan, herps.a)
        @error "state has NA values in herps"
        println(I)
    end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyBirmmals(max.(0.0, merged_state.a[50:256]))
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end