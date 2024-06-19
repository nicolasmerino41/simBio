using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
# Pkg.add(path = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\Dispersal.jl")
# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
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
############################ CLEAN #################################
######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{256, T}
    b::T
    
    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{256, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{256, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.zero(::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))
# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)
# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*256)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*256)
# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyStructs256
function Base.sum(structs::MyStructs256...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs256 instance with the summed results
    return MyStructs256(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256, b::MyStructs256)
    return MyStructs256(max.(a.a, b.a))
end
# Define maximum for MyStructs256 with a scalar
function Base.maximum(a::MyStructs256, b::AbstractFloat)
    return MyStructs256(max.(a.a, b))
end
# Define maximum for a scalar with MyStructs256
function Base.maximum(a::AbstractFloat, b::MyStructs256)
    return MyStructs256(max.(a, b.a))
end
function Base.zeros(dims::NTuple{2, Int}, type = nothing)
    if type == MyStructs256{Float64}
        return [MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end
################# MYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    α::Float64
end

abstract type AbstractKernelNeighborhood end

struct CustomDispersalKernel{N<:DG.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DG.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

# Define neighbors for custom kernel
function DynamicGrids.neighbors(kernel::CustomDispersalKernel, hood, center::MyStructs256, I)
    result_a = zero(center.a)
    for i in 1:256
        for (j, neighbor) in enumerate(hood)
            if center.a[i] > 0.0
            dist = distance(I, hood.coords[j])
            result_a += kernel.formulation(dist) * neighbor.a[i]
            end
        end
    end
    return MyStructs256(result_a)
end
# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{Float64}}, kernel::SVector{9, Float64})
    
    result_a = SVector{256, Float64}(fill(0.0, 256))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    return MyStructs256(SVector(result_a))
end
function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyStructs256{Float64}}, kernel::SVector{25, Float64})
    
    result_a = SVector{256, Float64}(fill(0.0, 256))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    # println(sum(result_a))
    # result_b = sum(result_a)
    return MyStructs256(SVector(result_a)) #, result_b)
end
############ DESERIALIZING DATA ############################
############################################################
A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}
full_A_matrix = A_matrix_list[3]
full_IM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
full_IM = full_IM_list[2][10]
########### VISUAL OUTPUTS ################
###########################################
# Visualization settings
DynamicGrids.to_rgb(scheme::ObjectScheme, obj::MyStructs256) = ARGB32(
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0)
)
DynamicGrids.to_rgb(scheme, obj::MyStructs256) = get(scheme, clamp(obj.b, 0.0, 1.0))
##################### MWE_TvdD ##################################
####################################################################
####################################################################
####################################################################
# For a heat map we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs256, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(x -> x > 1, mystruct.a), A)
    return Makie.convert_arguments(t, richness)
end
##################### UTMtoRASTER ##################################
####################################################################
####################################################################
####################################################################
utmraster = Raster(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio\\Rasters\\updated_utmraster.tif")) 

species_df = CSV.File("DFs\\Species_spain_df.csv") |> DataFrame

variables = species_df[!, 2:5]
rename!(variables, [:ID, :Value, :sum, :UTMCODE])

species = species_df[!, unique_species_in_web]

species_df = hcat(variables, species, makeunique=true)
species_df_matrix = Matrix(species_df)

# Is this necessary?
for i in axes(species_df_matrix, 1), j in axes(species_df_matrix, 2)
    if species_df_matrix[i, j] == 0
        species_df_matrix[i, j] = 0.0
    elseif species_df_matrix[i, j] == 1
        species_df_matrix[i, j] = 1.0
    end
end

utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)

DA = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            DA[j] = MyStructs256(SVector{256, Float64}(species_df_matrix[i, 5:260]))
        end
    end
end

###### Let's do a small example of the simBio with the actual ATLAS data ######
# DA_with_abundances = deepcopy(DA)
# for row in axes(DA, 1), col in axes(DA, 2)
#     if DA[row, col] != mMyStructs256(Vector{Float64}(fill(0.0, 256)))
#         new_a = Vector{Float64}([DA[row, col].a[i] != 0.0 ? 100 : DA[row, col].a[i] for i in 1:256])
#         DA_with_abundances[row, col] = mMyStructs256(new_a)
#     end
# end
# serialize("Objects\\DA_with_abundances.jls", DA_with_abundances)
# serialize("Objects\\DA_with_abundances_all100.jls", DA_with_abundances)
# serialize("Objects\\DA_with_abundances_all100_mMyStructs256.jls", DA_with_abundances)
# Load the serialized object from the file and cast to the specific type
DA_with_abundances = deserialize("Objects\\DA_with_abundances_all100.jls")::DimArray{MyStructs256{Float64},2}

# This is for visualising richness in the raster or for creating a boolmask
# DA_sum = falses(dims(DA_with_abundances))
# for i in 1:size(species_df, 1)
#     # println(perro_cropped.Value[i])
#     for j in 1:125*76
        
#         # println(utmraster_da[j])
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             DA_sum[j] =  true #species_df_matrix[i, 3]
#         # else
#         #     DA_sum[j] = false
#         end
#     end
# end
# serialize("Objects\\DA_sum.jls", DA_sum)
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

# random_raster = [MyStructs256(SVector{256, Float64}(rand([0.0, 1.0], 256))) for i in 1:125, j in 1:76]
# random_richness_raster = [rand([0.0, rand(10:100)], 256) for i in 1:125, j in 1:76]
# idx = []
# for i in 1:125, j in 1:76
#     idx = push!(idx, CartesianIndex(i, j))
# end

function richness_evaluation(array_output, DA_with_presences)
    matches = 0
    for i in idx 
        if any(isnan, array_output[i].a)
            throw(ArgumentError("NaN found in array_output"))
        end
        above_ten = [x > 10 ? 1.0 : 0.0 for x in array_output[i].a]
        matches += sum(above_ten .== DA_with_presences[i])
    end
    return matches/(length(idx)*256)
end
######################### NPP ####################################
npp_absolute = CSV.File("DFs\\npp_absolute_df.csv") |> DataFrame
rename!(npp_absolute, [:ID, :UTMCODE, :npp, :X, :Y]) 
npp_absolute_in_kg = deepcopy(npp_absolute)
npp_absolute_in_kg.npp = npp_absolute.npp .* 1000
npp_absolute_in_kg = npp_absolute_in_kg[:, [2, 3]]

species_df = leftjoin(species_df, npp_absolute_in_kg, on = :UTMCODE)
species_df_matrix = Matrix(species_df)

# npp_DA = DimArray(Raster("Rasters/npp_utmsize_kgC.tif"), (Dim{:a}(1:125), Dim{:b}(1:76)))

# highest_npp = maximum(filter(!isnan, npp_DA[:]))
# lowest_npp = minimum(filter(!isnan, npp_DA[:]))
# lowest_richness = minimum(filter(!iszero, DA_richness[:]))

# npp_DA_r = reverse(npp_DA, dims=1)
# npp_DA_p = permutedims(npp_DA, (2, 1))
# MK.plot(npp_DA);
# # fixing some NA's that should not be in npp_DA
# for i in axes(npp_DA, 1), j in axes(npp_DA, 2)
#     if isnan(npp_DA[i, j]) && !iszero(DA_random_with_abundances[i, j])
#         println(i, j)
#         neighbors = get_neighbors(npp_DA, i, j)
#         neighbors = filter(!isnan, neighbors)
#         npp_DA[i, j] = mean(neighbors)
#     end 
# end
# serialize("Objects\\npp_DA.jls", npp_DA)
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

random_iberian_interact = deepcopy(iberian_interact_NA)
for row in axes(iberian_interact_NA, 1), col in axes(iberian_interact_NA, 2)
    if iberian_interact_NA[row, col] != 0.0
        random_iberian_interact[col, row] = -1.0
    end
end
# new = turn_adj_into_inter(iberian_interact_NA, 100, 0.5)

# Initialize an empty OrderedDict to hold the resulting matrices
results = OrderedDict{Float64, OrderedDict{Float64, Matrix}}()
self_regulation = 0.001
# Iterate over epsilon values
@time for epsilon in [0.1, 0.5, 1.0, 1.5, 3.0]
    # Initialize an OrderedDict for the current epsilon
    epsilon_results = OrderedDict{Float64, Matrix}()
    
    # Iterate over sigma values from 1.0 to 0.01, decrementing by 0.001 each time
    for sigma in 0.001
        caca = deepcopy(iberian_interact_NA)
        
        # Call the turn_adj_into_inter function with the current sigma and epsilon values
        result_matrix = turn_adj_into_inter(caca, sigma, epsilon)
        
        # Append the result to the epsilon_results OrderedDict
        epsilon_results[sigma] = result_matrix
    end
    
    # Store the epsilon_results OrderedDict in the main results OrderedDict
    results[epsilon] = epsilon_results
end
##################### Temp&Prec DA #################################
####################################################################
####################################################################
####################################################################
temp_raster = RS.Raster("Rasters\\utm_raster_temp.tif")
MK.plot(temp_raster);
prec_raster = RS.Raster("Rasters\\utm_raster_prec.tif")
MK.plot(prec_raster);

temp_DA = DimArray(temp_raster, (Dim{:a}(1:125), Dim{:b}(1:76)))
temp_DA = parent(temp_raster)
prec_DA = DimArray(prec_raster, (Dim{:a}(1:125), Dim{:b}(1:76)))
prec_DA = parent(prec_raster)

matrix_abundances = Matrix(DA_with_abundances)
x = 0
y = 0
for index in CartesianIndices(matrix_abundances)
    if isnan(prec_DA[index]) && !iszero(matrix_abundances[index])
        x += 1
        neighbors = get_neighbors(prec_DA, index[1], index[2])
        if !isempty(neighbors)
            prec_DA[index] = mean(neighbors)
        end
    end

    if isnan(temp_DA[index]) && !iszero(matrix_abundances[index])
        y += 1
        neighbors = get_neighbors(temp_DA, index[1], index[2])
        if !isempty(neighbors)
            temp_DA[index] = mean(neighbors)
        end
    end
end

println("Number of changes made: $x")
println("Number of changes made: $y")

##################### CLIMATE-ONLY MODEL ###########################
####################################################################
####################################################################
####################################################################
######## TEMPERATURE ##############
####### Species_temp
species_temp = CSV.File("DFs\\species_temp.csv") |> DataFrame
species_temp = species_temp[!, 2:4]
# Find the sorting indices based on matching names_order to species column
order_indices = indexin(names(iberian_interact_df), species_temp[:, :species])
# Reorder rows based on these indices
species_temp = species_temp[order_indices, :]
####### Species_temp_range
species_temp_range = CSV.File("DFs\\species_temp_range.csv") |> DataFrame
species_temp_range = species_temp_range[!, 2:4]
# Find the sorting indices based on matching names_order to species column
order_indices = indexin(names(iberian_interact_df), species_temp_range[:, :species])
# Reorder rows based on these indices
species_temp_range = species_temp_range[order_indices, :]
######## PRECIPITATION ##############
####### Species_prec
species_prec = CSV.File("DFs\\species_prec.csv") |> DataFrame
species_prec = species_prec[!, 2:4]
# Reorder rows based on these indices
species_prec = species_prec[order_indices, :]
####### Species_prec_range
species_prec_range = CSV.File("DFs\\species_prec_range.csv") |> DataFrame
species_prec_range = species_prec_range[!, 2:4]
# Reorder rows based on these indices
species_prec_range = species_prec_range[order_indices, :]

##################### NEW NICHES ###########################
######## bio rasters  ##############
# bio5_DA = DimArray(Raster("Rasters/bio5.tif"), (Dim{:a}(1:125), Dim{:b}(1:76)))
# bio6_DA = DimArray(Raster("Rasters/bio6.tif"), (Dim{:a}(1:125), Dim{:b}(1:76)))
# bio12_DA = DimArray(Raster("Rasters/bio12.tif"), (Dim{:a}(1:125), Dim{:b}(1:76)))

# # fixing some NA's that should not be in bioDA (change number for each bioclim)
# for i in axes(bio12_DA, 1), j in axes(bio12_DA, 2)
#     if isnan(bio12_DA[i, j]) && !iszero(DA_random_with_abundances[i, j])
#         println(i, j)
#         neighbors = get_neighbors(bio12_DA, i, j)
#         neighbors = filter(!isnan, neighbors)
#         bio12_DA[i, j] = mean(neighbors)
#     end 
# end
# serialize("Objects\\bio12.jls", bio12_DA)
bio5_DA = deserialize("Objects\\bio5.jls")
bio6_DA = deserialize("Objects\\bio6.jls")
bio12_DA = deserialize("Objects\\bio12.jls")
######## niches_df  ##############
species_niches = CSV.File("DFs\\iberian_species_niches.csv", decimal = ',') |> DataFrame
order_indices = indexin(names(iberian_interact_df), species_niches[:, :Species])
species_niches = species_niches[order_indices, :]

function int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(self_regulation * (npp+0.1) .*
    SVector{256}((1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5))) .*
    SVector{256}((1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6))) .*
    SVector{256}((1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12))) .*
    state.a .* (1.0 - (state.b / ((npp+0.1)))))
end
# function int_Gr_with_range(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, temp::AbstractFloat, prec::AbstractFloat)
#     return MyStructs256(self_regulation * (npp+0.1) .*
#     SVector{256}((1 ./ (1 .+ abs.(prec .- species_prec_range.mean_Prec) ./ species_prec_range.sd_Prec))) .*
#     SVector{256}((1 ./ (1 .+ abs.(temp .- species_temp_range.mean_Temp) ./ species_temp_range.sd_Temp))) .*
#     state.a .* (1.0 .- (state.a ./ ((npp+0.1)))))
# end
dimensions = (125, 76)
idx_tupled = [(i, j) for i in 1:dimensions[1], j in 1:dimensions[2]]
function random_dimarray(dimensions::Tuple{Int64, Int64}; prevalence = 0.1)
    init_array = DimArray(zeros(dimensions, MyStructs256{Float64}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
    for i in idx
        init_array[i] = MyStructs256(100.0 .* SVector{256}(sample([1,0], Weights([prevalence, 1-prevalence]), 256)))
    end
    return init_array
end
DA_random_with_abundances = random_dimarray(dimensions; prevalence = 0.1)
#TODO At prevalence 0.277 or higher you get instability
# function int_Gr(state::mMyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, temp::AbstractFloat, prec::AbstractFloat)
#     return mMyStructs256(self_regulation * npp .*
#     (1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
#     (1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
#     state.a .* (1.0 .- (state.a ./ (npp))))
# end
climatic_niche_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
    end
    # prec_factor = (1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec))
    # temp_factor = (1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp))
    # println(MyStructs256(SVector{256}((self_regulation * 1000.0) .* state.a .* (1.0 .- (state.a ./ 1000.0)))))  
    return state + int_Gr(state, self_regulation, npp, bio5, bio6, bio12) 
    # + MyStructs256(self_regulation * 100.0 .*
    #  SVector{256}(1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
    #  SVector{256}(1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
    #  state.a .* (1.0 .- (state.a ./ (100.0))))
end
function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end
# DA_with_abundances[18, 1] + MyStructs256(self_regulation * 100.0 .*
# SVector{256}(1 ./ (1 .+ abs.(30.0 .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
# SVector{256}(1 ./ (1 .+ abs.(20.0 .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
# DA_with_abundances[18, 1].a .* (1.0 .- (DA_with_abundances[18, 1].a ./ (100.0))))
indisp = InwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.1),
    distancemethod=AreaToArea(30)
);
outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.1),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

npp_DA = deserialize("Objects\\npp_DA.jls")
npp_DA = npp_DA./100000
highest_npp = maximum(filter(!isnan, npp_DA[:]))
highest_richness = maximum(filter(!iszero, DA_richness[:]))
pepe = ( 
    state = Matrix(DA_random_with_abundances),
    npp = Matrix(npp_DA),
    richness = Matrix(DA_richness),
    bio5 = Matrix(bio5_DA),
    bio6 = Matrix(bio6_DA),
    bio12 = Matrix(bio12_DA)
)

array_output = ArrayOutput(
    pepe; tspan = 1:500,
    mask = Matrix(DA_sum)
)
@time r = sim!(array_output, Ruleset(climatic_niche_rule, outdisp; boundary = Reflect()))
sum(r[100].state).b ≈ sum(r[1].state).b
richness_evaluation(r[length(r)].state, DA_with_presences)

makie_output = MakieOutput(pepe, tspan = 1:500;
    fps = 50, ruleset = Ruleset(climatic_niche_rule, outdisp; boundary = Reflect()),
    mask = Matrix(DA_sum)) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:state_heatmap, :state_image, :npp, :richness, :bio5, :bio6, :bio12]
    titles = ["Abundance", "Simulated Richness", "Carrying capacity", "Real Richness", "bio5", "bio6", "bio12"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*3 + j]) for i in 1:2, j in 1:3]
        
    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :state_heatmap
            Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=:seismic, colorrange=(0, 2000))
        elseif key == :state_image
            Makie.image!(ax, frame[:state], colormap=:seismic, colorrange=(0, 256))
        elseif key == :richness
            Makie.heatmap!(ax, frame[key]; interpolate=false, colormap=:seismic, colorrange=(0, 168))
        else
            Makie.heatmap!(ax, frame[key]; interpolate=false, colormap=:seismic)
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.yreversed[] = true
    end
end

hello = 0