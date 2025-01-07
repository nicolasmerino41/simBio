using Pkg

Pkg.activate(joinpath(pwd(), "/simBio"))
cd(joinpath(pwd(), "/simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")

# Packages
using ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase
using ColorSchemes 
using Makie, WGLMakie

const DG, MK, AG, RS, Disp, DF = DynamicGrids, Makie, ArchGDAL, Rasters, Dispersal, DataFrames
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]

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
beta = 3.0

gbif_sizes = CSV.read("DFs\\gbif_sizes.csv", DataFrame)[:, 2:end]

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

web = CSV.read("DFs\\TetraEU_pairwise_interactions.csv", DataFrame)

web = DataFrame(predator = web.sourceTaxonName, prey = web.targetTaxonName)

web.predator = string.(web.predator)
web.prey = string.(web.prey)
web = web[:, [:predator, :prey]]

unique_predators = unique(web.predator)
unique_preys = unique(web.prey)

x = vcat(unique_predators, unique_preys)
unique_species = unique(x)

# Read the CSV file
diets = CSV.File("DFs\\TetraEU_generic_diet.csv") |> DataFrame

diets = hcat(diets.sourceTaxonName, diets.targetGenericItemName)

Amph = CSV.read("DFs/s", delim='\t', DataFrame)
Bird = CSV.read("DFs/DB_Birds_IP.txt", delim='\t', DataFrame)
Mamm = CSV.read("DFs/DB_Mammals_IP.txt", delim='\t', DataFrame)
Rept = CSV.read("DFs/DB_Reptiles_IP.txt", delim='\t', DataFrame)

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
function Base.maximum(a::MyStructs256)
    return max.(a.b)
end
function Base.maximum(a::Matrix{MyStructs256{Float64}})
    # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
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

##################### UTMtoRASTER ##################################
####################################################################
####################################################################
DA_with_abundances = deserialize("Objects\\DA_with_abundances_all100.jls")::DimArray{MyStructs256{Float64},2}

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
@time for epsilon in [1.0]
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
full_IM = results[1][1]
##################### CLIMATE-ONLY MODEL ###########################
####################################################################
####################################################################
####################################################################
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
order_indices = indexin(names(iberian_interact_df), species_niches[:, :Species])
species_niches = species_niches[order_indices, :]

lax_species_niches = CSV.File("DFs\\iberian_species_niches_withLaxNiche.csv", decimal = ',') |> DataFrame
order_indices = indexin(names(iberian_interact_df), lax_species_niches[:, :Species])
lax_species_niches = lax_species_niches[order_indices, :]

strict_species_niches = CSV.File("DFs\\iberian_species_niches_withVeryStrictNiche.csv", decimal = ',') |> DataFrame
order_indices = indexin(names(iberian_interact_df), strict_species_niches[:, :Species])
strict_species_niches = strict_species_niches[order_indices, :]

herbivore_names = CSV.File(joinpath(meta_path, "herbivore_names.csv")) |> DataFrame
herbivore_names = convert(Vector{String}, herbivore_names[:, 2])
binary_vector = [name in herbivore_names ? 1 : 0 for name in names(iberian_interact_df)]
opposite_binary_vector = [name in herbivore_names ? 0 : 1 for name in names(iberian_interact_df)]

function lax_int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(self_regulation * (npp+0.1) .*
    SVector{256}((1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5))) .*
    SVector{256}((1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6))) .*
    SVector{256}((1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12))) .*
    state.a .* (1.0 - (state.b / ((npp+0.1)))))
end
function int_Gr_for_biotic(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(SVector{256}(self_regulation * (npp .+ 0.1) .* 
    SVector{256}(1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5)) .* 
    SVector{256}(1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6)) .* 
    SVector{256}(1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12)) .* 
    state.a .* (1.0 .- (state.a ./ (npp .* binary_vector .+ npp .* opposite_binary_vector ./ 100.0)))))
end
function trophic_optimized(abundances, A_matrix)
    # Calculate the weighted interaction directly
    interaction = A_matrix * abundances.a
    return MyStructs256(SVector(interaction .* abundances.a))
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

biotic_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    merged_state = state + 
        int_Gr_for_biotic(state, self_regulation, npp, bio5, bio6, bio12)  +
        trophic_optimized(state, full_IM_used)
    return MyStructs256(max.(0.0, merged_state.a))
end
function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.2),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

npp_DA = deserialize("Objects\\npp_DA.jls")
npp_DA = npp_DA./1000000

function random_dimarray(dimensions::Tuple{Int64, Int64}; prevalence = 0.1)
    init_array = DimArray(zeros(dimensions, MyStructs256{Float64}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
    for i in idx
        init_array[i] = MyStructs256(SVector{256}(10.0 .* binary_vector .+ 10.0 .* opposite_binary_vector) .* SVector{256}(sample([1,0], Weights([prevalence, 1-prevalence]), 256)))
    end
    return init_array
end

DA_random_with_abundances = random_dimarray(dimensions; prevalence = 0.1)
#TODO At prevalence 0.277 or higher you get instability
pepe = ( 
    state = Matrix(DA_random_with_abundances),
    npp = Matrix(npp_DA),
    richness = Matrix(DA_richness),
    bio5 = Matrix(bio5_DA),
    bio6 = Matrix(bio6_DA),
    bio12 = Matrix(bio12_DA)
)

array_output = ArrayOutput(
    pepe; tspan = 1:100,
    mask = Matrix(DA_sum)
)
@time r = sim!(array_output, Ruleset(lax_climatic_niche_rule, outdisp; boundary = Wrap()))
rich_eval = richness_evaluation(r[length(r)].state, DA_with_presences)
