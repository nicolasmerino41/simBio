#################################################################################################
########################### LET'S GO ###############################
####################################################################
####################################################################
####################################################################
###################### remove_variables ############################
function remove_variable(var_name::Symbol)
    if isdefined(Main, var_name)
        eval(:(global $var_name = nothing))
        println("Variable $var_name removed from the environment.")
    else
        println("Variable $var_name is not defined.")
    end
end
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
#################### abundance_over_time##############################
######################################################################
function abundance_over_time(abundances)
    # Convert data to a matrix and transpose it
    transposed_abundances = hcat(abundances...)
    # Define timesteps
    timesteps = 1:size(transposed_abundances, 2)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Abundance Over Time", xlabel = "Timesteps", ylabel = "Abundance") 
    # Plot each individual's abundance
    for i in 1:size(transposed_abundances, 1)
        lines!(ax, timesteps, transposed_abundances[i, :], label = "Individual $i") 
    end
    display(fig)
end
#################### count_zeros_ones ##############################
######################################################################
function count_zeros_ones(DA::DimArray{Vector{Float64}, 2}, idx::Vector{CartesianIndex{2}})
    zero_count = 0
    one_count = 0
    
    for index in idx
        current_cell = DA[index]
        zero_count += count(x -> x == 0.0, current_cell)
        one_count += count(x -> x == 1.0, current_cell)
    end
    
    return zero_count, one_count
end

################### map_plot ######################
###################################################
function map_plot(plot::AbstractArray; type = nothing, lambda_DA = nothing, palette = nothing, legend = false, show_grid = true, flip = true, title = nothing, kw...)
    
    if isa(plot, DimArray)
        plot = Matrix(plot)
    end
    
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set the colormap
    pal = isnothing(palette) ? :thermal : palette

    # Initialize the variable for the plotting object
    plt_obj = nothing

    # Plot according to the specified type
    if type == "image"
        if isnothing(lambda_DA)
            @error("No lambda_DA being used, choose one.")
        end
        plt_obj = image!(ax, plot, lambda_DA; colormap = pal, kw...)
    elseif type == "plot"
        plt_obj = plot!(ax, plot; kw...)
    else
        plt_obj = heatmap!(ax, plot; colormap = pal, kw...)
    end

    # Reverse the y-axis if needed
    ax.yreversed = flip

    # Add legend if requested
    if legend
        non_na_values = filter(!isnan, plot)
        # Ensure that non_na_values is not empty to prevent errors
        if !isempty(non_na_values)
            # Remove colormap and limits from the Colorbar call
            Colorbar(fig[1, 2], plt_obj)
        else
            @warn "No valid data for Colorbar."
        end
    end

    # Add title if provided
    if !isnothing(title)
        ax.title = title
    end

    hidexdecorations!(ax; grid = show_grid)
    hideydecorations!(ax; grid = show_grid)

    fig
end

# Define the custom colormap
custom_palette = cgrad([colorant"black", colorant"blue", colorant"yellow", colorant"green", colorant"red"], [0.0, 0.000000000001, 0.33, 0.66, 1.0]);
# map_plot(DA_richness; palette = custom_palette, rev = true, type = "plot")

####################################################################################################
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
amphibian_names = names(Amph)[2:end]
reptile_names = names(Rept)[2:end]
mammal_names = names(Mamm)[2:end]
bird_names = names(Bird)[2:end]
herps_names = append!(deepcopy(amphibian_names), deepcopy( reptile_names))
birmmals_names = append!(deepcopy(mammal_names), deepcopy(bird_names))
spain_fauna = append!(deepcopy(herps_names), deepcopy(birmmals_names)) 
# bird_names_in_unique_species = length(intersect(bird_names, unique_species_in_web))
# println("The number of bird names found in unique_species_in_web is $bird_names_in_unique_species")
# mammal_names_in_unique_species = length(intersect(mammal_names, unique_species_in_web))
# println("The number of mammal names found in unique_species_in_web is $mammal_names_in_unique_species")
# reptile_names_in_unique_species = length(intersect(reptile_names, unique_species_in_web))
# println("The number of reptile names found in unique_species_in_web is $reptile_names_in_unique_species")
# amphibian_names_in_unique_species = length(intersect(amphibian_names, unique_species_in_web))
# println("The number of amphibian names found in unique_species_in_web is $amphibian_names_in_unique_species")

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
iberian_interact_df = DataFrame(iberian_interact_matrix, species_names)
function turn_adj_into_inter(adjacencyy, sigma, epsilon, self_regulation, beta)
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

############ DESERIALIZING DATA ############################
############################################################
A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}
full_A_matrix = A_matrix_list[3]
full_IM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
# full_IM = full_IM_list[2][10]
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
# # For a heatmap we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyBirmmals, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyHerps, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs256, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyBirmmals, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_birds[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyHerps, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_herps[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
# # WITH LAMBDA
# For MyStructs
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs256, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# For MyBirmmals
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyBirmmals, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector_birds[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# For MyHerps
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyHerps, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector_herps[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# FOR RASTER
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyBirmmals, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyHerps, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyStructs256, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyBirmmals, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_birds[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyHerps, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_herps[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
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
species_df.sum = [sum(species_df[i, 5:260]) for i in 1:size(species_df, 1)]
species_df_matrix = Matrix(species_df)

utmraster = Raster("Rasters\\updated_utmraster.tif")
utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)

# Initialize DA with MyStructs256 of zeroes
# DA = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
# # Iterate over species_df and populate DA
# for i in 1:size(species_df, 1)
#     for j in 1:125*76
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             if species_df.sum[i] == 54 
#                println(utmraster_DA[3]) 
#             end
#             # Convert species_df_matrix[i, 5:260] to SVector{256, Float64} before creating MyStructs256
#             DA[j] = MyStructs256(SVector{256, Float64}(species_df_matrix[i, 5:260]))
#         end
#     end
# end
# serialize("Objects\\DA.jls", DA)
DA = deserialize("Objects\\DA.jls")
# Save DA to a .jld2 file in the Objects1_9 folder
# @save "Objects1_9/DA.jld2" DA
# @load "Objects1_9/DA.jld2" DA
# # Initialize DA_herps with MyHerps instances filled with zeros using SVector
# DA_herps = DimArray(reshape([MyHerps(SVector{49, Float64}(fill(0.0, 49))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
# # Iterate over species_df and populate DA_herps
# for i in 1:size(species_df, 1)
#     for j in 1:125*76
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             # Assign the first 49 elements of DA[j].a to DA_herps[j] using SVector
#             DA_herps[j] = MyHerps(SVector{49, Float64}(DA[j].a[1:49]))
#         end
#     end
# end
# # Serialize the updated DA_herps object
# serialize("Objects\\DA_herps.jls", DA_herps)
DA_herps = deserialize("Objects\\DA_herps.jls")
# @load "Objects1_9/DA_herps.jld2" DA_herps
# # Initialize DA_birmmals with MyBirmmals instances filled with zeros using SVector
# DA_birmmals = DimArray(reshape([MyBirmmals(SVector{207, Float64}(fill(0.0, 207))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
# # Iterate over species_df and populate DA_birmmals
# for i in 1:size(species_df, 1)
#     for j in 1:125*76
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             # Assign the elements 50:256 of DA[j].a to DA_birmmals[j] using SVector
#             DA_birmmals[j] = MyBirmmals(SVector{207, Float64}(DA[j].a[50:256]))
#         end
#     end
# end
# # Serialize the updated DA_birmmals object
# serialize("Objects\\DA_birmmals.jls", DA_birmmals)
DA_birmmals = deserialize("Objects\\DA_birmmals.jls")
# @load "Objects1_9/DA_birmmals.jld2" DA_birmmals
##### Let's do a small example of the simBio with the actual ATLAS data ######
DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if DA[row, col] != MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        new_a = SVector{256, Float64}([DA[row, col].a[i] != 0.0 ? 10.0 : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = MyStructs256(new_a)
    end
end

initial_abundance = 0.41
DA_birmmals_with_abundances = deepcopy(DA_birmmals)
# Iterate over rows and columns
for row in axes(DA_birmmals, 1), col in axes(DA_birmmals, 2)
    current = DA_birmmals[row, col]
    empty_birmmals = MyBirmmals(SVector{207, Float64}(fill(0.0, 207)))
    
    if current != empty_birmmals
        new_a = SVector{207, Float64}([current.a[i] != 0.0 ? initial_abundance : current.a[i] for i in 1:207])
        DA_birmmals_with_abundances[row, col] = MyBirmmals(new_a)
    end
end

DA_herps_with_abundances = deepcopy(DA_herps)
# Iterate over rows and columns
for row in axes(DA, 1), col in axes(DA, 2)
    current = DA[row, col]
    empty_struct = MyHerps(SVector{49, Float64}(fill(0.0, 49)))
    
    if current != empty_struct
        new_a = SVector{49, Float64}([current.a[i] != 0.0 ? initial_abundance : current.a[i] for i in 1:49])
        DA_herps_with_abundances[row, col] = MyHerps(new_a)
    end
end
# serialize("Objects\\DA_with_abundances.jls", DA_with_abundances)
# serialize("Objects\\DA_with_abundances_all10.jls", DA_with_abundances)
# serialize("Objects\\DA_with_abundances_all100_mMyStructs256.jls", DA_with_abundances)
# Load the serialized object from the file and cast to the specific type
# DA_with_abundances = deserialize("Objects\\DA_with_abundances_all10.jls")::DimArray{MyStructs256{Float64},2}

# This is for visualising richness in the raster or for creating a boolmask
# DA_richness = zeros(dims(DA_with_abundances))
# for i in 1:size(species_df, 1)
#     # println(perro_cropped.Value[i])
#     for j in 1:125*76
        
#         # println(utmraster_da[j])
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             DA_richness[j] =  species_df_matrix[i, 3] #true
#         # else
#         #     DA_sum[j] = false
#         end
#     end
# end
# serialize("Objects\\DA_sum.jls", DA_sum)
DA_sum = deserialize("Objects\\DA_sum.jls")
# @load "Objects1_9/DA_sum.jld2" DA_sum
DA_sum_r = reverse(DA_sum, dims=1)
DA_sum_p = permutedims(DA_sum, (2, 1))

DA_with_abundances_r = reverse(DA_with_abundances, dims=1)
DA_with_abundances_p = permutedims(DA_with_abundances, (2, 1))
DA_with_abundances_p_masked = deepcopy(DA_with_abundances_p)
# raster_with_abundances_r = reverse(raster_with_abundances, dims=1)
# raster_with_abundances_p = permutedims(raster_with_abundances, (2, 1))
# DA_richness[34, 9] = Int(round(mean([105, 105, 114, 97, 73, 106, 84])))
# DA_richness[34, 10] = Int(round(mean([73, 97, 106, 84, 69, 124])))
# DA_richness[34, 11] = Int(round(mean([84, 69, 76, 106, 124, 86])))
# DA_richness[34, 12] = Int(round(mean([88, 90, 86, 124, 69, 76])))
# DA_richness[34, 13] = Int(round(mean([102, 94, 91, 90, 86, 88, 76])))
# serialize("Objects\\DA_richness.jls", DA_richness)

DA_richness = deserialize("Objects\\DA_richness.jls")::DimArray{Float64,2}
@load "Objects1_9/DA_richness.jld2" DA_richness
DA_richness[DA_richness .== 0.0] .= NaN
DA_richness_birmmals = deserialize("Objects/DA_richness_birmmals.jls")::DimArray{Float64,2}
@load "Objects1_9/DA_richness_birmmals.jld2" DA_richness_birmmals
DA_richness_birmmals[DA_richness_birmmals .== 0.0] .= NaN
DA_richness_herps = deserialize("Objects/DA_richness_herps.jls")::DimArray{Float64,2}
@load "Objects1_9/DA_richness_herps.jld2" DA_richness_herps
DA_richness_herps[DA_richness_herps .== 0.0] .= NaN

######################## RASTERISING DAs ################################
# See Rasterising DAs.jl for more details
raster_DA = deserialize("Objects\\raster_DA.jls")
raster_herps = deserialize("Objects\\raster_herps.jls")
raster_birmmals = deserialize("Objects\\raster_birmmals.jls")
raster_with_abundances = deserialize("Objects\\raster_with_abundances.jls")
raster_birmmals_with_abundances = deserialize("Objects\\raster_birmmals_with_abundances.jls")
raster_herps_with_abundances = deserialize("Objects\\raster_herps_with_abundances.jls")
raster_richness = deserialize("Objects\\raster_richness.jls")
raster_richness_birmmals = deserialize("Objects\\raster_richness_birmmals.jls")
raster_richness_herps = deserialize("Objects\\raster_richness_herps.jls")
raster_sum = deserialize("Objects\\raster_sum.jls")

########################## IDX #####################################
idx = findall(x -> x == 1.0, DA_sum)
DA_with_presences = DimArray([fill(0.0, 256) for _ in 1:125, _ in 1:76], (Dim{:a}(1:125), Dim{:b}(1:76)))
sum(DA_with_presences[18,1])
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
threshold = 0.061688
function presence_absence_prediction_accuracy(array_output, DA_with_presences, threshold)
    true_positives = 0
    true_negatives = 0
    total_presences = 0
    total_absences = 0
    
    for i in idx
        current_cell = array_output[i].a
        predicted_presence = current_cell .> threshold
        
        actual_presence = DA_with_presences[i]
        
        true_positives += sum((predicted_presence .== 1) .& (actual_presence .== 1))
        true_negatives += sum((predicted_presence .== 0) .& (actual_presence .== 0))
        
        total_presences += sum(actual_presence .== 1)
        total_absences += sum(actual_presence .== 0)
    end
    
    tpr = true_positives / total_presences
    tnr = true_negatives / total_absences
    
    if total_presences == 0 || total_absences == 0
        return NaN  # Return NaN if there are no presences or absences
    else
        return 2 * (tpr * tnr) / (tpr + tnr)  # Harmonic mean of TPR and TNR
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
# @load "Objects1_9/npp_DA.jld2" npp_DA
npp_DA = deserialize("Objects\\npp_DA.jls")
npp_raster = deepcopy(raster_richness)
# [npp_raster[i, j] = npp_DA[i, j] for i in axes(npp_raster, 1), j in axes(npp_raster, 2)]
# serialize("Objects\\npp_raster.jls", npp_raster)
npp_raster = deserialize("Objects\\npp_raster.jls")
################### EFFICIENT MATRIX FRAMEWORK #####################
####################################################################
####################################################################
####################################################################
# Load a DataFrame from a serialized file ('.jls' format).
iberian_interact_df = deserialize("Objects\\iberian_interact_df.jls")
iberian_interact_df = CSV.File("Objects1_9/iberian_interact_df.csv") |> DataFrame
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
# iberian_interact_NA_herps = deepcopy(iberian_interact_NA)
# iberian_interact_NA_birmmals = deepcopy(iberian_interact_NA)

# iberian_interact_NA_herps = iberian_interact_NA_herps[filter(row -> row ∈ herps_names, names(iberian_interact_NA_herps, 1)),
#  filter(col -> col ∈ herps_names, names(iberian_interact_NA_herps, 2))]
# iberian_interact_NA_birmmals = iberian_interact_NA_birmmals[filter(row -> row ∈ birmmals_names, names(iberian_interact_NA_birmmals, 1)),
#  filter(row -> row ∈ birmmals_names, names(iberian_interact_NA_birmmals, 2))]

# Initialize an empty OrderedDict to hold the resulting matrices
results = OrderedDict{Float64, OrderedDict{Float64, Matrix}}()
self_regulation = 1.0
# # Iterate over epsilon values
# @time for epsilon in [0.1, 0.5, 1.0, 1.5, 2.031]
#     # Initialize an OrderedDict for the current epsilon
#     epsilon_results = OrderedDict{Float64, Matrix}()
    
#     # Iterate over sigma values from 1.0 to 0.01, decrementing by 0.001 each time
#     for sigma in 0.1312
#         caca = deepcopy(iberian_interact_NA)
        
#         # Call the turn_adj_into_inter function with the current sigma and epsilon values
#         result_matrix = turn_adj_into_inter(caca, sigma, epsilon)
        
#         # Append the result to the epsilon_results OrderedDict
#         epsilon_results[sigma] = result_matrix
#     end
    
#     # Store the epsilon_results OrderedDict in the main results OrderedDict
#     results[epsilon] = epsilon_results
# end
# full_IM = results[2.031][0.1312]

# caca = deepcopy(iberian_interact_NA)
# sigma = 0.5
# epsilon = 1.0
# full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon))
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
# coco = 0
# cocolo = 0
for index in CartesianIndices(matrix_abundances)
    if isnan(prec_DA[index]) && !iszero(matrix_abundances[index])
        # coco += 1
        neighbors = get_neighbors(prec_DA, index[1], index[2])
        if !isempty(neighbors)
            prec_DA[index] = mean(neighbors)
        end
    end

    if isnan(temp_DA[index]) && !iszero(matrix_abundances[index])
        # cocolo += 1
        neighbors = get_neighbors(temp_DA, index[1], index[2])
        if !isempty(neighbors)
            temp_DA[index] = mean(neighbors)
        end
    end
end

# println("Number of changes made: $coco")
# println("Number of changes made: $cocolo")
##################### CLIMATE-ONLY MODEL ###########################
####################################################################
####################################################################
####################################################################
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
@load "Objects1_9/bio5_DA.jld2" bio5_DA
@load "Objects1_9/bio6_DA.jld2" bio6_DA
@load "Objects1_9/bio12_DA.jld2" bio12_DA
futurebio5_DA = bio5_DA .+ 1.0
futurebio6_DA = bio6_DA .+ 1.0
futurebio12_DA = bio12_DA .+ rand(Normal(0, 100), 125, 76)

########### bodymass vector  ##############
# Initialize an empty vector
body_mass_vector = Float64[]
# Loop through each species in spain_names and push the bodyMass into the vector
for i in spain_names
    # Find the bodyMass for the species and push it to the vector
    body_mass = gbif_sizes[gbif_sizes.species .== i, :bodyMass]
    if !isempty(body_mass)
        push!(body_mass_vector, body_mass[1])
    end
end
body_mass_vector_herps = body_mass_vector[1:49]
body_mass_vector_birds = body_mass_vector[50:256]
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

function int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(SVector{256, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

function int_Gr_for_biotic(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(SVector{256, Float64}(self_regulation * (npp .+ 0.1) .*
        1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5) .*
        1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6) .*
        1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12) .*
        state.a .* (1.0 .- (state.a ./ (npp .* binary_vector .+ npp .* opposite_binary_vector ./ 100.0)))))
end

function int_Gr_for_biotic_k(state::MyStructs256, self_regulation::AbstractFloat, k_DA::MyStructs256)
    return MyStructs256(SVector{256, Float64}(self_regulation .* (k_DA.a .+ 0.00001) .*
        state.a .* (1.0 .- (state.a ./ (k_DA.a .* binary_vector .+ k_DA.a .* opposite_binary_vector)))))
end

function lax_int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(SVector{256, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

function strict_int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs256(SVector{256, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- strict_species_niches.mean_bio5) ./ strict_species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- strict_species_niches.mean_bio6) ./ strict_species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- strict_species_niches.mean_bio12) ./ strict_species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

###################### GENERATE Ki(z) DimArray ################
###############################################################
# Initialize the empty DA arrays for each suitability method (same shape as k_DA)
DA_multiplicative = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_additive = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_geometric = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_min = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_harmonic = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

# Define the herbivore carnivore vector
herb_carv_svector = SVector{256, Float64}([name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names])
herb_carv_vector = [name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names]

# Loop through the axes of the DA arrays
for row in axes(DA_multiplicative, 1), col in axes(DA_multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate the scaled deviations for bio5, bio6, and bio12
        S_bio5 = 1 ./ (1 .+ abs.(bio5_DA[row, col] .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        S_bio6 = 1 ./ (1 .+ abs.(bio6_DA[row, col] .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)
        S_bio12 = 1 ./ (1 .+ abs.(bio12_DA[row, col] .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)

        # 1. Multiplicative Approach (Original)
        multiplicative_suitability = S_bio5 .* S_bio6 .* S_bio12
        DA_multiplicative[row, col] = MyStructs256(SVector{256, Float64}(multiplicative_suitability .* herb_carv_svector))

        # 2. Additive Approach
        additive_suitability = (S_bio5 .+ S_bio6 .+ S_bio12) / 3
        DA_additive[row, col] = MyStructs256(SVector{256, Float64}(additive_suitability .* herb_carv_svector))

        # 3. Geometric Mean Approach
        geometric_suitability = (S_bio5 .* S_bio6 .* S_bio12).^(1/3)
        DA_geometric[row, col] = MyStructs256(SVector{256, Float64}(geometric_suitability .* herb_carv_svector))

        # 4. Minimum Suitability Approach
        min_suitability = min(S_bio5, S_bio6, S_bio12)
        DA_min[row, col] = MyStructs256(SVector{256, Float64}(min_suitability .* herb_carv_svector))

        # 5. Harmonic Mean Approach
        harmonic_suitability = 3 ./ (1 ./ S_bio5 .+ 1 ./ S_bio6 .+ 1 ./ S_bio12)
        DA_harmonic[row, col] = MyStructs256(SVector{256, Float64}(harmonic_suitability .* herb_carv_svector))
    end
end

k_DA = (DA_multiplicative = DA_multiplicative, DA_additive = DA_additive, DA_geometric = DA_geometric, DA_min = DA_min, DA_harmonic = DA_harmonic)
a = maximum(k_DA.DA_multiplicative)
b = maximum(k_DA.DA_additive)
c = maximum(k_DA.DA_geometric)
d = maximum(k_DA.DA_min)
e = maximum(k_DA.DA_harmonic)
total_max = maximum([a, b, c, d, e]).b
serialize("Objects/k_DA.jls", k_DA)
k_DA = deserialize("Objects/k_DA.jls")
# k_raster = AbstractRaster[]
# for name in 1:length(k_DA)
#     prova = deepcopy(raster_with_abundances)
#     for row in axes(prova, 1), col in axes(prova, 2)
#     prova[row, col] = k_DA[name][row, col]
#     end
#     push!(k_raster, prova)
# end
# k_raster = (
#     raster_multiplicative = k_raster[1],
#     raster_additive = k_raster[2],
#     raster_geometric = k_raster[3],
#     raster_min = k_raster[4],
#     raster_harmonic = k_raster[5]
# )
# k_raster = serialize("Objects/k_raster.jls", k_raster)
k_raster = deserialize("Objects/k_raster.jls")
# map_plot(Matrix(k_DA.DA_multiplicative); type = "heatmap", palette = :thermal, title = "Geometric Suitability", colorrange = (0, total_max))
# map_plot(Matrix(k_DA.DA_additive); type = "heatmap", palette = :thermal, title = "Additive Suitability", colorrange = (0, total_max))
# map_plot(Matrix(k_DA.DA_geometric); type = "heatmap", palette = :thermal, title = "Geometric Suitability", colorrange = (0, total_max))
# map_plot(Matrix(k_DA.DA_min); type = "heatmap", palette = :thermal, title = "Minimum Suitability", colorrange = (0, total_max))
# map_plot(Matrix(k_DA.DA_harmonic); type = "heatmap", palette = :thermal, title = "Harmonic Suitability", colorrange = (0, total_max))

###################### GENERATE lambda_DA DimArray ################
###############################################################
# Function to calculate the lambda scalar for a given k_hat and NPP
function calculate_lambda_scalar(k_hat, NPP)
    lambda = NPP / sum(k_hat)
    return lambda
end

# Initialize a NamedTuple to store lambda_DA for each k_DA
lambda_DA = NamedTuple((
    multiplicative = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    additive = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    geometric = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    minimum = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    harmonic = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
))
# lambda_DA[1]

# Loop through the axes to calculate lambda_DA for each suitability method
for row in axes(lambda_DA.multiplicative, 1), col in axes(lambda_DA.multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate lambda for multiplicative suitability
        lambda_DA.multiplicative[row, col] = calculate_lambda_scalar(k_DA.DA_multiplicative[row, col].a, npp_DA[row, col])
        # Calculate lambda for additive suitability
        lambda_DA.additive[row, col] = calculate_lambda_scalar(k_DA.DA_additive[row, col].a, npp_DA[row, col])
        # Calculate lambda for geometric suitability
        lambda_DA.geometric[row, col] = calculate_lambda_scalar(k_DA.DA_geometric[row, col].a, npp_DA[row, col])
        # Calculate lambda for minimum suitability
        lambda_DA.minimum[row, col] = calculate_lambda_scalar(k_DA.DA_min[row, col].a, npp_DA[row, col])
        # Calculate lambda for harmonic suitability
        lambda_DA.harmonic[row, col] = calculate_lambda_scalar(k_DA.DA_harmonic[row, col].a, npp_DA[row, col])
    end
end

# lambda_raster = AbstractRaster[]
# for name in keys(lambda_DA)
#     prova = deepcopy(raster_richness)
#     for row in axes(lambda_DA[name], 1), col in axes(lambda_DA[name], 2)
#         prova[row, col] = lambda_DA[name][row, col]
#     end
#     push!(lambda_raster, prova)
# end
# lambda_raster = (
#     multiplicative = lambda_raster[1],
#     additive = lambda_raster[2],
#     geometric = lambda_raster[3],
#     min = lambda_raster[4],
#     harmonic = lambda_raster[5]
# )
# lambda_raster = serialize("Objects/lambda_raster.jls", lambda_raster)
lambda_raster = deserialize("Objects/lambda_raster.jls")
##########################################
##########################################
prop = [28/256, 57/256, 103/256, 68/256]
# Load and prepare the data
belonging = CSV.File("DFs\\block_per_species.csv") |> DataFrame
order_indices = indexin(names(iberian_interact_df), belonging[:, "Column1"])
belonging = belonging[order_indices, :]
species_block = belonging[:, 2]  # This array indicates the block each species belongs to

function trophic_optimized(abundances, full_IM)
    # Calculate the weighted interaction directly
    return MyStructs256(SVector{256, Float64}((full_IM * abundances.a) .* abundances.a))
end

dimensions = (125, 76)
idx_tupled = [(i, j) for i in 1:dimensions[1], j in 1:dimensions[2]]

climatic_niche_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    return state + int_Gr(state, self_regulation, npp, bio5, bio6, bio12) 
end

lax_climatic_niche_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    return state + lax_int_Gr(state, self_regulation, npp, bio5, bio6, bio12)
end

strict_climatic_niche_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    return state + strict_int_Gr(state, self_regulation, npp, bio5, bio6, bio12)
end

biotic_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    merged_state = state + 
        int_Gr_for_biotic(state, self_regulation, npp, bio5, bio6, bio12)  +
        trophic_optimized(state, full_IM)
    return MyStructs256(max.(0.0, merged_state.a))
end

biotic_rule_k = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return MyStructs256(
        SVector{256, Float64}(
            max.(
                0.0,
                (state +
                int_Gr_for_biotic_k(state, self_regulation, k_DA)  +
                trophic_optimized(state, full_IM)).a
            )
        )
    )
end

biotic_rule_k_herps = Cell{Tuple{:herps, :birmmals, :k_DA}, :herps}() do data, (herps, birmmals, k_DA), I
    # if any(isinf, birmmals.a) || any(isnan, birmmals.a)
    #     @error "state has NA values in birmmals"
    #     println(I)
    # end
    # if any(isinf, herps.a) || any(isnan, herps.a)
    #     @error "state has NA values in herps"
    #     println(I)
    # end
    # if I == study_cell
    #     push!(abundances, (deepcopy(herps) + deepcopy(birmmals)).a)
    # end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyHerps(SVector{49, Float64}(max.(0.0000000001, merged_state.a[1:49])))
end

biotic_rule_k_birmmals = Cell{Tuple{:herps, :birmmals, :k_DA}, :birmmals}() do data, (herps, birmmals, k_DA), I
    # if typeof(birmmals) != MyBirmmals{Float64}
    #     @error "birmmals is not MyBirmmals"
    # end
    # if typeof(herps) != MyHerps{Float64}
    #     @error "herps is not MyHerps"
    # end
    # if any(isinf, birmmals.a) || any(isnan, birmmals.a)
    #     @error "state has NA values in birmmals"
    #     println(I)
    # end
    # if any(isinf, herps.a) || any(isnan, herps.a)
    #     @error "state has NA values in herps"
    #     println(I)
    # end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyBirmmals(SVector{207, Float64}(max.(0.0000000001, merged_state.a[50:256])))
end

resource_rule = Cell{Tuple{:state, :npp, :bio5, :bio6, :bio12}, :state}() do data, (state, npp, bio5, bio6, bio12), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    merged_state = state + resource_int_Gr(state, self_regulation, npp, bio5, bio6, bio12)
    return MyStructs256(max.(0.0, merged_state.a))
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

alfa = 0.1
indisp = InwardsDispersal{:state, :state}(;
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30)
);

outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);
outdisp_birmmals = OutwardsDispersal{:birmmals, :birmmals}(;
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);
outdisp_herps = OutwardsDispersal{:herps, :herps}(;
    formulation=CustomKernel(alfa/2),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);
# full_IM_used = deepcopy(full_IM)
# full_IM_used = full_IM_used*0.1
highest_npp = maximum(filter(!isnan, npp_DA[:]))
highest_richness = maximum(filter(!iszero, DA_richness[:]))

# function random_dimarray(dimensions::Tuple{Int64, Int64}; prevalence = nothing)
#     init_array = DimArray(zeros(dimensions, MyStructs256{Float32}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
#     for i in idx
#         init_array[i] = MyStructs256((0.1 .* binary_vector .+ 0.1 .* opposite_binary_vector) .* sample([1,0], Weights([prevalence, 1-prevalence]), 256))
#     end
#     return init_array
# end

# DA_random_with_abundances = random_dimarray(dimensions; prevalence = 0.1)

# prevalence = 0.1
# function random_herp_dimarray(dimensions::Tuple{Int64, Int64}; prevalence = nothing)
#     init_array = DimArray(zeros(dimensions, MyHerps{Float32}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
#     for i in idx
#         init_array[i] = MyHerps((initial_abundance  .* fill(1.0, 49)) .* sample([1,0], Weights([prevalence, 1-prevalence]), 49))
#     end
#     return init_array
# end
# DA_random_herps_with_abundances = random_herp_dimarray(dimensions; prevalence = prevalence)

# function random_birmmals_dimarray(dimensions::Tuple{Int64, Int64}, prev)
#     init_array = DimArray(zeros(dimensions, MyBirmmals{Float32}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
#     for i in idx
#         init_array[i] = MyBirmmals((initial_abundance .* fill(1.0, 207)) .* sample([1,0], Weights([prev, 1-prev]), 207))
#     end
#     return init_array
# end
# DA_random_birmmals_with_abundances = random_birmmals_dimarray(dimensions, prevalence)

##### METRICS FOR THREADING QUICKLY #####
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

# # Iterate over the grid and create new SVector replacing 10.0 with 500.0
# triall = deserialize("Objects\\DA_with_abundances_all10.jls")::DimArray{MyStructs256{Float64},2}
# for row in axes(triall, 1), col in axes(triall, 2)
#     old_vector = triall[row, col].a
#     new_vector = SVector{256, Float64}(replace(old_vector, 10.0 => 500.0))  # Create a new SVector
#     triall[row, col] = MyStructs256(new_vector)  # Assign the new MyStructs256 with updated vector
# end
# richness_similarity(Matrix(triall), modified = true, caca = true)

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
# alive_predators(triall, modified = true, caca = true)

function total_biomass(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)).*lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state.*lambda_DA[position]
    end
    if caca
        combined_abundances = array_output.*lambda_DA[position]
    end
    return sum(combined_abundances).b
end