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
web = filter(row -> row.targetLifestageName == "adults", web) ## TODO THIS IS A KEY CHANGEEEEE

web = DataFrame(predator = web.sourceTaxonName, prey = web.targetTaxonName)

check_species_type = filter(row -> row.predator == "Pelobates cultripes", web) # Just add the species name to check if it's a predator

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

#### STABLISHING TROPHIC LEVELS AND HERBIVORE NAMES
TrophInd = CSV.File("DFs/TLs.csv") |> DataFrame
TrophInd = TrophInd[1:256, 2:3]
TrophInd[:, 1] = round.(TrophInd[:, 1], digits = 2)
# TrophInd[:, 2] = TrophInd[:, 2].-1
# TrophInd[256, 2] = 1.0 # For some reason last line had floating point error
rename!(TrophInd, Symbol("species") => :Species, Symbol("TrophInd") => :TL)
TrophInd[:, :TL] = round.(TrophInd[:, :TL].-1, digits = 2)
order_indices = indexin(spain_names, TrophInd[:, :Species])
TrophInd = TrophInd[order_indices, :]
TrophInd_vector = TrophInd[:, :TL]

# herbivore_names = TrophInd[TrophInd[:, :TL] .== 1, :Species]

# Turn iberian_interact_matrix into a DataFrame
# Convert the iberian_interact_matrix into a DataFrame with appropriate column and row names
iberian_interact_df = DataFrame(iberian_interact_matrix, spain_names)
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
iberian_interact_NA = NamedArray(Matrix(iberian_interact_df), (names(iberian_interact_df), names(iberian_interact_df)))
fill_diagonal!(iberian_interact_NA, 0.0)
# serialize("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\Networks\\DFs\\iberian_interact_NA.jls", iberian_interact_NA)

herbivore_names = []
for i in axes(iberian_interact_NA, 1)
    if all(x -> x == 0, iberian_interact_NA[i, :])
        push!(herbivore_names, names(iberian_interact_NA, 1)[i])
    end
end
predator_names = setdiff(spain_names, herbivore_names)

println("The number of 1s in the column Rupicabra pyrenaica is $(sum(iberian_interact_NA[:, "Rupicapra pyrenaica"] .== 1))")

non_zero_cols = names(iberian_interact_NA, 2)[findall(x -> x != 0, iberian_interact_NA["Coracias garrulus", :])]
println("The non-zero columns for Coracias garrulus are $non_zero_cols")

########### VISUAL OUTPUTS ################
###########################################
# Visualization settings
DynamicGrids.to_rgb(scheme::ObjectScheme, obj::MyStructs254) = ARGB32(
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0)
)
DynamicGrids.to_rgb(scheme, obj::MyStructs254) = get(scheme, clamp(obj.b, 0.0, 1.0))
##################### MWE_TvdD ##################################
####################################################################
####################################################################
####################################################################
# For a heatmap we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs254, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyBirmmals, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyHerps, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs254, 2})
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
# WITH LAMBDA
# For MyStructs
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs254, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] .* lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# PLOT
# For MyStructs
# function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs254, 2})
#     richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_DA.multiplicative)
#     return Makie.convert_arguments(t, richness)
# end
# MK.image(Matrix(des_file_to_try); colomap = custom_palette)
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
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyStructs254, 2})
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
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyStructs254, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
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
# species_df.sum = [sum(species_df[i, 5:260]) for i in 1:size(species_df, 1)] # OLD
species_df.sum = [sum(species_df[i, 5:258]) for i in 1:size(species_df, 1)] # NEW
species_df_matrix = Matrix(species_df)

utmraster = Raster("Rasters\\updated_utmraster.tif")
utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)

######################### NPP ####################################
npp_absolute = CSV.File("DFs\\npp_absolute_df.csv") |> DataFrame
rename!(npp_absolute, [:ID, :UTMCODE, :npp, :X, :Y]) 
npp_absolute_in_kg = deepcopy(npp_absolute)
npp_absolute_in_kg.npp = npp_absolute.npp .* 1000
npp_absolute_in_kg = npp_absolute_in_kg[:, [2, 3]]

species_df = leftjoin(species_df, npp_absolute_in_kg, on = :UTMCODE, makeunique = true)
species_df_matrix = Matrix(species_df)

npp_DA = deserialize("Objects\\npp_DA.jls")

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

##################### NEW NICHES ###########################
######## bio rasters  ##############
bio5_DA = deserialize("Objects\\bio5.jls")
bio6_DA = deserialize("Objects\\bio6.jls")
bio12_DA = deserialize("Objects\\bio12.jls")

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
body_mass_vector_birds = body_mass_vector[50:254]
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

# herbivore_names = CSV.File(joinpath(meta_path, "herbivore_names.csv")) |> DataFrame
# herbivore_names = convert(Vector{String}, herbivore_names[:, 2])
# binary_vector = [name in herbivore_names ? 1 : 0 for name in names(iberian_interact_df)]
# opposite_binary_vector = [name in herbivore_names ? 0 : 1 for name in names(iberian_interact_df)]

function int_Gr(state::MyStructs254, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs254(SVector{254, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

function int_Gr_for_biotic(state::MyStructs254, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs254(SVector{254, Float64}(self_regulation * (npp .+ 0.1) .*
        1 ./ (1 .+ abs.(bio5 .- species_niches.mean_bio5) ./ species_niches.sd_bio5) .*
        1 ./ (1 .+ abs.(bio6 .- species_niches.mean_bio6) ./ species_niches.sd_bio6) .*
        1 ./ (1 .+ abs.(bio12 .- species_niches.mean_bio12) ./ species_niches.sd_bio12) .*
        state.a .* (1.0 .- (state.a ./ (npp .* binary_vector .+ npp .* opposite_binary_vector ./ 100.0)))))
end

function int_Gr_for_biotic_k(state::MyStructs254, self_regulation::AbstractFloat, k_DA::MyStructs254)
    return MyStructs254(SVector{254, Float64}(self_regulation .* (k_DA.a .+ 0.00001) .*
        state.a .* (1.0 .- (state.a ./ (k_DA.a .* binary_vector .+ k_DA.a .* opposite_binary_vector)))))
end

function lax_int_Gr(state::MyStructs254, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs254(SVector{254, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

function strict_int_Gr(state::MyStructs254, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, bio5::AbstractFloat, bio6::AbstractFloat, bio12::AbstractFloat)
    return MyStructs254(SVector{254, Float64}(self_regulation * (npp + 0.1) .*
        (1 ./ (1 .+ abs.(bio5 .- strict_species_niches.mean_bio5) ./ strict_species_niches.sd_bio5)) .*
        (1 ./ (1 .+ abs.(bio6 .- strict_species_niches.mean_bio6) ./ strict_species_niches.sd_bio6)) .*
        (1 ./ (1 .+ abs.(bio12 .- strict_species_niches.mean_bio12) ./ strict_species_niches.sd_bio12)) .*
        state.a .* (1.0 - (state.b / (npp + 0.1)))))
end

###################### GENERATE Ki(z) DimArray ################
###############################################################
# Initialize the empty DA arrays for each suitability method (same shape as k_DA)
DA_multiplicative = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_additive = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_geometric = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_min = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_harmonic = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

# Define the herbivore carnivore vector
herb_carv_svector = SVector{254, Float64}([name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names])
# herb_carv_vector = [name in herbivore_names ? 1.0 : 0.00000001 for name in spain_names]

DA_sum = deserialize("Objects\\DA_sum.jls")
idx = findall(x -> x == 1.0, DA_sum)
# Loop through the axes of the DA arrays
for row in axes(DA_multiplicative, 1), col in axes(DA_multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate the scaled deviations for bio5, bio6, and bio12
        S_bio5 = 1 ./ (1 .+ abs.(bio5_DA[row, col] .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        S_bio6 = 1 ./ (1 .+ abs.(bio6_DA[row, col] .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)
        S_bio12 = 1 ./ (1 .+ abs.(bio12_DA[row, col] .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)

        # 1. Multiplicative Approach (Original)
        multiplicative_suitability = S_bio5 .* S_bio6 .* S_bio12
        DA_multiplicative[row, col] = MyStructs254(SVector{254, Float64}(multiplicative_suitability .* herb_carv_svector))

        # 2. Additive Approach
        additive_suitability = (S_bio5 .+ S_bio6 .+ S_bio12) / 3
        DA_additive[row, col] = MyStructs254(SVector{254, Float64}(additive_suitability .* herb_carv_svector))

        # 3. Geometric Mean Approach
        geometric_suitability = (S_bio5 .* S_bio6 .* S_bio12).^(1/3)
        DA_geometric[row, col] = MyStructs254(SVector{254, Float64}(geometric_suitability .* herb_carv_svector))

        # 4. Minimum Suitability Approach
        min_suitability = min(S_bio5, S_bio6, S_bio12)
        DA_min[row, col] = MyStructs254(SVector{254, Float64}(min_suitability .* herb_carv_svector))

        # 5. Harmonic Mean Approach
        harmonic_suitability = 3 ./ (1 ./ S_bio5 .+ 1 ./ S_bio6 .+ 1 ./ S_bio12)
        DA_harmonic[row, col] = MyStructs254(SVector{254, Float64}(harmonic_suitability .* herb_carv_svector))
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

# DA = DimArray(reshape([MyStructs254(SVector{254, Float64}(fill(0.0, 254))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
# # Iterate over species_df and populate DA
# for i in 1:size(species_df, 1)
#     for j in 1:125*76
#         if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
#             if species_df.sum[i] == 54 
#             #    println(utmraster_DA[3]) 
#             end
#             # Convert species_df_matrix[i, 5:260] to SVector{256, Float64} before creating MyStructs256
#             DA[j] = MyStructs254(SVector{254, Float64}(species_df_matrix[i, 5:258]))
#         end
#     end
# end
# serialize("Objects\\DA_254.jls", DA)
DA = deserialize("Objects\\DA_254.jls")

DA_birmmals = DimArray(reshape([MyBirmmals(SVector{205, Float64}(fill(0.0, 205))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in idx
    DA_birmmals[i] = MyBirmmals(SVector{205, Float64}(DA[i].a[50:end]))
end
