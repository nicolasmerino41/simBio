average_densities = CSV.File(joinpath("C:/Users", PC, "OneDrive/PhD/GitHub/Networks/DFs/average_densities.csv")) |> DataFrame
santini_names = average_densities[:, 2]

# Assuming `spain_names` and `santini_names` are vectors
matches_indices = findall(x -> x in santini_names, spain_names)

#######################################################################

# Assumptions: spain_names, gbif_sizes, and average_densities are already loaded
# spain_names: Vector of species names
# gbif_sizes: DataFrame with columns "species" and "bodymass"
# average_densities: DataFrame with columns "names" and "mean"

# Step 1: Filter spain_names to only include positions 50 to 256
mammals_birds_names = spain_names[50:254]

# Step 2: Prepare lookup dictionaries for bodymass and mean_density
bodymass_dict = Dict(row.species => row.bodyMass for row in eachrow(gbif_sizes))
mean_density_dict = Dict(row.names => row.mean for row in eachrow(average_densities))

# Step 3: Construct the final DataFrame
species = mammals_birds_names
mean_density = [get(mean_density_dict, name, missing) for name in mammals_birds_names]
bodyMass = [get(bodymass_dict, name, missing) for name in mammals_birds_names]

# Calculate biomass (product of mean_density and bodyMass), handle missing values
biomass = [(ismissing(md) || ismissing(bm)) ? missing : md * bm for (md, bm) in zip(mean_density, bodyMass)]

# Step 4: Create the DataFrame
birmmals_biomass = DataFrame(
    species = species,
    mean_density = mean_density,
    bodyMass = bodyMass,
    biomass = biomass
)

# Display the final DataFrame
println("Final DataFrame 'birmmals_biomass':")
println(birmmals_biomass)

##############################################################################
# Fixing missing species
birmmals_biomass_fixed = deepcopy(birmmals_biomass)
for i in axes(birmmals_biomass, 1)
    if ismissing(birmmals_biomass[i, :biomass])
        birmmals_biomass_fixed[i, :mean_density] = median(filter(!ismissing, birmmals_biomass[:, :mean_density]))
        birmmals_biomass_fixed[i, :biomass] = birmmals_biomass_fixed[i, :mean_density] * birmmals_biomass[i, :bodyMass]
    end
end

println(birmmals_biomass_fixed)

DA_birmmals_with_pi = deepcopy(DA_birmmals)
for i in idx
    vect = DA_birmmals[i].a
    vect1 = Vector(DA_birmmals[i].a)
    for j in 1:length(vect)
        if isone(vect[j])
            vect1[j] = birmmals_biomass_fixed[j, :biomass]
        end
    end
    DA_birmmals_with_pi[i] = MyBirmmals(SVector{205, Float64}(vect1))
end

# map_plot(DA_birmmals_with_pi; type = "heatmap", palette = :thermal)

DA_birmmals_with_pi[idx[20]].a

#### INVESTIGATING BIG BIOMASS DIFFERENCES ####
birmmals_biomass_only_herbivores = 
    birmmals_biomass_fixed[map(s -> s in herbivore_names, birmmals_biomass_fixed[:, :species]), :]

# valueee = 0
# for i in birmmals_biomass_only_herbivores[:, 1]
#     if i in herbivore_names
#         valueee += 1
#     end
#     println(valueee)
# end
