gbif_sizes = CSV.read("DFs/gbif_sizes.csv", DataFrame)[:, 2:end]

DA = deserialize("Objects/DA.jls")

DA_birmmals = deserialize("Objects/DA_birmmals.jls")

DA_sum = deserialize("Objects/DA_sum.jls")

DA_richness = deserialize("Objects/DA_richness.jls")::DimArray{Float64,2}

########################## IDX #####################################
idx = findall(x -> x == 1.0, DA_sum)

npp_DA = deserialize("Objects/npp_DA.jls")

################### iberian_interact_NA #####################
####################################################################
####################################################################
####################################################################

if true # I don't know why I did this but I think it messes things up
# Load a DataFrame from a serialized file ('.jls' format).
iberian_interact_df = deserialize("Objects/iberian_interact_df.jls")

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
end

##################### CLIMATE-ONLY MODEL ###########################
####################################################################
####################################################################
####################################################################
##################### NEW NICHES ###########################
######## bio rasters  ##############
bio5_DA = deserialize("Objects/bio5.jls")
bio6_DA = deserialize("Objects/bio6.jls")
bio12_DA = deserialize("Objects/bio12.jls")

######## niches_df  ##############
lax_species_niches = CSV.File("DFs/iberian_species_niches_withLaxNiche.csv", decimal = ',') |> DataFrame
order_indices = indexin(spain_names, lax_species_niches[:, :Species])
lax_species_niches = lax_species_niches[order_indices, :]

