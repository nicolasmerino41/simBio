gbif_sizes = CSV.read("DFs/gbif_sizes.csv", DataFrame)[:, 2:end]

########################## IDX #####################################
DA_sum = deserialize("Objects/DA_sum.jls")
idx = findall(x -> x == 1.0, DA_sum)

npp_DA = deserialize("Objects/npp_DA.jls")

utmraster = Raster("Rasters/updated_utmraster.tif")
utmraster_DA = DimArray(utmraster)
# CSV.write("DFs/species_df_254.csv", species_df)
species_df = CSV.read("DFs/species_df_254.csv", DataFrame)
species_df_matrix = Matrix(species_df)

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
# serialize("Objects/DA_254.jls", DA)
DA = deserialize("Objects/DA_254.jls")

DA_birmmals = DimArray(reshape([MyBirmmals(SVector{205, Float64}(fill(0.0, 205))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in idx
    DA_birmmals[i] = MyBirmmals(SVector{205, Float64}(DA[i].a[50:end]))
end

DA_richness = deserialize("Objects/DA_richness.jls")::DimArray{Float64,2}

################### iberian_interact_NA #####################
####################################################################
####################################################################
####################################################################
# CSV.write("DFs/spain_names.csv", DataFrame(spain_names = spain_names))
spain_names = CSV.read("DFs/spain_names.csv", DataFrame)[:, :spain_names]
spain_names = convert(Vector{String}, spain_names)

# dedo = Matrix(iberian_interact_NA)
# dedo_df = DataFrame(dedo, :auto)
# CSV.write("DFs/dedo_matrix.csv", dedo_df)

dedoo = CSV.File("DFs/dedo_matrix.csv", decimal = '.') |> DataFrame
iberian_interact_NA = NamedArray(Matrix(dedoo), (spain_names, spain_names), ("A", "B")) 
# iberian_interact_NA = deserialize("Objects/iberian_interact_NA.jls") #TODO this does not work due to julia version

herbivore_names = []
for i in axes(iberian_interact_NA, 1)
    if all(x -> x == 0, iberian_interact_NA[i, :])
        push!(herbivore_names, names(iberian_interact_NA, 1)[i])
    end
end
predator_names = setdiff(spain_names, herbivore_names)
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

