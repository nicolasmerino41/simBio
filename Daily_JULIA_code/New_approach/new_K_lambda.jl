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
        DA_multiplicative[row, col] = MyStructs256(SVector{256, Float64}(multiplicative_suitability))

        # 2. Additive Approach
        additive_suitability = (S_bio5 .+ S_bio6 .+ S_bio12) / 3
        DA_additive[row, col] = MyStructs256(SVector{256, Float64}(additive_suitability))

        # 3. Geometric Mean Approach
        geometric_suitability = (S_bio5 .* S_bio6 .* S_bio12).^(1/3)
        DA_geometric[row, col] = MyStructs256(SVector{256, Float64}(geometric_suitability))

        # 4. Minimum Suitability Approach
        min_suitability = min(S_bio5, S_bio6, S_bio12)
        DA_min[row, col] = MyStructs256(SVector{256, Float64}(min_suitability))

        # 5. Harmonic Mean Approach
        harmonic_suitability = 3 ./ (1 ./ S_bio5 .+ 1 ./ S_bio6 .+ 1 ./ S_bio12)
        DA_harmonic[row, col] = MyStructs256(SVector{256, Float64}(harmonic_suitability))
    end
end

new_k_DA = (DA_multiplicative = DA_multiplicative, DA_additive = DA_additive, DA_geometric = DA_geometric, DA_min = DA_min, DA_harmonic = DA_harmonic)
a = maximum(new_k_DA.DA_multiplicative)
b = maximum(new_k_DA.DA_additive)
c = maximum(new_k_DA.DA_geometric)
d = maximum(new_k_DA.DA_min)
e = maximum(new_k_DA.DA_harmonic)
total_max = maximum([a, b, c, d, e]).b
serialize("Objects/new_k_DA.jls", new_k_DA)
new_k_DA = deserialize("Objects/new_k_DA.jls")
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
new_lambda_DA = NamedTuple((
    multiplicative = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    additive = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    geometric = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    minimum = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76))),
    harmonic = DimArray(reshape([NaN64 for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
))

# Loop through the axes to calculate new_lambda_DA for each suitability method
for row in axes(new_lambda_DA.multiplicative, 1), col in axes(new_lambda_DA.multiplicative, 2)
    if isone(DA_sum[row, col])
        # Calculate lambda for multiplicative suitability
        new_lambda_DA.multiplicative[row, col] = calculate_lambda_scalar(k_DA.DA_multiplicative[row, col].a, npp_DA[row, col])
        # Calculate lambda for additive suitability
        new_lambda_DA.additive[row, col] = calculate_lambda_scalar(k_DA.DA_additive[row, col].a, npp_DA[row, col])
        # Calculate lambda for geometric suitability
        new_lambda_DA.geometric[row, col] = calculate_lambda_scalar(k_DA.DA_geometric[row, col].a, npp_DA[row, col])
        # Calculate lambda for minimum suitability
        new_lambda_DA.minimum[row, col] = calculate_lambda_scalar(k_DA.DA_min[row, col].a, npp_DA[row, col])
        # Calculate lambda for harmonic suitability
        new_lambda_DA.harmonic[row, col] = calculate_lambda_scalar(k_DA.DA_harmonic[row, col].a, npp_DA[row, col])
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
