# Define the size and bounds for the raster
size_x, size_y = 100, 100
x_range = X(LinRange(-180.0, 180.0, size_x))  # 100 evenly spaced points between -180 and 180
y_range = Y(LinRange(-90.0, 90.0, size_y))    # 100 evenly spaced points between -90 and 90

# Create a function to generate MyStructs256 objects
function create_mystructs256()
    data = SVector{256, Float64}(rand(Float64, 256))
    return MyStructs256(data)
end

# Generate raster data
raster_data = [create_mystructs256() for _ in 1:(size_x * size_y)]

# Create a Raster with MyStructs256
raster = Raster(
    reshape(raster_data, size_x, size_y),
    dims = (x_range, y_range)
)

# Load the original raster and convert it to a boolean mask
utmraster = Raster("Rasters\\updated_utmraster.tif")
utmraster_data = map(x -> isnothing(x) || isnan(x) ? false : true, Matrix(utmraster))

# Extract the dimensions from the original raster
x_dim = dims(utmraster, X)
y_dim = dims(utmraster, Y)

# Create a Raster for the utmraster_da using the extracted dimensions
utmraster_r = Raster(
    utmraster_data,
    dims = (x_dim, y_dim)
)

species_df = CSV.File("DFs\\Species_spain_df.csv") |> DataFrame

variables = species_df[!, 2:5]
rename!(variables, [:ID, :Value, :sum, :UTMCODE])

species = species_df[!, unique_species_in_web]
species = species[:, spain_names]

species_df = hcat(variables, species, makeunique=true)
species_df_matrix = Matrix(species_df)

# Extract species and variables data
variables = species_df[!, 2:5]
species = species_df[!, 6:end]

# Convert species DataFrame to Matrix
species_df_matrix = Matrix(species_df)

# Create a Raster for DA with MyStructs256
raster_data_DA = [MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76]
raster_DA = Raster(reshape(raster_data_DA, 125, 76), dims=dims(utmraster))

# Populate the raster_DA
for i in 1:size(species_df, 1)
    for j in 1:(125*76)
        if Float32(species_df.Value[i]) == Float32(utmraster[j])
            raster_DA[j] = MyStructs256(SVector{256, Float64}(species_df_matrix[i, 5:260]))
        end
    end
end

# Save the raster
# serialize("Objects\\raster_DA.jls", raster_DA)
raster_DA = deserialize("Objects\\raster_DA.jls")

# Initialize raster for DA_herps
raster_data_herps = [MyHerps(SVector{49, Float64}(fill(0.0, 49))) for _ in 1:125*76]
raster_herps = Raster(reshape(raster_data_herps, 125, 76), dims=dims(utmraster))

# Populate raster_herps
for i in 1:size(species_df, 1)
    for j in 1:(125*76)
        if Float32(species_df.Value[i]) == Float32(utmraster[j])
            raster_herps[j] = MyHerps(SVector{49, Float64}(raster_DA[j].a[1:49]))
        end
    end
end

# Save the raster
# serialize("Objects\\raster_herps.jls", raster_herps)
raster_herps = deserialize("Objects\\raster_herps.jls")

# Initialize raster for DA_birmmals
raster_data_birmmals = [MyBirmmals(SVector{207, Float64}(fill(0.0, 207))) for _ in 1:125*76]
raster_birmmals = Raster(reshape(raster_data_birmmals, 125, 76), dims=dims(utmraster))

# Populate raster_birmmals
for i in 1:size(species_df, 1)
    for j in 1:(125*76)
        if Float32(species_df.Value[i]) == Float32(utmraster[j])
            raster_birmmals[j] = MyBirmmals(SVector{207, Float64}(raster_DA[j].a[50:256]))
        end
    end
end

# Save the raster
# serialize("Objects\\raster_birmmals.jls", raster_birmmals)
raster_birmmals = deserialize("Objects\\raster_birmmals.jls")

# Adjust abundances for raster_DA
raster_with_abundances = deepcopy(raster_DA)
for row in axes(raster_DA, 1), col in axes(raster_DA, 2)
    if raster_DA[row, col] != MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        new_a = SVector{256, Float64}([raster_DA[row, col].a[i] != 0.0 ? initial_abundance : raster_DA[row, col].a[i] for i in 1:256])
        raster_with_abundances[row, col] = MyStructs256(new_a)
    end
end
# serialize("Objects\\raster_with_abundances.jls", raster_with_abundances)
raster_with_abundances = deserialize("Objects\\raster_with_abundances.jls")

# Adjust abundances for raster_birmmals
raster_birmmals_with_abundances = deepcopy(raster_birmmals)
for row in axes(raster_birmmals, 1), col in axes(raster_birmmals, 2)
    if raster_birmmals[row, col] != MyBirmmals(SVector{207, Float64}(fill(0.0, 207)))
        new_a = SVector{207, Float64}([raster_birmmals[row, col].a[i] != 0.0 ? initial_abundance : raster_birmmals[row, col].a[i] for i in 1:207])
        raster_birmmals_with_abundances[row, col] = MyBirmmals(new_a)
    end
end
# serialize("Objects\\raster_birmmals_with_abundances.jls", raster_birmmals_with_abundances)
raster_birmmals_with_abundances = deserialize("Objects\\raster_birmmals_with_abundances.jls")

# Adjust abundances for raster_herps
raster_herps_with_abundances = deepcopy(raster_herps)
for row in axes(raster_herps, 1), col in axes(raster_herps, 2)
    if raster_herps[row, col] != MyHerps(SVector{49, Float64}(fill(0.0, 49)))
        new_a = SVector{49, Float64}([raster_herps[row, col].a[i] != 0.0 ? initial_abundance : raster_herps[row, col].a[i] for i in 1:49])
        raster_herps_with_abundances[row, col] = MyHerps(new_a)
    end
end
# serialize("Objects\\raster_herps_with_abundances.jls", raster_herps_with_abundances)
raster_herps_with_abundances = deserialize("Objects\\raster_herps_with_abundances.jls")

raster_richness = deepcopy(utmraster)
for i in 1:(125*76)
    raster_richness[i] = DA_richness[i]
end
# serialize("Objects\\raster_richness.jls", raster_richness)
raster_richness = deserialize("Objects\\raster_richness.jls")

raster_richness_birmmals = deepcopy(utmraster)
for i in 1:(125*76)
    raster_richness_birmmals[i] = DA_richness_birmmals[i]
end
# serialize("Objects\\raster_richness_birmmals.jls", raster_richness_birmmals)
raster_richness_birmmals = deserialize("Objects\\raster_richness_birmmals.jls")

raster_richness_herps = deepcopy(utmraster)
for i in 1:(125*76)
    raster_richness_herps[i] = DA_richness_herps[i]
end
# serialize("Objects\\raster_richness_herps.jls", raster_richness_herps)
raster_richness_herps = deserialize("Objects\\raster_richness_herps.jls")

# Extract the dimensions from the original raster
x_dim = dims(utmraster, X)
y_dim = dims(utmraster, Y)
# Create a new boolean raster with the same dimensions as utmraster
raster_sum_data = reshape(Bool.(DA_sum), (length(x_dim), length(y_dim)))
# Create the new Raster object with boolean data
raster_sum = Raster(raster_sum_data, dims=(x_dim, y_dim))
# serialize("Objects\\raster_sum.jls", raster_sum)
raster_sum = deserialize("Objects\\raster_sum.jls")