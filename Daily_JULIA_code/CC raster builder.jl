bounds = X(Rasters.Between(-10, 4)), Y(Rasters.Between(35.5, 44))
futurebio5_raster = RS.resample(Raster("C:/Users/MM-1/Downloads/CHELSA_bio5_2071-2100_gfdl-esm4_ssp585_V.2.1.tif", lazy = true)[bounds...], to = utmraster)
futurebio5_raster = Rasters.map(x -> Float32(x), futurebio5_raster)

bio5_raster = Raster("Rasters/bio5.tif")

MK.heatmap(futurebio5_raster, colormap = :thermal)
MK.heatmap(bio5_raster, colormap = :thermal)

maximum(filter(!isnan, bio5_raster))
minimum(filter(!isnan, bio5_raster))
for i in axes(DA_sum, 1), j in axes(DA_sum, 2)
    if iszero(DA_sum[i, j])
        futurebio5_raster[i, j] = NaN
    else
        futurebio5_raster[i, j] = futurebio5_raster[i, j] / 10 - 273.15
    end
end
t = maximum(filter(!isnan, futurebio5_raster))
o = minimum(filter(!isnan, futurebio5_raster))

f = Figure(resolution = (600, 400))
ax1 = Axis(f[1, 1])
ax2 = Axis(f[1, 2])
MK.heatmap!(ax1, bio5_raster; colormap=custom_palette, colorrange=(18.52, 36.41))
MK.heatmap!(ax2, futurebio5_raster; colormap=custom_palette, colorrange=(o, t))
f

function create_raster_series(origin_raster::AbstractRaster, dest_raster::AbstractRaster, num_steps::Int, start_date::Date, interval::Period)
    # Ensure the dimensions of origin and destination rasters are the same
    if size(origin_raster) != size(dest_raster)
        error("Origin and destination rasters must have the same dimensions.")
    end

    # Preallocate the raster list with a specific concrete type
    raster_list = Raster[]  # Use a more concrete type than AbstractRaster

    # Generate dates for the series based on the start date, number of steps, and interval type
    dates = [start_date + i * interval for i in 0:num_steps-1]

    # Generate interpolated rasters
    for step in 0:num_steps-1
        fraction = step / (num_steps - 1)
        # Linear interpolation between origin and destination rasters
        new_raster = origin_raster * (1 - fraction) + dest_raster * fraction
        push!(raster_list, new_raster)
    end

    # Create a RasterSeries with given time steps
    raster_series = RasterSeries(raster_list, Ti(dates))

    return raster_series
end

# Example usage:
# Assuming origin_raster and dest_raster are previously loaded or created rasters
num_steps = 10
start_date = Date(2024, 1, 1)
interval = Month(1)

# Call the function with origin_raster and dest_raster (bio5_raster and futurebio5_raster)
raster_series = create_raster_series(bio5_raster, futurebio5_raster, num_steps, start_date, interval)

# Optionally, you can combine the raster series
combined_raster = Rasters.combine(raster_series, Ti)[bounds...]
