bounds = (X(Rasters.Between(-10, 4)), Y(Rasters.Between(35.5, 44)))

# Step 2: Load your initial and final rasters (replace with your actual data paths)
bio5_initial_raster = Raster("CHELSA/BioClim/bio5_2010.tif", lazy = true)[bounds...]
bio5_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio5_2100.tif", lazy = true)[bounds...]

bio6_initial_raster = Raster("CHELSA/BioClim/bio6_2010.tif", lazy = true)[bounds...]
bio6_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio6_2100.tif", lazy = true)[bounds...]

bio12_initial_raster = Raster("CHELSA/BioClim/bio12_2010.tif", lazy = true)[bounds...]
bio12_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio12_2100.tif", lazy = true)[bounds...]

# Step 3: Resample all rasters to match utmraster
bio5_initial_raster = Rasters.resample(bio5_initial_raster, to = utmraster)
bio5_final_raster = Rasters.resample(bio5_final_raster, to = utmraster)

bio6_initial_raster = Rasters.resample(bio6_initial_raster, to = utmraster)
bio6_final_raster = Rasters.resample(bio6_final_raster, to = utmraster)

bio12_initial_raster = Rasters.resample(bio12_initial_raster, to = utmraster)
bio12_final_raster = Rasters.resample(bio12_final_raster, to = utmraster)

# Step 4: Apply the scaling factors and offsets
# For bio5 and bio6: Temperature in °C = (Value * Scale) + Offset
# Scale = 0.1, Offset = -273.15

bio5_initial_raster = map(x -> (Float32(x) * 0.1) + (-273.15), bio5_initial_raster)
bio5_final_raster = map(x -> (Float32(x) * 0.1) + (-273.15), bio5_final_raster)

bio6_initial_raster = Rasters.map(x -> (Float32(x) * 0.1) + (-273.15), bio6_initial_raster)
bio6_final_raster = Rasters.map(x -> (Float32(x) * 0.1) + (-273.15), bio6_final_raster)

# For bio12: Precipitation in mm/year = (Value * Scale) + Offset
# Scale = 0.1, Offset = 0

bio12_initial_raster = Rasters.map(x -> Float32(x) * 0.1, bio12_initial_raster)
bio12_final_raster = Rasters.map(x -> Float32(x) * 0.1, bio12_final_raster)

# Step 5: Define the time span and dates
t_start = Date(2010, 1, 1)
t_end = Date(2100, 1, 1)
tspan = t_start:tspan_intervals:t_end

# Step 6: Calculate the interpolation factors
total_days = Dates.value(t_end - t_start)
factors = [Dates.value(d - t_start) / total_days for d in tspan]

# Step 7: Generate the interpolated rasters for each variable
bio5_raster_list = Raster[]
bio6_raster_list = Raster[]
bio12_raster_list = Raster[]

for f in factors
    # Interpolate bio5 raster
    bio5_interp = (1 - f) .* bio5_initial_raster + f .* (static ? bio5_initial_raster : bio5_final_raster)
    push!(bio5_raster_list, bio5_interp)
    
    # Interpolate bio6 raster
    bio6_interp = (1 - f) .* bio6_initial_raster + f .* (static ? bio6_initial_raster : bio6_final_raster)
    push!(bio6_raster_list, bio6_interp)
    
    # Interpolate bio12 raster
    bio12_interp = (1 - f) .* bio12_initial_raster + f .* (static ? bio12_initial_raster : bio12_final_raster)
    push!(bio12_raster_list, bio12_interp)
end

# Step 8: Create RasterSeries for each variable
bio5_raster_series = RasterSeries(bio5_raster_list, Ti(tspan))
bio6_raster_series = RasterSeries(bio6_raster_list, Ti(tspan))
bio12_raster_series = RasterSeries(bio12_raster_list, Ti(tspan))

########## PLOTTING THE DATA ############
# q = min(minimum(bio5_raster_series[Ti(1)]), minimum(bio5_raster_series[Ti(end)]))
# w = max(maximum(bio5_raster_series[Ti(1)]), maximum(bio5_raster_series[Ti(end)]))
# e = min(minimum(bio12_raster_series[Ti(1)]), minimum(bio12_raster_series[Ti(end)]))
# r = max(maximum(bio12_raster_series[Ti(1)]), maximum(bio12_raster_series[Ti(end)]))
# # Create a figure with 2 axes
# fig = Figure(resolution = (800, 400))

# # Plot the first raster in the bio5 series
# ax1 = Axis(fig[1, 1])
# ax2 = Axis(fig[1, 2])
# heatmap!(ax1, bio5_raster_series[Ti(25)], colorrange = (q, w), colormap = :thermal)
# heatmap!(ax2, bio5_raster_series[Ti(75)], colorrange = (q, w), colormap = :thermal)

# # Display the figure
# fig