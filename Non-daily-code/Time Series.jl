# aux=(; climate_aux_name=climatedata)

# ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio"

# raster = getraster(CHELSA{BioClim}, 1)

# series = Raster[]
# for i in 1:5
#     raster = Raster(getraster(CHELSA{BioClim}, 2))
#     raster_resampled = RS.resample(raster, to = utmraster)
#     push!(series, raster_resampled) 
# end
# MK.plot(series[1])
# series1 = RS.RasterSeries(series, (125, 76))  # Assuming these are the dimensions of the rasters

path = "CHELSA/BioClim/"
filelist = readdir(path)
# MK.plot(Raster(path .* filelist[1], lazy = true)[bounds...], colormap = :thermal)
dates = collect(Date(2023,1):Month(1):Date(2025, 6))
bounds = X(Rasters.Between(-10, 4)), Y(Rasters.Between(35.5, 44)) # bounding box around IP

raster_list = [RS.resample(Raster(i, lazy = true)[bounds...], to = utmraster) for i in path .* filelist]
raster_series = RasterSeries(raster_list, Ti(dates))
combined_raster = Rasters.combine(raster_series, Ti)[bounds...]

# MK.plot(raster_series[Ti(1)], colormap = :thermal)
