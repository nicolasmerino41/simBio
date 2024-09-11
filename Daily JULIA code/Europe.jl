europe = RS.Raster("Rasters\\Europe_bio1.tif")
MK.plot(europe, colormap = custom_palette)

europe_resampled = RS.resample(europe, to = utmraster)
MK.plot(europe_resampled, colormap = custom_palette)
RS.write("Rasters\\Europe_bio1_resampled.tif", europe_resampled)
