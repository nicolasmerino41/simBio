europe = RS.Raster("Rasters\\Europe_bio1.tif")
# MK  .plot(europe, colormap = custom_palette)

europe_resampled = RS.resample(europe, to = utmraster)
MK.plot(europe_resampled, colormap = custom_palette)
RS.write("Rasters\\Europe_bio1_resampled.tif", europe_resampled)

europe = RS.Raster("masked_raster_resampled.tif")
MK.plot(europe)
hola = reverse(parent(europe), dims = 2)
europe[13:22, 1:2] .= NaN
# europe[152, 351:352] .= NaN
europe_trial = deepcopy(europe)
europe_trial[457:470, 139:155] .= 50000
hola[457:470, 139:155] .= 50000
MK.plot(europe_trial)
MK.plot(hola)