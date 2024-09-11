setwd("C:/Users/nicol/OneDrive/PhD/JuliaSimulation/simBio")
library(raster)
library(terra)

europe = rast("Rasters/Europe_bio1.tif")

europe = project(utmraster, europe, )
utmraster = rast("Rasters/utmraster.tif")

europe_resampled = raster::resample(europe, utmraster)
plot(europe_resampled)

spain_res <- res(utmraster)
# Resample the Europe raster using the resolution of the Spain raster
resampled_europe_raster <- resample(europe, utmraster, method="bilinear")
