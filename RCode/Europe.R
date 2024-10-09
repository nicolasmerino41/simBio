library(raster)
library(terra)
library(sf)
library(rnaturalearth)
library(lwgeom)
europe_bio1 = rast("Rasters/Europe_bio1.tif")
utmraster = rast("Rasters/utmraster.tif")

# Get Europe data
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Manually remove island countries/territories by their ISO_A3 codes
# These codes are based on your data input
islands_to_remove <- c("ISL", "MLT", "CYP", "IMN", "JEY", "GGY", "ALA", "FRO", "VAT", "SMR", "MCO")

# Filter the europe dataframe to exclude these island countries
europe_mainland <- europe[!europe$iso_a3 %in% islands_to_remove, ]

EEA = vect("Europe/euro-global-map-shp/NUTS_3.shp")
EEA_df = as.data.frame(EEA)
unique(EEA$ICC)
countries_to_remove <- c("GB", "IE", "ND", "MT", "IS")

# Filter out these country codes from the EEA object
EEA_filtered <- EEA[!EEA$ICC %in% countries_to_remove, ]
extent_bbox <- ext(-10, 50, 32, 75)
# Crop the SpatVector to the specified bounding box
EEA_cropped <- crop(EEA_filtered, extent_bbox)
EEA_cropped_df = as.data.frame(EEA_cropped)
plot(EEA_cropped)
eacrs = crs(EEA, proj = TRUE)
crs(nuts2_trial_projected, proj = TRUE)
nuts2_trial_projected = project(nuts2_trial, eacrs)
####### NUTS
nuts = vect("NUTS/NUTS_RG_20M_2024_3035.shp")
plot(nuts)
nuts_df = as.data.frame(nuts)
nuts2 = nuts[nuts_df$LEVL_CODE %in% c(1,2)]
nuts2_df = as.data.frame(nuts2)
nuts2_trial <- nuts2[!(nuts2$NAME_LATN %in% island_regions | nuts2$CNTR_CODE %in% c("IE", "TR")), ]
nuts2_trial_projected = project(nuts2_trial, eacrs)
nuts2_trial = crop(nuts2_trial_projected, extent_bbox)
mediterranean = vect("Mediterranean/World_Seas_IHO_v3.shp")
mediterranean = mediterranean[mediterranean$ID %in% c("28A", "28B"),]
mediterranean_df = as.data.frame(mediterranean)
plot(mediterranean)
plot(nuts2_trial)
# Convert to sf object if necessary
# Replace 'geometry_column' with the actual geometry column name if different
nuts2_final <- st_as_sf(nuts2_trial, wkt = "geometry_column")
# Inspect the structure of nuts2_final
str(nuts2_final)

# View the first few rows
head(nuts2_final)


island_regions <- c(
  "K\xfdpros",
  "Voreio Aigaio",
  "Notio Aigaio",
  "Kriti",
  "Illes Balears",
  "Canarias",
  "\xcdsland",
  "Sicilia",
  "Sardegna",
  "Malta",
  "\xc5land",
  "Corse",
  "Guadeloupe",
  "Martinique",
  "La R\xe9union",
  "Mayotte",
  "Ionia Nisia",
  "Svalbard og Jan Mayen",
  "Regi\xe3o Aut\xf3noma dos A\xe7ores",
  "Regi\xe3o Aut\xf3noma da Madeira",
  "Isole",
  "RUP FR ? R\xe9gions Ultrap\xe9riph\xe9riques Fran\xe7aises"
)

# Get CRS of raster
crs_raster <- crs(europe_bio1)

# Get CRS of vector
crs_vector <- st_crs(nuts2_final)$wkt

# Compare CRS
if (crs_vector != crs_raster) {
  # Transform vector to match raster CRS
  nuts2_final_transformed <- st_transform(nuts2_final, crs = crs_raster)
  message("Vector data transformed to match raster CRS.")
} else {
  nuts2_final_transformed <- nuts2_final
  message("CRS of vector and raster match. No transformation needed.")
}
# Convert sf object to SpatVector
land_vector <- vect(nuts2_final_transformed)

# Inspect the SpatVector
print(land_vector)
plot(land_vector, border = "red", lwd = 2)

# Perform the mask operation
masked_raster <- mask(europe_bio1, land_vector)

# Inspect the masked raster
print(masked_raster)
plot(masked_raster, main = "Masked Raster (Land Areas Only)")


####################################
####################################
# 1. Load your rasters
# Replace with your actual file paths
masked_raster <- rast("Rasters/Europe_bio1.tif")  # Replace with your masked_raster path
utmraster <- rast("Rasters/utmraster.tif")          # Replace with your utmraster path

# 2. Inspect CRS and resolutions
print("CRS of masked_raster:")
print(crs(masked_raster))  # Should show EPSG:4326

print("CRS of utmraster:")
print(crs(utmraster))       # Currently empty

print("Resolution of masked_raster:")
print(res(masked_raster))  # 0.1666667, 0.1666667

print("Resolution of utmraster:")
print(res(utmraster))       # 9927.768, 9927.768

# 3. Assign CRS to utmraster if known (optional)
# Example: If utmraster is supposed to be in UTM Zone 33N (EPSG:32633)
# crs(utmraster) <- "EPSG:32633"  # Uncomment and replace with correct EPSG if applicable

# Since you want to keep masked_raster's CRS, we'll proceed without assigning a CRS to utmraster

# 4. Adjust masked_raster's resolution
current_res <- res(masked_raster)
print(paste("Original Resolution:", current_res[1], current_res[2]))  # 0.1666667, 0.1666667

# Desired resolution: ~0.0833333 degrees (~9.26 km at equator)
factor <- 2

# Disaggregate masked_raster
masked_raster_new <- disagg(masked_raster, fact = factor, method = "near")

# Verify the new resolution
new_res <- res(masked_raster_new)
print(paste("New Resolution:", new_res[1], new_res[2]))  # ~0.0833333, ~0.0833333

# Inspect new dimensions
print(dim(masked_raster_new))  # Should be roughly double in each dimension (444 x 726)

# 5. Visualize the original and resampled rasters
par(mfrow = c(1, 2))  # Set up a 1x2 plotting area

# Plot original masked_raster
plot(masked_raster, main = "Original Masked Raster (0.1666667°)")

# Plot resampled masked_raster_new
plot(masked_raster_new, main = "Resampled Masked Raster (~0.0833333°)")
masked_raster_new = mask(masked_raster_new, land_vector)
extent_bbox <- ext(-10, 32, 32, 75)
masked_raster_new = crop(masked_raster_new, extent_bbox)
par(mfrow = c(1, 1))  # Reset plotting area
plot(masked_raster_new)
# 6. (Optional) Save the resampled raster
output_path <- "masked_raster_resampled.tif"  # Replace with your desired path
writeRaster(masked_raster_new, filename = output_path, overwrite = TRUE)
message("Resampled masked raster saved to ", output_path)


