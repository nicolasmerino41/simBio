library(readr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)
library(geodata)
library(sdm)
library(usdm)
library(terra)
library(raster)
library(sf)
library(progress)
library(tidyterra)

setwd("C:/Users/MM-1/OneDrive/PhD/JuliaSimulation/simBio")
##### ATLAS DATA #####
# Loading species data from ATLAS data
iberian_interact_matrix = read.csv("Objects/iberian_interact_matrix.csv")[,2:257]
Species_spain <- vect("Objects/Species_spain_vector.shp")
species_names = gsub("\\.", " ", colnames(iberian_interact_matrix))
names(Species_spain) = species_names 
Species_spain_df <- as.data.frame(Species_spain) # Transform the vector into a df

spanish_points <- centroids(Species_spain)
crs(spanish_points) = "+proj=longlat +datum=WGS84 +no_defs"
##### BIO VARIABLES #####
species_names = species_names
bio5 <- rast("wc2.1_5m/wc2.1_5m_bio_5.tif")
bio6 <- rast("wc2.1_5m/wc2.1_5m_bio_6.tif")
bio12 <- rast("wc2.1_5m/wc2.1_5m_bio_12.tif")
bio5_cropped <- crop(bio5, Species_spain)
bio6_cropped <- crop(bio6, Species_spain)
bio12_cropped <- crop(bio12, Species_spain)

preds = brick(c(bio5_cropped, bio6_cropped, bio12_cropped))

# ex = raster::extract(preds, as(spanish_points, "Spatial"), ID = F)
# species_extracts_list = list()
# for(i in names(Species_spain)) {
#   df <- cbind(as.data.frame(ex), Species_spain_df[,i])
#   colnames(df) = c("bio5", "bio6", "bio12", "presence")
#   species_extracts_list[[i]] <- df
# }
# coordi = geom(spanish_points)[,3:4]
# coordi_for_each_species = list()
# i = species_names[[1]]
# 
# for(i in species_names) {
#   d = gsub(" ", ".", i)
#   df = cbind(coordi, species_extracts_list[[d]][,4])
#   colnames(df) = c("x", "y", "species")
#   df = as.data.frame(df)
#   coordinates(df) = c("x", "y")
#   coordi_for_each_species[[i]] = df
# }
# 
# sdmData_list = list()
# for(i in species_names) {
#   sdmData_list[[i]] = sdmData(formula = species ~ .,
#                               train = coordi_for_each_species[[i]],
#                               predictors = preds)
# }
# 
# sdm_list = list()
# for (i in species_names) {
#   sdm_list[[i]] = sdm(species ~ ., data = sdmData_list[[i]],
#                       methods = c("glm", "brt", "rf"),
#                       replication = "boot", n = 1,
#                       parallelSettings = list(ncores = 6, method = "parallel")
#   )
#   # Using match() to find the position
#   position <- match(i, species_names)
#   cat("sdm", i, "done, nº", position, "\n")
# }

niches <- data.frame(Species = species_names,
                     mean_bio5 = rep(0.0, length(species_names)),
                     sd_bio5 = rep(0.0, length(species_names)),
                     mean_bio6 = rep(0.0, length(species_names)),
                     sd_bio6 = rep(0.0, length(species_names)),
                     mean_bio12 = rep(0.0, length(species_names)),
                     sd_bio12 = rep(0.0, length(species_names)))

# Function to rescale values to [0, 1]
rescale_to_unit <- function(valid_values, max_value, min_value) {
  (valid_values - min_value) / (max_value - min_value)
}

for (species in species_names) {
  file_name <- paste0("Ensembles/", species, ".tif")
  
  # Perform ensemble modeling
  en <- raster(file_name)
  species = gsub(" ", ".", species)
  # Extract probability and bio variables only for presence points
  absence_data <- spanish_points[spanish_points[[species]] == 0, ]
  probs <- extract(en, as(absence_data, "Spatial"))
  
  # Find the highest probability with absence
  max_prob_absence <- max(probs)
  
  thresholded_ensemble = en 
  thresholded_ensemble[]<-ifelse(thresholded_ensemble[] < max_prob_absence | is.na(thresholded_ensemble[]),-0.01, thresholded_ensemble[])
  # Extract values that are not NA or negative for scaling
  all_values <- values(thresholded_ensemble)  # Extract all values
  valid_indices <- which(all_values >= 0.0)  # Indices of values that are >= -0.01
  valid_values <- all_values[valid_indices]  # Extract those values
  
  # Calculate min and max from these valid values
  max_value <- max(valid_values)
  min_value <- min(valid_values)
  
  valid_cells <- which(values(thresholded_ensemble) > 0.0)
  valid_coords <- xyFromCell(thresholded_ensemble, valid_cells)
  
  # Calculate means and standard deviations
  val5 = extract(bio5_cropped, valid_coords)
  mean_bio5 <- mean(val5$wc2.1_5m_bio_5, na.rm = TRUE)
  sd_bio5 <- sd(val5$wc2.1_5m_bio_5, na.rm = TRUE)
  val6 = extract(bio6_cropped, valid_coords)
  mean_bio6 <- mean(val6$wc2.1_5m_bio_6, na.rm = TRUE)
  sd_bio6 <- sd(val6$wc2.1_5m_bio_6, na.rm = TRUE)
  val12 = extract(bio12_cropped, valid_coords)
  mean_bio12 <- mean(val12$wc2.1_5m_bio_12, na.rm = TRUE)
  sd_bio12 <- sd(val12$wc2.1_5m_bio_12, na.rm = TRUE)
  
  species = gsub("\\.", " ", species)
  # Assign the calculated values to the niches data frame
  niches[niches$Species == species, "mean_bio5"] <- mean_bio5
  niches[niches$Species == species, "sd_bio5"] <- sd_bio5
  niches[niches$Species == species, "mean_bio6"] <- mean_bio6
  niches[niches$Species == species, "sd_bio6"] <- sd_bio6
  niches[niches$Species == species, "mean_bio12"] <- mean_bio12
  niches[niches$Species == species, "sd_bio12"] <- sd_bio12
  
  # Using match() to find the position
  position <- match(species, species_names)
  cat("Ensemble", species, "done, nº", position, "\n")
}

write.csv2(niches, "iberian_species_niches_withVeryStrictNiche.csv", row.names = F)
















