using Pkg
# Desktop PC
# Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")
# cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 
# Laptop
Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\GitHub\\simBio")
cd("C:\\Users\\nicol\\OneDrive\\PhD\\GitHub\\simBio")


using Rasters, ArchGDAL, Plots

using CSV, DataFrames
# meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop
# using PlotShapefiles, Colors

spain_vector = Shapefile.Table(joinpath(meta_path, "ATLAS_data", "UTMgrid_IberianPeninsula", "UTM_GRID.shp")) |> DataFrame 
# plot(spain_vector, color = :white, fill = true, alpha = 0.5)

Amph = CSV.read(joinpath(meta_path, "Atlas_data", "DB_Amphibians_IP.txt"), delim='\t', DataFrame) |> DataFrame

spain_vector_with_amph = leftjoin(spain_vector, Amph, on = :UTMCODE)
spain_vector_with_amph.rowsum = sum.(eachrow(spain_vector_with_amph[:, 3:end])) 

rowsum = spain_vector_with_amph.rowsum

# Plot the shapefile with rowsum values mapped to colors
plot(spain_vector_with_amph.geometry, color = spain_vector_with_amph.rowsum, fill = true, alpha = 0.5)
plot(spain_vector_with_amph.geometry, color = spain_vector_with_amph.rowsum, seriestype = :shape, fillalpha = 0.5, legend = false)
centroids = Shapefile.Table(joinpath(meta_path, "iberian_centroids_shp", "iberian_centroids.shp")) 
geoms = Shapefile.shapes(centroids)

bioclim_5 = Raster(joinpath(meta_path, "Rasters", "iberian_temperature.tif"))
plot(bioclim_5)

plot(centroids, color = :black, alpha = 0.1)
