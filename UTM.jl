using Pkg
# Desktop PC
# Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 
# cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")
# Laptop
Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

# meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using Rasters, NCDatasets, Plots, Shapefile, ArchGDAL
using CSV, DataFrames
#################################################################################################

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

# Create a DataFrame to hold the X and Y coordinates of the first two points from each geometry
geom_points_df = DataFrame(X = Float64[], Y = Float64[])

# Iterate over the geometries and extract the first two points
for geom in geoms
    push!(geom_points_df, [geom.x, geom.y])
end

# Plot geom_points_df with column X on the x-axis and column Y on the y-axis
plot(geom_points_df.X, geom_points_df.Y, seriestype = :scatter, xlabel = "X", ylabel = "Y")

bioclim_5 = Raster(joinpath(meta_path, "Rasters", "iberian_temperature.tif"))
plot(bioclim_5, title = "")
bioclim_5 = bioclim_5[X(-10 .. 40), Y(36 .. 46)]

centroids

p = plot(bioclim_5, title = "")
plot!(p, [geom.x for geom in geoms], [geom.y for geom in geoms], 
legend = false, color = :black, alpha = 0.5, markersize = 2, seriestype = :scatter)

plot!(p, [geom.x for geom in geoms], [geom.y for geom in geoms], legend = false, color = :black, 
alpha = 0.5, markersize = 2, seriestype = :scatter) # They're the same, just two different ways to do it




