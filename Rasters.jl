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

# Load the raster data
raster = ArchGDAL.read("Rasters\\iberian_temperature.tif")
raster_band = ArchGDAL.getband(raster, 1)

plot(raster_band)

spain_vector = open_shapefile("ATLAS_data/UTMgrid_IberianPeninsula/UTM_GRID.shp")
spain = open_shapefile("spain.shp")
PlotShapefiles.plotshape(spain)

cucu = Shapefile.open("spain.shp")
plot(cucu)
cucu = Shapefile.Table("spain.shp")

pepe = read("spain.shp")
typeof(pepe)
plot(pepe) # This runs but not properly
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"
# bio5_path = "/wc2.1_10m_bio5.tif"
# bio5_path = RasterDataSources.getraster(WorldClim{BioClim}, 5)
bioclim_paths = RasterDataSources.getraster(WorldClim{BioClim}, (5,7,8,12))

bioclim_stack = RasterStack(WorldClim{BioClim}, (5, 7, 8, 12), res="5m")
bioclim_stack = RasterStack(bioclim_paths)

plot(bioclim_stack)
bioclim_5 = Raster(WorldClim{BioClim}, 5, res="5m")
spain_bioclim = bioclim_stack[X(-9 .. 3), Y(36 .. 43)]
spain_bioclim5 = bioclim_stack[X(10 .. 30), Y(36 .. 43)]
madagascar_bioclim = bioclim_stack[X(-180.0 .. -166.5), Y(-27.0 .. -11.0)]
cropped_raster_stack = Rasters.crop(bioclim_5, raster, Touches = true)
plot(spain_bioclim)
nz_bounds = X(At(165 .. 178)), Y(At(-50 .. -32))
nz_bioclim = bioclim_5[nz_bounds...]

plot(bioclim_5)
dims(bioclim_5)
parent(bioclim_5)
crs(bioclim_5)
messed_up_raster = rotl90(bioclim_5)

plot(messed_up_raster)
const AG = ArchGDAL
dataset = ArchGDAL.read("world.tif")
bioclim_5_AG = AG.read("WorldClim/BioClim/wc2.1_5m_bio_5.tif")
bioclim_5_band = AG.getband(bioclim_5_AG, 1)
plot(bioclim_5_band)
plot(bioclim_stack)
band = ArchGDAL.getband(dataset, 1)
bandd = AG.getband(bioclim_5_AG, 1)
banddd = rotl90(bandd)
bandddd = permutedims(banddd)
bioclim_stack_AG = AG.read(bioclim_paths)
bioclim_paths

plot(bandddd)
bandd
values(bandd)
maximum(skipmissing(bandd))
















using Rasters, RasterDataSources, ArchGDAL, Plots
using CairoMakie
CairoMakie.activate!()
# st = RasterStack(WorldClim{Climate}; month=1);
st_paths = RasterDataSources.getraster(WorldClim{Climate}; month=1)
st = RasterStack(st_paths)

africa = st[X(-20 .. 60), Y(-40 .. 35)]
a = plot(africa)

aus = st[X(100.0 .. 160.0), Y(-50.0 .. -10.0)]
b = plot(aus)

# Combine with mosaic
mos = mosaic(first, aus, africa)
c = plot(mos)

savefig(a, "docs/build/mosaic_example_africa.png")
savefig(b, "docs/build/mosaic_example_aus.png")
savefig(c, "docs/build/mosaic_example_combined.png")
nothing
# output

using Rasters, RasterDataSources
const RS = Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{BioClim}, 5, res="5m")
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04))


### Load and plot the file
# Set the environment variable for raster data sources to a default path
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio"

bioclim_5 = Raster(WorldClim{BioClim}, 5, res="5m")

iberian_temp = Raster(joinpath(meta_path, "Rasters", "iberian_temperature.tif"))
iberian_temp = iberian_temp[X(-10 .. 4), Y(35 .. 44)]

plot(iberian_temp)

cucu = Rasters.crop(bioclim_5, to = iberian_temp)

pepe = DataFrame(iberian_temp)
non_missing_line_count = count(isequal(true), .!ismissing.(pepe))

plot(iberian_temp)
using DataFrames

countries = Shapefile.Table("C:\\Users\\MM-1\\Downloads\\country_shapes\\country_shapes.shp") |> DataFrame
plot(countries.geometry, color = :black)

countries
Portugal_borders = filter(x -> x.cntry_name in ["Portugal"], countries).geometry[1]
Spain_borders = filter(x -> x.cntry_name in ["Spain"], countries).geometry[1]

portugal_raster = Rasters.rasterize(Portugal_borders, res = 0.1, missingval = 0, fill=1, boundary = :touches)
spain_raster = Rasters.rasterize(Spain_borders, res = 0.1, missingval = 0, fill=1, boundary = :touches)

portugal_raster_cropped = Rasters.crop(portugal_raster, to = iberian_temp)
spain_raster_cropped = Rasters.crop(spain_raster, to = iberian_temp)

IP_raster = mosaic(portugal_raster_cropped, spain_raster_cropped)
plot(spain_raster_cropped)

plot(IP_borders.geometry)
spain_border_rasterised

utm = Shapefile.Table(joinpath(meta_path, "ATLAS_data", "UTMgrid_IberianPeninsula", "UTM_GRID.shp")) |> DataFrame
plot(utm.geometry)