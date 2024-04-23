using Pkg

Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")

using Rasters, ArchGDAL, NCDatasets, GeoData, Plots, Shapefile
using PlotShapefiles, HypothesisTests, Distributions
using GeoInterface
using RasterDataSources # This is needed to download common data sources
using IntervalSets
using DimensionalData, CoordinateTransformations

meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"

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
