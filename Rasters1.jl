using Pkg

Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 

ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio"

using Rasters, Dates, Plots

lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti)) 
# this generates random numbers with the dimensions given

lon = Rasters.lookup(ras, X) # if X is longitude
lat = Rasters.lookup(ras, Y) # if Y is latitude

ras[Ti=At(DateTime(2001))]

using Rasters, RasterDataSources
using ArchGDAL
using DynamicGrids
using Plots

world = Raster(WorldClim{BioClim}, 5)
iberian_peninsula = view(world, X(-10 .. 3), Y(36 .. 43))
plot(iberian_peninsula)

world[X(-10 .. 3)]
world[X(-10 .. 3), Y(36 .. 43)]

Between(10, 20)