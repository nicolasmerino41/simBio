using Pkg

Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 

meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"

ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio"

using Dates, Plots, DataFrames

lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti)) 
# this generates random numbers with the dimensions given

lon = Rasters.lookup(ras, X) # if X is longitude
lat = Rasters.lookup(ras, Y) # if Y is latitude

ras[Ti=At(DateTime(2001))]

using Rasters, RasterDataSources
using ArchGDAL
# using DynamicGrids

world = Raster(WorldClim{BioClim}, 5)
iberian_peninsula = view(world, X(-10 .. 3), Y(35 .. 44))
plot(iberian_peninsula)

world[X(-10 .. 3)]
world[X(-10 .. 3), Y(35 .. 44)]

shrubs_lazy_raster = Raster(EarthEnv{LandCover}, :shrubs; lazy = true)
plot(shrubs_lazy_raster)

typeof(parent(shrubs_lazy_raster)) # LazyArray{Float32, 3}
broadcasted_lazy_raster = replace_missing(shrubs_lazy_raster) .*0.01
plot(broadcasted_lazy_raster)

biodf = DataFrame(replace_missing(iberian_peninsula)) 
filter(x -> !ismissing(x.bio5), biodf)

using Shapefile, GeoInterface
spain_vector = Shapefile.Table(joinpath(meta_path, "ATLAS_data", "UTMgrid_IberianPeninsula", "UTM_GRID.shp")) |> DataFrame

plot(spain_vector.geometry)
utm_rasterised = Rasters.rasterize(spain_vector.geometry[1:5927], res=0.1, missingval = 0,
fill = 1, boundary = :touches) 

plot(utm_rasterised)
caca = spain_vector.geometry[1]
plot(caca)

polygon_list = []

# Extracting the centroids of each geometry in the spain_vector DataFrame
spain_centroids = [centroid(geom) for geom in spain_vector.geometry]

