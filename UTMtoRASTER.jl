using Pkg
# Desktop PC
# Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

# meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
spain_vector = Shapefile.Table(joinpath(meta_path, "ATLAS_data", "UTMgrid_IberianPeninsula", "UTM_GRID.shp"))  
Plots.plot(spain_vector, color = :black, fill = true, alpha = 0.5)

bbox = SH.rect(spain_vector)

utmraster = Raster("C:/Users/MM-1/OneDrive/PhD/JuliaSimulation/simBio/updated_utmraster.tif") 
pr = parent(utmraster)
MK.plot(utmraster)
RS.values(utmraster)

perro = CSV.File("Species_spain_df.csv") |> DataFrame
perro = select(perro, Not([:UTMCODE, :Value, :sum, :ID]), Not(in(names(perro), Ref(unique_species_in_web))))

variables = perro[!, [:UTMCODE, :Value, :ID, :sum]]
variables.Value = map(x -> x == "NA" ? 0.0 : parse(Float64, x), variables.Value)

rest = perro[!, unique_species_in_web]
for i in axes(perro_cropped_matrix, 1), j in axes(perro_cropped_matrix, 2)
    if perro_cropped_matrix[i, j] == 0
        perro_cropped_matrix[i, j] = 0.0
    elseif perro_cropped_matrix[i, j] == 1
        perro_cropped_matrix[i, j] = 1.0
    end
end

perro_cropped = hcat(variables, rest, makeunique=true)
perro_cropped_matrix = Matrix(perro_cropped)
utmraster_da = DimArray(utmraster)
D = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))


for i in 1:size(perro_cropped, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        # println(utmraster_da[j])
        if Float32(perro_cropped.Value[i]) == Float32(utmraster_da[j])
            D[j] = MyStructs256(SVector{256, Float64}(perro_cropped_matrix[i, 5:260]))
        end
    end
end

MK.plot(D, color = :black, fill = true, alpha = 0.5)

D_with_abundances = deepcopy(D)
for row in axes(D, 1), col in axes(D, 2)
    if D[row, col] != MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        new_a = SVector{256, Float64}([D[row, col].a[i] != 0.0 ? rand(10:100) : D[row, col].a[i] for i in 1:256])
        D_with_abundances[row, col] = MyStructs256(new_a)
    end
end

