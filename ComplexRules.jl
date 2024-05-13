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
using Dates, Distributions
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{256, T}
    b::T
    
    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{256, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{3, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

function growth(abundance::AbstractFloat, self_regulation::AbstractFloat, K::AbstractFloat)
    # Assuming growth is a simple function of abundance, self-regulation, and carrying capacity (K)
    return self_regulation * K * abundance * (1 - abundance / K)
end
function growth_mystruct(abundance::MyStructs256, self_regulation::AbstractFloat, K::AbstractFloat)
    return MyStructs256(self_regulation .* K .* abundance.a .* (1.0 .- (abundance.a ./ K)))
end

init = MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0))
@time growth_mystruct(init, 0.01, 1000.0)


A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}
full_A_matrix = A_matrix_list[3]

full_IM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
full_IM = fullIM_list[3][10]
growth()
function ode{MyStructs256}(abundances, self_regulation, K, num_steps = nothing)
    
    abundances = data[I...]
    K = npp_array[I...]
    
    growth = abundances .+ growth(abundances, K, self_regulation)
    return MyStructs256(abundances)
end

function trophic(abundances, A_matrix)
    
end

