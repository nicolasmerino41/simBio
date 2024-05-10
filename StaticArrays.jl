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

struct MyStructss{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{3, T}
    b::T
    
    # Custom constructor
    function MyStructss(a::SVector{3, T}) where T <: AbstractFloat
        new{T}(a, sum(a))  # Correctly specify the type parameter with `new`
    end
end

# Ensure MyStructss is immutable and behaves correctly in grid operations
Base.zero(::Type{MyStructss{T}}) where T = MyStructss(SVector{3, T}(zero(T), zero(T), zero(T)))
Base.oneunit(::Type{MyStructss{T}}) where T = MyStructss(SVector{3, T}(oneunit(T), oneunit(T), oneunit(T)))

Base.isless(a::MyStructss, b::MyStructss) = isless(a.b, b.b)  # Comparison based on the total abundance 'b'

# Arithmetic operations should also account for immutability
Base.:*(x::MyStructss, scalar::Real) = MyStructss(x.a .* scalar)  # Element-wise multiplication for vector
Base.:/(x::MyStructss, scalar::Real) = MyStructss(x.a ./ scalar)  # Element-wise division for vector
Base.:+(x1::MyStructss, x2::MyStructss) = MyStructss(x1.a .+ x2.a)  # Element-wise addition for vectors
Base.:-(x1::MyStructss, x2::MyStructss) = MyStructss(x1.a .- x2.a)  # Element-wise subtraction for vectors

# Visualization function to reflect total abundance 'b' in grayscale
DynamicGrids.to_rgb(::ObjectScheme, obj::MyStructss) = begin
    intensity = clamp(obj.b / 1000, 0.0, 1.0)  # Normalize and clamp 'b' to the range [0,1]
    RGB(intensity, intensity, intensity)       # Create a grayscale color based on intensity
end


rule = Cell{:grid1}() do data, state, I
        2state
end
init = (grid1 = MyStructss(SVector(150.0, 200.0, 100.0)))  # Initialize with three species automatically calculating 'b'
output = ArrayOutput(init; tspan=1:3)ERROR: The constructor for MyStructss{Float64}(::Float64, ::Float64) 






# Example to create an instance of MyStruct
abundances = SVector(150.0, 200.0, 100.0)  # Example abundances for three species
abundances
cell = MyStructss(abundances)  # b is automatically set to the sum of abundances
cell.b

using StaticArrays
using Base: @kwdef

using StaticArrays
using Base: @kwdef

@kwdef struct MyStructss{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{3, T}
    b::T

    # Primary constructor
    function MyStructss(a::SVector{3, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor for direct initialization
    function MyStructss(a::SVector{3, T}, b::T) where T <: AbstractFloat
        new{T}(a, b)
    end

    # Debug constructor for incorrect calls
    function MyStructss(a::Float64, b::Float64)
        @warn "Incorrect constructor call for MyStructss with (Float64, Float64)"
        new{Float64}(SVector{3, Float64}(a, b, 0), a + b)
    end
end

# Define operations
import Base: zero, oneunit, isless, *, +, -, size, getindex

size(::MyStructss) = (2,)
getindex(x::MyStructss, i::Int) = i == 1 ? x.a : x.b

zero(::Type{MyStructss{T}}) where {T <: AbstractFloat} = MyStructss(SVector{3, T}(zero(T), zero(T), zero(T)))
oneunit(::Type{MyStructss{T}}) where {T <: AbstractFloat} = MyStructss(SVector{3, T}(one(T), one(T), one(T)))

function *(x::MyStructss, scalar::Real)
    new_a = x.a * scalar
    MyStructss(new_a)  # Automatically recalculates 'b'
end

function +(x1::MyStructss, x2::MyStructss)
    new_a = x1.a + x2.a
    MyStructss(new_a)  # Automatically recalculates 'b'
end

function -(x1::MyStructss, x2::MyStructss)
    new_a = x1.a - x2.a
    MyStructss(new_a)  # Automatically recalculates 'b'
end


# Simulation example
rule = Cell{:grid1}() do data, state, I
    2 * state  # Applies multiplication, 'b' recalculates automatically
end

# Initialize with direct values for 'a' and 'b' (only if absolutely necessary)
init = (grid1 = MyStructss(SVector(150.0, 200.0, 100.0)))  # Direct initialization
output = ArrayOutput(init; tspan=1:3)

