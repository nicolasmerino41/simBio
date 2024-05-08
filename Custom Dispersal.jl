using Pkg
# Desktop PC
Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
# Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
# meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
#################################################################################################
####################################################################
######## Non-weighted Dispersal
struct CustomKernel <: KernelFormulation
    α::Float64
end
(kernel::CustomKernel)(distance) = exp(-distance / (2 * kernel.α^2))
# Define your custom dispersal kernel struct
Base.@kwdef struct CustomDispersalKernell{P} <: KernelFormulation
    α::P = Param(1.0, bounds=(0.0, Inf))  # Let's ensure α is positive
end
(f::CustomDispersalKernell)(d) = CustomKernel(f.α)

####### WEIGHTED DISPERSAL KERNEL
# Define the dispersal kernel considering two properties from different rasters
struct CustomKernel <: KernelFormulation
    α::Float64  # Dispersion coefficient
end

# Implement the dispersal kernel function that uses abundance and carrying capacity
function (kernel::CustomKernel)(abundance, carrying_capacity, distance)
    proportion = abundance / max(carrying_capacity, 1)  # Avoid division by zero
    return exp(-distance / (2 * kernel.α^2)) * proportion
end

struct CustomDispersal{K<:KernelFormulation}
    kernel::K
    distancemethod::Function
end

function applyrule!(data::AbstractSimData, rule::CustomDispersal, state, I::Tuple{Int, Int})
    # Assuming `state` refers to a structured data form containing multiple grids
    abundance = state.a[I...]
    carrying_capacity = state.b[I...]

    # Determine how neighbors will affect the current cell
    # You need to define `neighbors` and make sure it is appropriately initialized
    total_dispersal = 0.0
    for offset in neighbors  # `neighbors` should be a predefined list of neighbor offsets
        J = (I[1] + offset[1], I[2] + offset[2])
        if checkbounds(Bool, state.a, J)
            distance = rule.distancemethod(I, J)
            dispersal_effect = rule.kernel(distance)
            # Adjust how dispersal effect is applied based on abundance and carrying capacity
            total_dispersal += dispersal_effect * (abundance / max(carrying_capacity, 1))
        end
    end

    # Update the state with the calculated total dispersal
    data.a[I...] += total_dispersal
end


custom_dispersal = CustomDispersal(
    kernel = CustomKernel(1.0),
    distancemethod = AreaToArea(10)
)