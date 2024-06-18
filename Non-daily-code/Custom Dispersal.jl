using Pkg
# Desktop PC
Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
# Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

# meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
# meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
struct CustomKernel <: KernelFormulation
    α::Float64
end

# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{Float64}}, kernel::SVector{9, Float64})
    result_a = hood[5].a
    new_vector = MArray(SVector{256, Float64}(fill(0.0, 256)))
    # println(hood[5].a)
    for i in 1:256
        if result_a[i] > 0.0
            new_vector[i] = 1.0
        end
    end
    for (i, k) in enumerate(kernel)
        if i != 5.0
            result_a += (hood[i].a .* k) 
            result_a = result_a .* new_vector
            # println("result_a: ", sum(result_a), " hood: ", hood[i].b, " k: ", k)
            # println(typeof(archetypes))
        end
    end
      
    return MyStructs256(result_a)
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

#### Here we implement group_specific dispersal
taxonomy = CSV.File(joinpath(meta_path, "Metaweb_data", "species_codes_and_taxonomy.csv")) |> DataFrame
taxonomy = taxonomy[:,4:7]

taxonomy = taxonomy[in.(taxonomy.Species, Ref(names(species))), :]
taxonomy = taxonomy[:, [:Class, :Species]]

# Reorder the taxonomy DataFrame based on the order of species names
taxonomy = taxonomy[sortperm([findfirst(x -> x == tax, names(species)) for tax in taxonomy.Species]), :]
taxonomy[!, :archetype] .= float(0.0)
for i in 1:256
    if taxonomy[i, :Class] == "Aves" 
        taxonomy[i, :archetype] = 0.0
    else
        taxonomy[i, :archetype] = 0.0
    end
end
archetypes = SVector{256, Float64}(taxonomy.archetype) 

init = fill(0.0, (9,9)) 
init[5,5] = 100.0

output = ArrayOutput(init; tspan=1:100)
ruleset = Ruleset(indisp)
save = sim!(output, ruleset)
o = save[3 ]

init = fill(MyStructs256(SVector{256, Float64}(fill(100.0, 256))), (9,9)) 
init[2,2] = MyStructs256(SVector{256, Float64}(fill(100.0, 256)))

output = ArrayOutput(init; tspan=1:10)

save = sim!(output, indisp)
o = save[2]

struct DispersalKernel{R,N,L,H<:Neighborhood{R,N,L},K,F,C,D<:DistanceMethod} <: AbstractKernelNeighborhood{R,N,L,H}
    neighborhood::H
    kernel::K
    formulation::F
    cellsize::C
    distancemethod::D

    function DispersalKernel(
        hood::H, kernel, formulation::F, cellsize::C, distancemethod::D
    ) where {H<:Neighborhood{R,N,L},F,C,D<:DistanceMethod} where {R,N,L}
        if hood isa DispersalKernel
            hood
        else
            # Build the kernel matrix
            newkernel = scale(buildkernel(hood, formulation, distancemethod, cellsize))
            new{R,N,L,H,typeof(newkernel),F,C,D}(
                hood, newkernel, formulation, cellsize, distancemethod
            )
        end
    end

    function DispersalKernel{R,N,L,H,K,F,C,D}(
        hood::H, kernel::K, formulation::F, cellsize::C, distancemethod::D
    ) where {R,N,L,H,K,F,C,D}
        new{R,N,L,H,K,F,C,D}(hood, kernel, formulation, cellsize, distancemethod)
    end
end

function DispersalKernel(;
    radius=1,
    neighborhood=Moore(radius),
    formulation=CustomKernel(1.0),  # Use CustomKernel by default
    cellsize=1.0,
    distancemethod=CentroidToCentroid(),
)
    DispersalKernel(neighborhood, nothing, formulation, cellsize, distancemethod)
end

DispersalKernel{R}(; radius=R, kw...) where R = DispersalKernel(; radius=radius, kw...)

constructorof(::Type{<:DispersalKernel}) = DispersalKernel

function setwindow(n::DispersalKernel{R,N,L,<:Any,K,F,C,D}, buffer) where {R,N,L,K,F,C,D}
    newhood = DG.setwindow(neighborhood(n), buffer)
    DispersalKernel{R,N,L,typeof(newhood),K,F,C,D}(
        newhood, kernel(n), formulation(n), cellsize(n), distancemethod(n)
    )
end

cellsize(hood::DispersalKernel) = hood.cellsize
distancemethod(hood::DispersalKernel) = hood.distancemethod
formulation(hood::DispersalKernel) = hood.formulation

# Create an instance of the custom DispersalKernel
custom_kernel = CustomKernel(0.5)  # Define a custom alpha value
dispersal_kernel = DispersalKernel(neighborhood=Moore(1), formulation=custom_kernel, cellsize=1.0, distancemethod=CentroidToCentroid())
dispersal_kernel.neighborhood
# Use this kernel in your dispersal model
indisp = InwardsDispersal{}(;
    formulation=dispersal_kernel.formulation,
    distancemethod=dispersal_kernel.distancemethod,
    neighborhood=Moore(2)
)

ruleset = Ruleset(indisp)
full_IM = results[1][0.001]

