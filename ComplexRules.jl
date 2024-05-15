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
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{256, T}
    b::T
    
    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{256, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{256, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# # Define zero and oneunit for MyStructs256
# Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
# Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))

# # Comparison based on 'b' field
# Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)

# # Element-wise arithmetic operations ensuring 'b' is recalculated correctly
# Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
# Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
# Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
# Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
# Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*256)
# Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*256)

# Define the zero method for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(ntuple(_ -> zero(T), 256)), zero(T))

# Define the oneunit method for MyStructs256
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(ntuple(_ -> oneunit(T), 256)), oneunit(T))

# Define isless method for MyStructs256
Base.isless(a::MyStructs256, b::MyStructs256) = isless(a.b, b.b)

# Ensure arithmetic operations recalculate 'b'
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar)
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar)
Base.:+(x1::MyStructs256, x2::MyStructs256) = MyStructs256(x1.a .+ x2.a)
Base.:-(x1::MyStructs256, x2::MyStructs256) = MyStructs256(x1.a .- x2.a)

function growth(abundance::AbstractFloat, self_regulation::AbstractFloat, K::AbstractFloat)
    # Assuming growth is a simple function of abundance, self-regulation, and carrying capacity (K)
    return self_regulation * K * abundance * (1 - abundance / K)
end

# Adding a method in the sum function for MyStructs256
function Base.sum(structs::MyStructs256...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs256 instance with the summed results
    return MyStructs256(summed_a, summed_b)
end
# Example instances
init1 = MyStructs256(SVector{256, Float64}(fill(1.0, 256)))
init2 = MyStructs256(SVector{256, Float64}(fill(2.0, 256)))
# Sum the instances
result = sum(init1, init2, init1)
init1*0.1

A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}
full_A_matrix = A_matrix_list[3]

full_IM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
full_IM = full_IM_list[3][10]

function intrinsic_growth_256(abundance::MyStructs256, self_regulation::AbstractFloat, K::AbstractFloat)
    return MyStructs256(self_regulation .* K .* abundance.a .* (1.0 .- (abundance.a ./ K)))
end

function trophic(abundances, A_matrix)
    return sum(abundances.a * abundances.a' .* A_matrix, dims=2)
end
function trophic_optimized(abundances, A_matrix)
    # Calculate the weighted interaction directly
    interaction = A_matrix * abundances.a
    return MyStructs256(SVector(interaction .* abundances.a))
end
function merge_intr_troph(intr, troph)
    return MyStructs256(max.(0.0, SVector(intr.a .+ troph.a)))
end
############ RULES ####################
int_gr = Cell{}() do data, state, I
    return intrinsic_growth_256(state, 0.01, 1000.0)
end
troph = Cell{}() do data, state, I
    return trophic_optimized(state, full_IM)
end
merge_rule = Cell{}() do data, state, I
    return merge_intr_troph(intrinsic_growth_256(state, 0.01, 1000.0), trophic_optimized(state, full_IM))
end
##################### DISPERSAL ##############################
dispersal_rule = Neighbors{}() do data, neighborhood, cell, I
    return MyStructs256((cell + sum(neighborhood)*0.1).a - cell.a*0.1*length(neighborhood))
end
indisp = InwardsDispersal{}(;
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
ruleset = Ruleset(merge_rule, indisp, proc = SingleCPU())
######### ARRAY OUTPUT ################ Everything works fine until here
init = (grid1 = [MyStructs256(SVector{256}(rand(Float64, 256).*100)) for _ in 1:77, _ in 1:77])
output = ArrayOutput(init; tspan=1:100)
@time sim!(output, ruleset)




# Custom Kernel for MyStructs256
struct MyStructs256KernelFormulation <: KernelFormulation
    λ::Float64  # Example parameter for the kernel
end

# Define the functor method for MyStructs256KernelFormulation
(k::MyStructs256KernelFormulation)(distance) = exp(-distance / k.λ)

#= Here we define a custom dispersal kernel that applies the kernel formulation 
to each element of the SVector in MyStructs256 =#
abstract type AbstractKernelNeighborhood end
struct MyStructs256DispersalKernel{N<:DG.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function MyStructs256DispersalKernel(; 
    neighborhood::DG.Neighborhood=Moore(1), 
    formulation::KernelFormulation=MyStructs256KernelFormulation(1.0)
)
    MyStructs256DispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

# Define radius method for MyStructs256DispersalKernel
DynamicGrids.radius(kernel::MyStructs256DispersalKernel) = radius(kernel.neighborhood)

# Define unsafe_neighbors method for MyStructs256DispersalKernel
function DynamicGrids.unsafe_neighbors(kernel::MyStructs256DispersalKernel, data::DynamicGrids.GridData, I::CartesianIndex)
    hood = DynamicGrids.neighbors(kernel.neighborhood, data, I)
    neighbors(kernel, hood, data[I], I)
end

# Extend neighbors function for dispersal
function DynamicGrids.neighbors(kernel::MyStructs256DispersalKernel, hood, center, I)
    result = zero(center)
    for i in 1:256
        for (j, neighbor) in enumerate(hood)
            dist = distance(I, hood.coords[j])
            result.a[i] += kernel.formulation(dist) * neighbor.a[i]
        end
    end
    result.b = sum(result.a)
    return result
end

# Define a rule using the custom kernel
rule = Neighbors{}(MyStructs256DispersalKernel()) do data, hood, (s1,), I
    sum(hood)
end

# Define to_rgb for visualizing MyStructs256
DynamicGrids.to_rgb(::ObjectScheme, obj::MyStructs256) = begin
    intensity = clamp(obj.b / 25600, 0.0, 1.0)
    RGB(intensity, 0, 0)
end

# Create initial state with MyStructs256 objects
init = [MyStructs256(SVector{256, Float32}(rand(10.0:100.0, 256))) for _ in 1:1000, _ in 1:1000]

# Setup the output
output = ArrayOutput(init; tspan=1:3)

# Run the simulation
sim!(output, rule)

# Visualization
ruleset = Ruleset(rule)

output_gif = GifOutput(
    init; tspan = 1:100,
    filename = "test.gif",
    ruleset = ruleset,
    scheme = ObjectScheme()
)

# Run the simulation to generate GIF
DynamicGrids.sim!(output_gif, ruleset)







###################### Custom Kernel for MyStructs256 #######################
#############################################################################
#############################################################################
# Kernel formulation
struct MyStructs256KernelFormulation <: Dispersal.KernelFormulation
    λ::Float64
end
(k::MyStructs256KernelFormulation)(distance) = exp(-distance / k.λ)
abstract type AbstractKernelNeighborhood end
# Custom dispersal kernel
struct MyStructs256DispersalKernel{N<:DG.Neighborhood, F<:Dispersal.KernelFormulation} <: DG.Neighborhood
    neighborhood::N
    formulation::F
end

function MyStructs256DispersalKernel(; 
    neighborhood::Neighborhood=Moore(1), 
    formulation::Dispersal.KernelFormulation=MyStructs256KernelFormulation(1.0)
)
    MyStructs256DispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

# Extend neighbors function for dispersal
function DynamicGrids.neighbors(kernel::MyStructs256DispersalKernel, hood, center, I)
    result = zero(center)
    for i in 1:256
        for (j, neighbor) in enumerate(hood)
            dist = Dispersal.distance(I, hood.coords[j])
            result.a[i] += kernel.formulation(dist) * neighbor.a[i]
        end
    end
    result.b = sum(result.a)
    return result
end

# Define a rule using the custom kernel
rule = Neighbors{Tuple{:grid1}, :grid1}(MyStructs256DispersalKernel()) do data, hood, (s1,), I
    sum(hood)
end


