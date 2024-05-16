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

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*256)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*256)

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

# Adding a method to maximum
# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256, b::MyStructs256)
    return MyStructs256(max.(a.a, b.a))
end
# Define maximum for MyStructs256 with a scalar
function Base.maximum(a::MyStructs256, b::AbstractFloat)
    return MyStructs256(max.(a.a, b))
end

# Define maximum for a scalar with MyStructs256
function Base.maximum(a::AbstractFloat, b::MyStructs256)
    return MyStructs256(max.(a, b.a))
end


# Example instances
init1 = MyStructs256(SVector{256, Float64}(fill(1.0, 256)))
init2 = MyStructs256(SVector{256, Float64}(fill(2.0, 256)))
# Sum the instances
result = sum(init1, init2, init1)

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
    return MyStructs256(SVector(intr.a .+ troph.a))
end
############ RULES ####################
int_gr = Cell{}() do data, state, I
    return intrinsic_growth_256(state, 0.01, 1000.0)
end
troph = Cell{}() do data, state, I
    return trophic_optimized(state, full_IM)
end
# Merge rule with non-negative constraint for MyStructs256
merge_rule = Cell{}() do data, state, I
    merged_state = state + merge_intr_troph(intrinsic_growth_256(state, 0.01, 100.0), trophic_optimized(state, full_IM))
    if any(merged_state.a .< 0.0) || any(isnan.(merged_state.a))
        println("merged_state is NA", state, merged_state)
    end
    return MyStructs256(max.(merged_state.a, 0.0))  # Ensure non-negative values
end
##################### DISPERSAL ##############################
dispersal_rule = Neighbors{}() do data, neighborhood, cell, I
    return MyStructs256((cell + sum(neighborhood)*0.1).a - cell.a*0.1*length(neighborhood))
end
indisp = InwardsDispersal{}(;
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
ruleset = Ruleset(merge_rule, indisp)
######### ARRAY OUTPUT ################ Everything works fine until here
init = (grid1 = [MyStructs256(SVector{256}(fill(10.0, 256))) for _ in 1:77, _ in 1:77])
output = ArrayOutput(init; tspan=1:100)
@time sim!(output, ruleset)

output = GifOutput(init; tspan=1:100, ruleset=ruleset, 
    filename="mystruct256.gif", scheme=scheme, 
    fps = 50,
    minval=0, maxval=500
)
@time sim!(output, ruleset)

struct CustomKernel <: KernelFormulation
    α::Float64
end
(kernel::CustomKernel)(distance) = exp(-distance / (2 * kernel.α^2))

# Define MyStructs256
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{256, T}
    b::T
    
    function MyStructs256(a::SVector{256, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{256, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define kernel formulation
struct CustomKernel <: KernelFormulation
    α::Float64
end

(kernel::CustomKernel)(distance) = exp(-distance / kernel.α)

# Define custom dispersal kernel
struct CustomDispersalKernel{N<:DG.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DG.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

# Define neighbors for custom kernel
function DynamicGrids.neighbors(kernel::CustomDispersalKernel, hood, center::MyStructs256, I)
    result_a = zero(center.a)
    for i in 1:256
        for (j, neighbor) in enumerate(hood)
            dist = distance(I, hood.coords[j])
            result_a = result_a + kernel.formulation(dist) * neighbor.a[i]
        end
    end
    return MyStructs256(result_a)
end

# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{Float64}}, kernel::SVector{9, Float64})
    result_a = zero(SVector{256, Float64})
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    result_b = sum(result_a)
    return MyStructs256(result_a, result_b)
end

# Ensure necessary arithmetic operations
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar)
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar)
Base.:+(x1::MyStructs256, x2::MyStructs256) = MyStructs256(x1.a .+ x2.a)
Base.:-(x1::MyStructs256, x2::MyStructs256) = MyStructs256(x1.a .- x2.a)
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))
Base.isless(a::MyStructs256, b::MyStructs256) = isless(a.b, b.b)

# Create initial state with MyStructs256 objects
init = [MyStructs256(SVector{256, Float64}([rand([0.0, rand(10:100)]) for _ in 1:256])) for _ in 1:77, _ in 1:77]
init = [MyStructs256(SVector{256, Float64}(rand(10.0:100.0, 256))) for _ in 1:77, _ in 1:77]

# Setup the output
output = ArrayOutput(init; tspan=1:100)

# Define the rule and ruleset
indisp = InwardsDispersal{}(;
    formulation=CustomKernel(0.0125),
    distancemethod=AreaToArea(30)
)
ruleset = Ruleset(indisp)
chained_ruleset = Chain(indisp, merge_rule)
# Run the simulation
@time sim!(output, ruleset)

# Visualization
max_value = maximum(init).b
min_value = minimum(init).b
# Define the maximum value based on your data
const MAX_VALUE = 15686.0
const MIN_VALUE = 12619.0
value_random = 20
DynamicGrids.to_rgb(scheme::ObjectScheme, obj::MyStructs256) = ARGB32(
    clamp(obj.b, 0.0, 1.0),
    clamp(obj.b, 0.0, 1.0),
    clamp(obj.b, 0.0, 1.0)
)
DynamicGrids.to_rgb(scheme, obj::MyStructs256) = get(scheme, obj.b)

ruleset = Ruleset(indisp)
chained_ruleset = Chain(merge_ruleindisp)
gif_output = GifOutput(
    init; tspan = 1:50,
    filename = "mystruct256_withgoodisp.gif",
    ruleset = merge_rule,
    scheme = ObjectScheme()
)
@time sim!(gif_output, merge_rule)

# output = MakieOutput(init, tspan=1:100; ruleset=ruleset) do (; layout, frame)
#     ax = Axis(layout[1, 1])
#     image!(ax, frame; interpolate=false, colormap=:inferno)
# end

# Going back to simulating it in a Raster
gif_output_raster = GifOutput(
    transposed_init; tspan = 1:50,
    filename = "mystruct256_raster.gif",
    ruleset = merge_rule,
    scheme = ObjectScheme()
)
@time sim!(gif_output_raster, merge_rule)