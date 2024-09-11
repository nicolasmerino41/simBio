using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
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
# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.zero(::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
# Create an instance of MyStructs256
u = MyStructs256(SVector{256, Float64}(fill(0.0, 256)), 2700.0)
# Get the zero value for MyStructs256
z = zero(u)
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))
# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)
# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*256)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*256)
# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

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
function Base.zeros(dims::NTuple{2, Int}, type = nothing)
    if type == MyStructs256{Float64}
        return [MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end
zeros((10, 10), MyStructs256{Float64})
################# MYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    α::Float64
end

abstract type AbstractKernelNeighborhood end

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
            if center.a[i] > 0.0
            dist = distance(I, hood.coords[j])
            result_a += kernel.formulation(dist) * neighbor.a[i]
            end
        end
    end
    return MyStructs256(result_a)
end
# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{Float64}}, kernel::SVector{9, Float64})
    
    result_a = fill(0.0, 256)
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    return MyStructs256(result_a)
end
function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyStructs256{Float64}}, kernel::SVector{25, Float64})
    
    result_a = fill(0.0, 256)
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    # println(sum(result_a))
    # result_b = sum(result_a)
    return MyStructs256(result_a) #, result_b)
end
#################### FUNCTIONS ############################
################################################################
################################################################
function growth(abundance::AbstractFloat, self_regulation::AbstractFloat, K::AbstractFloat)
    # Assuming growth is a simple function of abundance, self-regulation, and carrying capacity (K)
    return self_regulation * K * abundance * (1 - abundance / K)
end
function intrinsic_growth_256(abundance::MyStructs256, self_regulation::AbstractFloat, K::Union{AbstractFloat, AbstractArray})
    return MyStructs256((self_regulation * K) .* abundance.a .* (1.0 .- (abundance.a ./ K)))
end
function trophic(abundances, A_matrix)
    return sum(abundances.a * abundances.a' .* A_matrix, dims=2)
end
function trophic_optimized(abundances, A_matrix)
    # Calculate the weighted interaction directly
    interaction = A_matrix * abundances.a
    return MyStructs256(SVector(interaction .* abundances.a))
end
function competition(abundances, A_matrix)
    interaction = A_matrix * abundances.a
    return MyStructs256(SVector(interaction .* abundances.a))
end
function merge_intr_troph(intr, troph)
    return MyStructs256(SVector(intr.a .+ troph.a))
end
########################## RULES #################################
##################################################################
int_gr = Cell{}() do data, state, I
    return intrinsic_growth_256(state, 0.01, 1000.0)
end
troph = Cell{}() do data, state, I
    return trophic_optimized(state, full_IM)
end
# Merge rule with non-negative constraint for MyStructs256
merge_rule = Cell{}() do data, state, I
    if -0.1 < state.b
        
        merged_state = state + MyStructs256(SVector{256}(
            intrinsic_growth_256(state, 0.01, 100.0).a .+ #= You can switch the 100-0 by
            transposed_npp[I...] if needed =#
            trophic_optimized(state, full_IM).a)
        )
    # if any(merged_state.a .< 0.0) || any(isnan.(merged_state.a))
    #     println("merged_state is NA", state, merged_state)
    # end
    return MyStructs256(max.(merged_state.a, 0.0))
    end
end
##################### DISPERSAL RULES ##############################
####################################################################
dispersal_rule = Neighbors{}() do data, neighborhood, cell, I
    return MyStructs256((cell + sum(neighborhood)*0.1).a - cell.a*0.1*length(neighborhood))
end
indisp = InwardsDispersal{}(;
    formulation=CustomKernel(0.0002),
    distancemethod=AreaToArea(30)
)
ruleset = Ruleset(merge_rule, indisp)

############ DESERIALIZING DATA ############################
############################################################
A_matrix_list = deserialize(joinpath(meta_path, "A_matrix_list.jls"))::Vector{Matrix{Float64}}
full_A_matrix = A_matrix_list[3]
full_IM_list = deserialize(joinpath(meta_path, "fullIM_list.jls"))::Vector{Vector{Matrix{Float64}}}
full_IM = full_IM_list[2][10]

######### ARRAY OUTPUT ################ 
#######################################
init = (grid1 =  [MyStructs256(SVector{256, Float64}([rand([0.0, rand(10:100)]) for _ in 1:256])) for _ in 1:77, _ in 1:77])
init[80:200] = [MyStructs256(SVector{256, Float64}(fill(NaN, 256))) for _ in 1:121]
ext= Extent(xmin=0, xmax=76, ymin=0, ymax=76)
# masklayer = boolmask(init, to = )
output = ArrayOutput(init; tspan=1:100, mask = boolmask(init));
@time sim!(output, ruleset)

########### VISUAL OUTPUTS ################
###########################################
# Visualization settings
DynamicGrids.to_rgb(scheme::ObjectScheme, obj::MyStructs256) = ARGB32(
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0),
    clamp(obj.b/25600, 0.0, 1.0)
)
DynamicGrids.to_rgb(scheme, obj::MyStructs256) = get(scheme, clamp(obj.b, 0.0, 1.0))

##### GifOutput #####
ruleset = Chain(indisp)
# chained_ruleset = Chain(indisp, merge_rule)
permutedims_D = permutedims(D, (2, 1)) # If you put in init it works fine
gif_output = GifOutput(
    permutedims_D; tspan = 1:50,
    filename = "mystruct256_withgoodisp.gif",
    ruleset = ruleset,
    minval=0, maxval=5000,
    fps =10,
    mask = boolmask(transposed_npp), 
    scheme=ColorSchemes.inferno,
    zerocolor=RGB24(0.0)
);
# Run the simulation
result = sim!(gif_output, ruleset);

############# MyStructs256Raster DimensionalArray and MakieOutput ################
################ It works with GifOutput but not MakieOutput #####################
##################################################################################
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"
bioclim_paths = RasterDataSources.getraster(WorldClim{BioClim}, (5,7,8,12))
# bioclim_stack = RasterStack(WorldClim{BioClim}, (5, 7, 8, 12), res="10m")
bioclim_stack = RasterStack(bioclim_paths)
bioclim_5 = bioclim_stack[:bio5]
spain = bioclim_5[X(-10 .. 4), Y(36 .. 45)]
spain = replace_missing(spain, 0)

# Create a DimArray with 84 rows and 54 columns filled with MyStructs256 objects
D = DimArray(reshape([MyStructs256(SVector{256, Float64}(rand([0.0, rand(10:100)], 256))) for _ in 1:84*54], 84, 54), (Dim{:a}(1:84), Dim{:b}(1:54)))
typeof(D)
for row in axes(spain, 1), col in axes(spain, 2)
    if spain[row, col] == 0 #&& (row in 90:100) && (col in 90:100)
        D[row, col] = MyStructs256(SVector{256, Float64}(rand(Float64, 256)).*0.0)
    end
end
reversed_D = reverse(D, dims=2)
MyStructs256Plot(reversed_D);
 
# Create a MakieOutput object
makie_output = MakieOutput(reversed_D, tspan = 1:100; fps = 10, ruleset = Ruleset(merge_rule)) do (; layout, frame)
    ax = Axis(layout[1, 1])
    hideydecorations!(ax)
    hidexdecorations!(ax)
    # Plot the frame data
    image!(ax, frame; interpolate=false, colormap=:inferno)
end

