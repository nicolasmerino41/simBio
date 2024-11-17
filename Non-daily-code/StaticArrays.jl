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
############### STRUCT ###################
struct MyStructs{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{3, T}
    b::T
    
    # Custom constructor for automatic sum calculation
    function MyStructs(a::SVector{3, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs(a::SVector{3, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs
Base.zero(::Type{MyStructs{T}}) where {T <: AbstractFloat} = MyStructs(SVector{3, T}(zero(T), zero(T), zero(T)), zero(T))
Base.oneunit(::Type{MyStructs{T}}) where {T <: AbstractFloat} = MyStructs(SVector{3, T}(oneunit(T), oneunit(T), oneunit(T)), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs, y::MyStructs) = isless(x.b, y.b)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs, y::MyStructs) = MyStructs(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs, y::MyStructs) = MyStructs(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs, scalar::Real) = MyStructs(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs, scalar::Real) = MyStructs(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs, scalar::Real) = MyStructs(x.a .- scalar, x.b - scalar*3)
Base.:+(x::MyStructs, scalar::Real) = MyStructs(x.a .+ scalar, x.b + scalar*3)


a = MyStructs(SVector(150.0, 200.0, 100.0))  # Uses the auto-sum constructor
b = MyStructs(SVector(150.0, 200.0, 100.0), 300.0)  # Uses the explicit constructor
############ RULES ####################
rule = Cell{}() do data, state, I
    MyStructs(state.a .* 1.1)  # This will double each element of 'a' and recalculate 'b'
end
rule1 = Neighbors() do data, neighborhood, cell, I
    
    return MyStructs(data[I...].a .* 2)  
end
######### ARRAY OUTPUT ################
init = (grid1 = [MyStructs(SVector(rand(10.0:100.0), rand(10.0:100.0), rand(10.0:100.0))) for _ in 1:1000, _ in 1:1000])
output = ArrayOutput(init; tspan=1:3)
sim!(output, rule)
############## GIF ####################
ruleset = Ruleset(rule)

DynamicGrids.to_rgb(scheme::ObjectScheme, obj::MyStructs) = ARGB32(clamp(obj.b , 0.0, 1.0), clamp(obj.b , 0.0, 1.0), clamp(obj.b , 0.0, 1.0))
DynamicGrids.to_rgb(scheme, obj::MyStructs) = get(scheme, obj.b)

init = (grid1 = [MyStructs(SVector(rand(10.0:100.0), rand(10.0:100.0), rand(10.0:100.0))) for _ in 1:10, _ in 1:10])
scheme = ObjectScheme()  # Assuming ObjectScheme is already tailored for use with MyStructs

output = GifOutput(init; tspan=1:100, ruleset=Ruleset(rule), filename="grey_mystruct.gif", scheme=scheme, minval=0, maxval=1000, fps = 1)
sim!(output, rule)
############# MAKIE and multiple grids ####################
cell = Cell{Tuple{:a, :b}, :a}() do data, (a,b), I
    return MyStructs(a .+ growth.(a.a, self_regulation, b))
end
indisp_sarray = InwardsDispersal(
    formulation = ExponentialKernel(λ=0.0125),
    distancemethod = AreaToArea(30)
)
pepe_for_makie = (a = reverse(inits[1], dims=2), b = reversed_npp)

init_my_struct = deepcopy(inits[1])
for row in axes(init_my_struct, 1), col in axes(init_my_struct, 2)
    if !isnan(init_my_struct[row, col]) && (row in 45:50) && (col in 45:50)
        init_my_struct[row, col] = MyStructs(SVector(rand(10.0:100.0), rand(10.0:100.0), rand(10.0:100.0)))
    else
        init_my_struct[row, col] = MyStructs(SVector(0.0,0.0, 0.0))
    end
end

output = MakieOutput(pepe_for_makie;
        tspan=1:100,
        ruleset=Ruleset(indisp_sarray),
        fps=1,
        mask=masklayer_for_makie
) do (; layout, frame)
    ax1 = Axis(layout[1, 1])
    ax2 = Axis(layout[1, 2])
    image!(ax1, frame.a)
    image!(ax2, frame.b)
end

# Assuming `inits[1]` gives the dimensions or similarly accessible raster data
dims = size(inits[1])  # Extract dimensions from an existing raster for continuity

# Create an array of MyStructs to serve as the raster
init_my_struct = Array{MyStructs{Float64}, 2}(undef, dims...)

# Populate the raster
for row in 1:dims[1], col in 1:dims[2]
    if row in 45:50 && col in 45:50
        # Assign random MyStructs values for specific cells within a specified range
        init_my_struct[row, col] = MyStructs(SVector(rand(10.0:100.0), rand(10.0:100.0), rand(10.0:100.0)))
    else
        # Assign zero MyStructs values outside that range
        init_my_struct[row, col] = MyStructs(SVector(0.0, 0.0, 0.0))
    end
end
MK.plot(init_my_struct)
using Makie

@MK.recipe(MyStructsPlot, MyStructs) do scene
    Theme(
        colormap = :viridis,
        colorrange = nothing
    )
end

parent_init = parent(inits[1])
DimensionalData.dims(parent_init, que)

dims(inits[1])

matrix_struct = Matrix{MyStructs{Float32}}(undef, 84, 54)
for row in axes(parent_init, 1), col in axes(parent_init, 2)
    if !isnan(parent_init[row, col]) && (row in 45:50) && (col in 45:50)
        matrix_struct[row, col] = MyStructs(SVector(rand(Float32(10.0):Float32(100.0)),
            rand(Float32(10.0):Float32(100.0)), rand(Float32(10.0):Float32(100.0))))
    else
        matrix_struct[row, col] = MyStructs(SVector(Float32(0.0), Float32(0.0), Float32(0.0)))
    end
end
lookup(inits[1], X)
# Create a new Raster using the data and metadata from the reference raster
custom_raster = DimArray(matrix_struct)

DimensionalData..Dimensions(inits[1])
function Makie.plot!(plot::Matrix{MyStructs{Float32}})
    dims = size(plot[:input][1])
    values = [plot[:input][1][row, col].b for row in 1:dims[1], col in 1:dims[2]]
    image!(plot, values; plot[:attributes]...)
end
##################### SAME FOR 256 species ################
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

# Define zero and oneunit for MyStructs
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(ntuple(_ -> zero(T), 256)))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(ntuple(_ -> oneunit(T), 256)))

# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*3)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*3)

init = (grid1 = [MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0)) for _ in 1:84, _ in 1:54])

rule = Cell{}() do data, state, I
    MyStructs256(state.a .* 1.1)  # This will double each element of 'a' and recalculate 'b'
end

@time output = ArrayOutput(init; tspan=1:100)
@time sim!(output, rule)

############## RECIPE FOR MAKIE
@MK.recipe(MyStructs256Plot, mystructs_matrix) do scene
    Attributes(
        colormap = :viridis,  # Choose a colormap
        colorrange = nothing
    )
end

function Makie.plot!(plot::MyStructs256)
    # Extract the matrix from the plot object
    mystructs_matrix = plot[1][]
    nrows, ncols = size(mystructs_matrix)

    # Prepare a matrix to hold the scalar 'b' values for visualization
    values = Matrix{Float64}(undef, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        values[i, j] = mystructs_matrix[i, j].b  # Assuming you want to plot the scalar 'b'
    end

    # Create a heatmap plot using the values extracted
    heatmap!(plot, values; plot.attributes...)
    return plot
end

using StaticArrays

# Example matrix of MyStructs256
matrix_of_structs = [MyStructs256(SVector{256, Float64}(rand(256)), rand()) for i in 1:10, j in 1:10]

# Plot using the custom recipe
scene = Scene();
mystructs256plot = MK.plot(matrix_of_structs)
display(scene)
