######## IMPROVING LogisticGrowth and trophic_optimized #########
using StaticArrays, Serialization, DynamicGrids, Dispersal, Rasters, Dates, WGLMakie

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
Base.zero(x::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(fill(oneunit(T), 256), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar * 256)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar * 256)

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

# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs256
function Base.maximum(a::Matrix{MyStructs256{AbstractFloat}})
    # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

########### CREATING DATA ##########
######### raster_with_abundances ########
# Define the size of the raster
size_x, size_y = 125, 76
# Define the x and y dimensions directly
x_dim = X(1:size_x)
y_dim = Y(1:size_y)
# A 125x76 grid of random MyStructs256
raster_matrix = reshape([MyStructs256(SVector{256, Float64}(fill(10.0, 256))) for _ in 1:(size_x * size_y)], size_x, size_y)
raster_with_abundances = Raster(raster_matrix, dims=(x_dim, y_dim))
######### raster_sum ########
# Create a 125x76 grid of Bool values with first and last rows as false, others true
sum_matrix = [i == 1 || i == size_x ? false : true for i in 1:size_x, j in 1:size_y]
# Create the Bool raster with the specified dimensions
raster_sum = Raster(sum_matrix, dims=(x_dim, y_dim))

######### HERE STARTS THE RELEVANT PART #############
struct CustomKernel <: KernelFormulation
    α::AbstractFloat
end

abstract type AbstractKernelNeighborhood end

struct CustomDispersalKernel{N<:DynamicGrids.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DynamicGrids.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

outdisp = OutwardsDispersal{:raster_with_abundances, :raster_with_abundances}(;
    formulation=CustomKernel(1.0),
    distancemethod=AreaToArea(30),
    radius=1
);

neighbors_rule = Neighbors{:raster_with_abundances, :raster_with_abundances}(Moore(1)) do data, neighborhood, cell, I
    # Get the distances for this neighborhood
    dist = distances(neighborhood)
    
    # Apply the custom dispersal kernel based on the distance
    # (This is equivalent to the kernel formulation)
    kernel = CustomKernel(1.0)  # Example with α = 1.0
    weighted_sum = zero(cell.a)  # Initialize the weighted sum as an SVector
    
    for (neighbor, d) in zip(neighbors(neighborhood), dist)
        weighted_sum += kernel(d) * neighbor.a  # Apply the kernel to each neighbor's `a` field
    end
    
    # Return the weighted sum for this cell, encapsulated in a new MyStructs256
    return MyStructs256(weighted_sum)
end


init = (; 
    raster_with_abundances = raster_with_abundances
)

array_output = ArrayOutput(init, tspan = 1:10;
    fps = 10, ruleset = Ruleset(neighbors_rule))
@time s = sim!(array_output, Ruleset(outdisp))
@time r = sim!(array_output, Ruleset(neighbors_rule))

s[end].raster_with_abundances ≈ r[end].raster_with_abundances
makie_output = MakieOutput(init, tspan = 1:10;
fps = 10, ruleset = Ruleset(neighbors_rule; boundary = Reflect()),
    mask = raster_sum) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:raster_with_abundances]
    titles = ["raster_with_abundances"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:1, j in 1:1]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :raster_with_abundances
            Makie.heatmap!(ax, frame[:raster_with_abundances]; interpolate=false)
        end
    end
end

using Stencils
struct OutwardsDispersalRemix{R,W,S<:Stencils.AbstractKernelStencil, M} <: SetNeighborhoodRule{R,W}
    stencil::S
    maskbehavior::M
end

# Constructors for OutwardsDispersalRemix
function OutwardsDispersalRemix{R,W}(stencil::S; maskbehavior::Union{Dispersal.CheckMaskEdges, Dispersal.IgnoreMaskEdges}=IgnoreMaskEdges()) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersalRemix{R,W,S,typeof(maskbehavior)}(stencil, maskbehavior)
end

function OutwardsDispersalRemix{R,W}(; maskbehavior::Union{Dispersal.CheckMaskEdges, Dispersal.IgnoreMaskEdges}=IgnoreMaskEdges(), kw...) where {R,W}
    stencil = DispersalKernel(; kw...)
    OutwardsDispersalRemix{R,W,typeof(stencil),typeof(maskbehavior)}(stencil, maskbehavior)
end

@inline function applyrule!(data, rule::OutwardsDispersalRemix{R,W}, N, I) where {R,W}
    # Check if body_mass_vector is defined
    if !isdefined(Main, :body_mass_vector)
        error("Error: `body_mass_vector` is not defined in the environment. Please define it before running the simulation.")
    end
    
    # Retrieve the body_mass_vector from the environment
    body_mass_vector = Main.body_mass_vector

    N == zero(N) && return nothing

    # Check if the current cell is masked, skip if it is
    mask_data = if rule.maskbehavior === IgnoreMaskEdges() nothing else DynamicGrids.mask(data) end
    if !isnothing(mask_data) && !mask_data[I...]
        return nothing
    end

    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        (target_mod, inbounds) = DynamicGrids.inbounds(data, target)

        if inbounds && (isnothing(mask_data) || mask_data[target_mod...])
            # Adjust kernel by the body_mass_vector
            adjusted_kernel = k .* body_mass_vector

            @inbounds propagules = N .* adjusted_kernel  # Apply the adjusted kernel
            @inbounds add!(data[W], propagules, target_mod...)  # Add to neighboring cell
            sum += propagules
        end
    end
    
    # Subtract the sum of dispersal from the current cell
    @inbounds sub!(data[W], sum, I...)
    
    return nothing
end

