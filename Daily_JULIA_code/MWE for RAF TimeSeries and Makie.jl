######## TRYING TIME SERIES AND MWE FOR RAF #########
using StaticArrays, DynamicGrids, Dispersal, Rasters, Dates, Makie, WGLMakie

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

# # Define what a NaN is for MyStructs256
# Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

# # Adding a method in the sum function for MyStructs256
# function Base.sum(structs::MyStructs256...)
#     # Sum the 'a' vectors
#     summed_a = sum([s.a for s in structs])

#     # Sum the 'b' values
#     summed_b = sum([s.b for s in structs])

#     # Create a new MyStructs256 instance with the summed results
#     return MyStructs256(summed_a, summed_b)
# end

# # Adding a method to maximum
# # Define maximum for MyStructs256
# function Base.maximum(a::MyStructs256, b::MyStructs256)
#     return MyStructs256(max.(a.a, b.a))
# end

# # Define maximum for MyStructs256 with a scalar
# function Base.maximum(a::MyStructs256, b::AbstractFloat)
#     return MyStructs256(max.(a.a, b))
# end

# # Define maximum for a scalar with MyStructs256
# function Base.maximum(a::AbstractFloat, b::MyStructs256)
#     return MyStructs256(max.(a, b.a))
# end

# # Define maximum for MyStructs256
# function Base.maximum(a::MyStructs256)
#     return maximum(a.a)
# end

# # Define maximum for a matrix of MyStructs256
# function Base.maximum(a::Matrix{MyStructs256{AbstractFloat}})
#     # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
#     return maximum(map(x -> x.b, a))
# end

########### CREATING DATA ##########
######### raster_with_abundances ########
# Define the size of the raster
size_x, size_y = 125, 76
# Define the x and y dimensions directly
x_dim = X(1:size_x)
y_dim = Y(1:size_y)
# A 125x76 grid of random MyStructs256
raster_matrix = reshape([MyStructs256(SVector{256, Float64}(rand(256))) for _ in 1:(size_x * size_y)], size_x, size_y)
raster_with_abundances = Raster(raster_matrix, dims=(x_dim, y_dim))
######### raster_sum ########
# Create a 125x76 grid of Bool values with first and last rows as false, others true
sum_matrix = [i == 1 || i == size_x ? false : true for i in 1:size_x, j in 1:size_y]
# Create the Bool raster with the specified dimensions
raster_sum = Raster(sum_matrix, dims=(x_dim, y_dim))
######## combined_raster ########
time_dim = Ti([Date(2023, i, 1) for i in 1:8])  # 8 time points
# Create 8 individual rasters with random Int16 values
raster_layers = [
    Raster(rand(Int16, size_x, size_y), dims=(x_dim, y_dim)) for _ in 1:8
]
# Set the first and last rows to -32768 to represent missing values
for raster in raster_layers
    raster[1, :] .= -32768  # First row
    raster[end, :] .= -32768  # Last row
end
# Combine the rasters into a 3D raster
combined_raster = Raster(cat(raster_layers..., dims=3), dims=(x_dim, y_dim, time_dim))

####################################
aux = (; combined_raster=combined_raster)
tspan = Date(2023, 1):Month(1):Date(2023, 8)

init_for_timeseries = (;
    state = raster_with_abundances
)

self_regulation = 0.01
function int_Gr_for_timeseries(state::MyStructs256, self_regulation::AbstractFloat, combined_raster)
    return MyStructs256(SVector{256, Float64}(state.a .* self_regulation .* combined_raster))
end

timeseries_rule = let combined_raster_aux=Aux{:combined_raster}()
    Cell{:state, :state}() do data, state, I
        combined_raster = get(data, combined_raster_aux, I) # Important bit here!!
        merged_state = state + int_Gr_for_timeseries(state, self_regulation, combined_raster)
        return MyStructs256(SVector{256, Float64}(merged_state.a))
    end
end

# array_output = ResultOutput(
#     init_for_timeseries; tspan = tspan,
#     aux = aux,
#     mask = raster_sum
# )
# @time s = sim!(array_output, Ruleset(timeseries_rule))

###### HERE THE RELEVANT PART Makie convert_arguments #############
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyStructs256, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > 0.1, 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end
Makie.heatmap(raster_with_abundances) # This works fine
Makie.image(raster_with_abundances) # And this too

makie_output = MakieOutput(init_for_timeseries, tspan = tspan;
    fps = 10, ruleset = Ruleset(timeseries_rule),
    aux = aux, mask = raster_sum) do (; layout, frame)

    plot_keys = [:biomass, :richness]
    titles = ["Biomass", "Simulated Richness"]

    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:1, j in 1:2]
    
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            # Makie.heatmap!(ax, frame[:state]; interpolate=false)
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state]; interpolate=false)
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title
        ax.titlegap[] = 5
        ax.titlesize[] = 12
        # ax.titlecolor[] = RGBA(0, 0, 0, 1)
        ax.yreversed[] = false
    end
end