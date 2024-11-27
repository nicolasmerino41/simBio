######## IMPROVING LogisticGrowth and trophic_optimized #########
using StaticArrays, Serialization, DynamicGrids, Dispersal, Rasters, Dates

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

growth_layer = deepcopy(k_DA_hf_multiplicative)
for i in axes(growth_layer, 1), j in axes(growth_layer, 2)
    growth_layer[i, j] = MyStructs256(SVector{256, Float64}(growth_layer[i, j].a .* self_regulation))
end
map_plot(Matrix(growth_layer););
######### HERE STARTS THE RELEVANT PART #############
aux = (; carrycap = k_DA_hf_multiplicative, rate = growth_layer)
tspan = Date(2023, 1):Month(1):Date(2023, 8)

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

log_growth = let k_raster_aux=Aux{:k_raster}()
    Cell{:state, :state}() do data, state, I
        k_raster = get(data, k_raster_aux, I) # Important bit here!!
        return MyStructs256(SVector{256, Float64}(state.a .* k_raster))
    end
end

rate = rand([0.5, 1], (125, 76))
rate = deepcopy(utmraster)
carrycap = deepcopy(utmraster)
state = deepcopy(utmraster)

raster_with_abundances

for i in 1:125, j in 1:76
    rate[i, j] = rand([0.5, 1.0])
    carrycap[i, j] = rand()
    state[i, j] = rand([10.0, 20.0, 30.0])
end
MK.plot(state);
pepe = (;
    state = raster_with_abundances
)
##############################################
struct LogisticGrowthForMyStructs256{R,W,GR,CC,TS,S} <: Dispersal.GrowthRule{R,W}
    rate::GR
    carrycap::CC
    timestep::TS
    nsteps::S
end
function LogisticGrowthForMyStructs256{R,W}(;
    rate=INTRINSICRATE_PARAM,
    carrycap=CARRYCAP_PARAM,
    timestep,
    nsteps_type=Float64,
) where {R,W}
LogisticGrowthForMyStructs256{R,W}(rate, carrycap, timestep, zero(nsteps_type))
end

modifyrule(rule::LogisticGrowthForMyStructs256, data) = precalc_nsteps(rule, data)

@inline function applyrule(data, rule::LogisticGrowthForMyStructs256, I)
    # Get the population at index I
    N = get(data, :state, I)

    N.b > zero(N.b) || return zero(N)

    # Get rate and carrycap from data using the index I
    rt = get(data, rule.rate, I...) * rule.nsteps
    k = get(data, rule.carrycap, I...)

    new_a = if rt > zero(rt)
        # Apply logistic growth formula to each element in the SVector
        (N.a .* k) ./ (N.a .+ (k - N.a) .* exp(-rt))
    else
        N.a .* exp(rt)
    end

    # Create a new MyStructs256 object with updated .a and recalculated sum .b
    new_N = MyStructs256(SVector{256, Float64}(new_a), sum(new_a))

    # Return the result, bounded between zero and carrying capacity
    return MyStructs256(min.(max.(zero(new_N.a), new_N.a), k), sum(min.(max.(zero(new_N.a), new_N.a), k)))
end


@inline function applyrule(data, rule::LogisticGrowthForMyStructs256, Ns::AbstractArray, I)
    # Get the state from the grid at index I
    N = get(data, :state, I)
    
    # Get growth rates and carrying capacities from data at index I
    rts = get(data, rule.rate, I...) .* rule.nsteps
    ks = get(data, rule.carrycap, I...)

    # Apply logistic growth to each element in Ns
    new_Ns = map(Ns, rts, ks) do N::MyStructs256, rt, k
        new_a = if rt > zero(rt)
            # Apply logistic growth element-wise to SVector in N.a
            (N.a .* k) ./ (N.a .+ (k .- N.a) .* exp(-rt))
        else
            N.a .* exp(rt)
        end
        # Create new MyStructs256 with updated .a and recalculated sum .b
        MyStructs256(SVector{256, Float64}(new_a), sum(new_a))
    end

    # Return the result, bounded between 0 and carrying capacity
    return map(new_Ns, ks) do new_N, k
        # Apply bounds (0 <= N <= k) element-wise
        bounded_a = min.(max.(zero(new_N.a), new_N.a), k)
        MyStructs256(SVector{256, Float64}(bounded_a), sum(bounded_a))
    end
end

@inline function applyrule(data, rule::LogisticGrowthForMyStructs256, I::CartesianIndex)
    # Get the MyStructs256 object at index I from the state grid
    N = get(data, :state, I)

    # Ensure we have a valid population (skip if zero)
    N.b > zero(N.b) || return zero(N)

    # Get the rate and carrying capacity at index I
    rt = get(data, rule.rate, I...) * rule.nsteps
    k = get(data, rule.carrycap, I...)

    # Apply logistic growth formula to each element in the SVector (N.a)
    new_a = if rt > zero(rt)
        (N.a .* k) ./ (N.a .+ (k - N.a) .* exp(-rt))
    else
        N.a .* exp(rt)
    end

    # Create a new MyStructs256 object with updated .a and recalculated sum .b
    new_N = MyStructs256(SVector{256, Float64}(new_a), sum(new_a))

    # Return the result, bounded between zero and carrying capacity
    return MyStructs256(min.(max.(zero(new_N.a), new_N.a), k), sum(min.(max.(zero(new_N.a), new_N.a), k)))
end


log_growth_rule = LogisticGrowthForMyStructs256{:state, :state}(; 
    rate = Aux(:rate),
    carrycap = Aux(:carrycap),
    timestep = 1,
    nsteps_type = Int64
)

array_output = ResultOutput(
    pepe; tspan = 1:100,
    aux = aux,
    mask = raster_sum
)

# Run the simulation with the updated rule
@time s = sim!(array_output, Ruleset(log_growth_rule))

makie_output = MakieOutput(pepe, tspan = 1:2000;
    fps = 10, ruleset = Ruleset(log_growth_rule),
    aux = aux, mask = raster_sum) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:state]
    titles = ["state"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:1, j in 1:1]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :state
            Makie.plot!(ax, frame[:state]; interpolate=false, colormap=custom_palette)
        end
    end
end