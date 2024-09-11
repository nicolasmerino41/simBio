include("One-click code.jl")

include("human_footprint.jl")

#TODO At prevalence 0.277 or higher you get instability
pepe = (
    birmmals = Matrix(DA_birmmals_with_abundances),
    herps = Matrix(DA_herps_with_abundances),
    k_DA = Matrix(k_DA_hf_additive),
    birmmals_richness = Matrix(DA_richness_birmmals),
    herps_richness = Matrix(DA_richness_herps)
)
DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)
pepe_state = (
    state = Matrix(DA_with_abundances),
    k_DA = Matrix(k_DA_hf_additive),
    npp_DA = Matrix(npp_DA),
    state_richness = Matrix(DA_richness)
)

caca = deepcopy(iberian_interact_NA)
self_regulation = 1.0
sigma = 0.00001
epsilon = 1.0
full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon, self_regulation))
remove_variable(:alpha)
alpha = 0.1
exp(-(1^2) / (2*(alpha^2)))
m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

function trophic_optimized(abundances, full_IM)
    # Calculate the weighted interaction directly
    return MyStructs256(SVector{256, Float64}((full_IM * abundances.a) .* abundances.a))
end
function int_Gr_for_biotic_k(state::MyStructs256, self_regulation::AbstractFloat, k_DA::MyStructs256)
    return MyStructs256(SVector{256, Float64}(self_regulation .* (k_DA.a .+ 0.00001) .*
        state.a .* (1.0 .- (state.a ./ (k_DA.a .* binary_vector .+ k_DA.a .* opposite_binary_vector)))))
end
function GLV(state::MyStructs256, k_DA::MyStructs256)
    return MyStructs256(
        SVector{256, Float64}(
            state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a)) 
        )
    )
end

biotic_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return GLV(state, k_DA)
end

biotic_rule_k = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    if any(isinf, state.a) || any(isnan, state.a)
        @error "state has NA values"
        println(I)
    end
    return MyStructs256(
        SVector{256, Float64}(
            max.(
                0.0,
                (state +
                int_Gr_for_biotic_k(state, self_regulation, k_DA)  +
                trophic_optimized(state, full_IM)).a
            )
        )
    )
end
neighbors_rule = Neighbors{:state, :state}(Moore(1)) do data, neighborhood, cell, I
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

scaling_vector = fill(1.0, 49)
scaling_vector = append!(scaling_vector, fill(1.0, 207))

neighbors_rule_size_based = Neighbors{:state, :state}(Moore(1)) do data, neighborhood, cell, I
    # Get the distances for this neighborhood
    dist = distances(neighborhood)
    
    # Initialize an SVector to store the weighted sum
    weighted_sum = zero(cell.a)  # Efficient initialization

    # Custom dispersal kernel function (as a function of distance)
    kernel = CustomKernel(1.0)

    # Iterate over neighbors and apply dispersal based on the kernel and body mass
    for (neighbor_offset, d) in zip(neighbors(neighborhood), dist)
        # Calculate kernel value based on distance
        kernel_value = kernel(d)
        
        # Adjust kernel values for each species based on body mass
        adjusted_kernel_values = kernel_value .* body_mass_vector  # Adjust the kernel by the body mass vector
        
        # Calculate the destination index for the neighbor
        dest, is_inbounds = inbounds(data, I .+ neighbor_offset)
        
        # Apply the inbounds check
        if is_inbounds
            # Accumulate weighted sum for the current neighbor with broadcasting (element-wise operations)
            weighted_sum .+= adjusted_kernel_values .* neighbor.a .- adjusted_kernel_values .* cell.a
        end
    end
    
    # Return the weighted sum encapsulated in a new MyStructs256
    return MyStructs256(weighted_sum)
end

# Example usage of OutwardsDispersalRemix in a simulation
remix_outdisp = OutwardsDispersalRemix{:state, :state}(
    formulation=CustomKernel(alpha),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(alpha),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

indisp = InwardsDispersal{:state, :state}(;
    formulation=ExponentialKernel(),
    distancemethod=AreaToArea(30)
);

##### MAKIE STATE #####
# array_output = ResultOutput(
#     pepe_state; tspan = 1:3,
#     mask = Matrix(DA_sum)
# )
# array_output2 = ResultOutput(
#     pepe_state; tspan = 1:3,
#     mask = Matrix(DA_sum)
# )
# @time a = sim!(array_output, Ruleset(biotic_rule_k))[end].state
# @time b = sim!(array_output2, Ruleset(biotic_GLV))[end].state
# a ≈ b
# MK.image(Matrix(DA_with_abundances), lambda_DA.multiplicative; colormap = :thermal, colorrange = (0, total_max))
# map_plot(Matrix(r[end].state); lambda_DA = lambda_DA.multiplicative, type = "image", palette = :thermal, colorrange = (0, total_max))
# map_plot(Matrix(DA_with_abundances); lambda_DA = lambda_DA.multiplicative, type = "image", palette = :thermal, colorrange = (0, total_max))

MakieOutput(pepe_state, tspan = 1:100;
    fps = 10, ruleset = Ruleset(biotic_GLV, remix_outdisp; boundary = Reflect()),
    mask = Matrix(DA_sum)) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:biomass, :simulated_richness, :npp, :real_richness]
    titles = ["Biomass", "Simulated Richness", "NPP", "Real Richness"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state], lambda_DA.multiplicative; colormap=custom_palette, colorrange = (0, 256))
        elseif key == :npp
            Makie.heatmap!(ax, frame[:npp_DA]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
        elseif key == :real_richness
            Makie.heatmap!(ax, frame[:state_richness]; interpolate=false, colormap=custom_palette, colorrange = (0, 256))
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.titlegap[] = 5  # Adjust the title gap to make it smaller
        ax.titlesize[] = 12 # Set the title font size
        ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
        ax.yreversed[] = true
    end
end

##### LAX NICHE #####
array_output = ArrayOutput(
    pepe; tspan = 1:100,
    mask = Matrix(DA_sum)
)
@time p = sim!(array_output, Ruleset(biotic_rule_k_herps, biotic_rule_k_birmmals, outdisp_birmmals, outdisp_herps; proc = ThreadedCPU(), boundary = Reflect()))

# include("../Non-daily-code/MetricsExploration.jl")
# include("../Daily JULIA code/New_metrics.jl")
# zeros_ones = count_zeros_ones(DA_with_presences, idx)
##### MAKIE BIRMMALS/HERPS #####
makie_output = MakieOutput(pepe, tspan = 1:1000;
    fps = 10, ruleset = Ruleset(biotic_rule_k_herps, biotic_rule_k_birmmals, outdisp_birmmals, outdisp_herps;
    boundary = Reflect()),
    mask = Matrix(DA_sum)) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:birmmals_biomass, :herps_biomass, :birmmals_richness, :herps_richness, :real_birmmals_richness, :real_herps_richness]
    titles = ["Birds+Mammals Biomass", "Herps Biomass", "Simulated Birds+Mammals Richness", "Simulated Herps Richness", "Real Birmmals Richness", "Real Herps Richness"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*3 + j]) for i in 1:2, j in 1:3]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :birmmals_biomass
            Makie.heatmap!(ax, frame[:birmmals]; interpolate=false, colormap=custom_palette)
        elseif key == :herps_biomass
            Makie.heatmap!(ax, frame[:herps]; interpolate=false, colormap=custom_palette)
        elseif key == :birmmals_richness
            Makie.image!(ax, frame[:birmmals]; colormap=custom_palette)
        elseif key == :herps_richness
            Makie.image!(ax, frame[:herps]; colormap=custom_palette)
        elseif key == :real_birmmals_richness
            Makie.heatmap!(ax, frame[:birmmals_richness]; interpolate=false, colormap=custom_palette)
        elseif key == :real_herps_richness
            Makie.heatmap!(ax, frame[:herps_richness]; interpolate=false, colormap=custom_palette)
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.titlegap[] = 5  # Adjust the title gap to make it smaller
        ax.titlesize[] = 12 # Set the title font size
        ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
        ax.yreversed[] = true
    end
end

######## TRYING TIME SERIES AND MWE FOR RAF #########
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

######### HERE STARTS THE RELEVANT PART #############
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

array_output = ResultOutput(
    init_for_timeseries; tspan = tspan,
    aux = aux,
    mask = raster_sum
)
@time s = sim!(array_output, Ruleset(timeseries_rule))
