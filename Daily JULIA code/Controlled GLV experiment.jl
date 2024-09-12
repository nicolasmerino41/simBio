######## PARAMETERS ########
num_species = 100
# connectance = 0.11
# sigma = 0.1
self_regulation = -1.0
carrying_capacity = 1.0
presence_threshold = 0.1
# alpha = 0.1
# herbivore_proportion = 0.3
include("loading_pkg (remove if not needed).jl")
include("kernels.jl")
# Now run first 21 lines from One-click code.jl
########### CREATING DATA ##########
######### raster_with_abundances ########
# Define the size of the raster
size_x, size_y = 20, 20
# Define the x and y dimensions directly
x_dim = X(1:size_x)
y_dim = Y(1:size_y)
# A 125x76 grid of random MyStructs256
raster_matrix = reshape([MyStructs256(SVector{num_species, Float64}(rand(Float64, num_species))) for _ in 1:(size_x * size_y)], size_x, size_y)
# Set first and last rows and columns to zero using the defined Base.zero function
for j in 1:size_y
    raster_matrix[1, j] = zero(MyStructs256{Float64})    # First row
    raster_matrix[end, j] = zero(MyStructs256{Float64})  # Last row
end
for i in 1:size_x
    raster_matrix[i, 1] = zero(MyStructs256{Float64})    # First column
    raster_matrix[i, end] = zero(MyStructs256{Float64})  # Last column
end
raster_with_abundances = Raster(raster_matrix, dims=(x_dim, y_dim))
######### raster_sum ########
# Create a 125x76 grid of Bool values with first and last rows as false, others true
sum_matrix = [i == 1 || i == size_x || j == 1 || j == size_y ? false : true for i in 1:size_x, j in 1:size_y]
# Create the Bool raster with the specified dimensions
raster_sum = Raster(sum_matrix, dims=(x_dim, y_dim))
########## creating k_DA ##########
raster_matrix_k_DA = reshape([MyStructs256(SVector{num_species, Float64}(rand(Float64, num_species))) for _ in 1:(size_x * size_y)], size_x, size_y)
raster_k_DA = Raster(raster_matrix, dims=(x_dim, y_dim))
######## For plotting ##########
# For a heat map we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs256, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > presence_threshold, 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
custom_palette = cgrad([colorant"black", colorant"blue", colorant"yellow", colorant"green", colorant"red"], [0.0, 0.000000000001, 0.33, 0.66, 1.0]);
###### Adjacency Matrix Maker ######
function Adjacency_Matrix_Maker(dim::Int, connectance::Float64, herbivore_proportion::Float64)
    # Initialize an empty adjacency matrix of zeros
    adjacency_matrix = zeros(Float64, dim, dim)

    # Determine the number of herbivores based on the proportion
    num_herbivores = round(Int, herbivore_proportion * dim)
    herbivores = randperm(dim)[1:num_herbivores]  # Randomly select herbivore indices

    # Total number of possible interactions
    total_interactions = dim * (dim - 1)  # Excluding self-interactions

    # Number of non-zero interactions based on connectance
    num_interactions = round(Int, connectance * total_interactions)

    # Randomly choose predator-prey interactions
    for _ in 1:num_interactions
        i, j = rand(1:dim), rand(1:dim)
        
        # Avoid self-interactions
        while i == j
            i, j = rand(1:dim), rand(1:dim)
        end

        if i in herbivores
            # If i is a herbivore, it can only be prey, so it gets a negative interaction (preyed upon)
            adjacency_matrix[j, i] = -1  # Species j eats herbivore i
        else
            # If i is not a herbivore, handle regular predator-prey interactions
            adjacency_matrix[i, j] = 1  # Species i preys on species j
            adjacency_matrix[j, i] = -1 # Species j is prey to species i
        end
    end
    
    return adjacency_matrix, herbivores
end
# A_matrix = Adjacency_Matrix_Maker(num_species, connectance, herbivore_proportion)[1]
##### Interaction Matrix Maker #####
function Interaction_Matrix_Maker(matrix::Array{Float64,2}, sigma::Float64)
    dim = size(matrix, 1)
    # Iterate over the upper triangle of the matrix (since interaction is symmetric)
    for i in 1:dim
        for j in i+1:dim
            if matrix[i, j] != 0
                # Sample a value from a normal distribution with mean 0 and sd sigma
                interaction_value = rand(Normal(0, sigma))
                # Assign the values to both predator-prey pairs
                matrix[i, j] = interaction_value
                matrix[j, i] = -interaction_value
            end
        end
    end
    return matrix
end
# I_matrix = Interaction_Matrix_Maker(A_matrix, sigma)

##################### FILL DIAGONAL ###################################
function fill_diagonal!(mat, val)
    for i in 1:min(size(mat)...)
        mat[i, i] = val
    end
    return mat
end
# I_matrix = fill_diagonal!(I_matrix, self_regulation)
# ##################### RULES ###################################
# A_matrix = Adjacency_Matrix_Maker(num_species, connectance, herbivore_proportion)[1]
# I_matrix = Interaction_Matrix_Maker(A_matrix, sigma)
# I_matrix = fill_diagonal!(I_matrix, self_regulation)
# function GLV(state::MyStructs256, k_DA::MyStructs256)
#     return MyStructs256(
#         SVector{num_species, Float64}(
#             state.a + (state.a .* (k_DA.a - state.a) + ((I_matrix * state.a) .* state.a)) 
#         )
#     )
# end

# biotic = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
#     # if any(isinf, state.a) || any(isnan, state.a)
#     #     @error "state has NA values"
#     #     println(I)
#     # end
#     return MyStructs256(SVector{num_species, Float64}(max.(0.0, GLV(state, k_DA).a)))
# end

# disp = OutwardsDispersal{:state, :state}(
#     formulation = CustomKernel(alpha),
#     distancemethod = AreaToArea(30),
#     maskbehavior = Dispersal.CheckMaskEdges(),
# )

# pepe = (
#     state = Matrix(raster_with_abundances),
#     k_DA = Matrix(raster_k_DA),
# )

# array_output = ResultOutput(
#     pepe, tspan = 1:100;
#     mask = Matrix(raster_sum),
# )

# @time a = sim!(array_output, Ruleset(biotic, disp; boundary = Reflect(), proc = ThreadedCPU()))

# Makie.heatmap(a[end].state)
# Makie.image(a[end].state, colormap=custom_palette, colorrange = (0, num_species))

# MakieOutput(pepe, tspan = 1:100;
#     fps = 10, ruleset = Ruleset(biotic, disp; boundary = Reflect()),
#     mask = Matrix(raster_sum)) do (; layout, frame)

#     # Setup the keys and titles for each plot
#     plot_keys = [:Biomass, :Richness]
#     titles = ["Biomass", "Richness"]

#     # Create axes for each plot and customize them
#     axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:1, j in 1:2]

#     # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
#     for (ax, key, title) in zip(axes, plot_keys, titles)
#         if key == :Biomass
#             Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette, colorrange = (0, sum(num_species)))
#         elseif key == :Richness
#             Makie.image!(ax, frame[:state]; interpolate=false, colormap=custom_palette, colorrange = (0, num_species))
#         end
#         hidexdecorations!(ax; grid=false)
#         hideydecorations!(ax; grid=false)
#         ax.title = title  # Set the title for each axis
#         ax.titlegap[] = 5  # Adjust the title gap to make it smaller
#         ax.titlesize[] = 12 # Set the title font size
#         ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
#         ax.yreversed[] = true
#     end
# end

# ruleset = Ruleset(Life())
# # And the time-span for it to run
# tspan = 1:100
# # Create our own plots with Makie.jl
# output = MakieOutput(rand(Bool, 200, 300); tspan, ruleset) do (; layout, frame)
#     image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
# end