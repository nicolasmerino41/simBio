# Load Required Packages
using StaticArrays, Random, Makie, Distributions, DynamicGrids, Rasters, WGLMakie, Dispersal

######## PARAMETERS ########
num_species = 100
connectance = 0.11
sigma = 0.1
self_regulation = -1.0
carrying_capacity = 1.0
presence_threshold = 0.1
alfa = 0.2
herbivore_proportion = 0.6
mortality_rate = 0.1  # m_i for all species
consumption_rate = 0.01  # u_i for all species
efficiency = 0.1
# Define the number of herbivores and predators
num_herbivores = round(Int, herbivore_proportion * num_species)
num_predators = num_species - num_herbivores
herbivores = 1:num_herbivores
predators = num_herbivores+1:num_species

######## STRUCT DEFINITION ########
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end
########### CREATING DATA ##########
######### raster_with_abundances ########
# Define the size of the raster
size_x, size_y = 100, 100

# Define the x and y dimensions directly
x_dim = X(1:size_x)
y_dim = Y(1:size_y) 

# Initialize abundance grid with random values
initial_abundance = [MyStructs256(@SVector(rand(Float64, num_species))) for _ in 1:(size_x * size_y)]
raster_matrix = reshape(initial_abundance, size_x, size_y)

raster_with_abundances = Raster(raster_matrix, dims=(x_dim, y_dim))

######### raster_sum ########
# Create a grid indicating active cells (inner cells are active)
sum_matrix = [i == 1 || i == size_x || j == 1 || j == size_y ? false : true for i in 1:size_x, j in 1:size_y]
raster_sum = Raster(sum_matrix, dims=(x_dim, y_dim))

########## creating g_i ##########
# Initialize g_i grid with random suitability values between 0 and 1
g_i_grid = [@SVector(rand(Uniform(0.0, 1.0), num_species)) for _ in 1:(size_x * size_y)]
raster_g_i = reshape(g_i_grid, size_x, size_y)

########## creating npp_grid ##########
# Initialize NPP grid with a constant or random renewal rates
npp_grid = rand([10.0, 20.0, 30.0], size_x, size_y)  # Constant NPP for simplicity
########## creating resource_grid ##########
# Initialize resource grid with a constant or random renewal rates
resource_grid = rand([100.0, 200.0], size_x, size_y)  # Constant NPP for simplicity

########## creating u_i ##########
# Consumption rates for herbivores (random between 0.05 and 0.15)
u_i = @SVector(fill(consumption_rate, num_species))

########## creating e_i ##########
# Efficiency of converting consumed prey into predator biomass (random between 0.1 and 0.3)
e_i = @SVector([efficiency for _ in 1:num_species])

########## Interaction Matrix Maker ##########
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

# Create interaction matrices
A_matrix = Adjacency_Matrix_Maker(num_species, connectance, herbivore_proportion)[1]
I_matrix = Interaction_Matrix_Maker(A_matrix, sigma)

########## DEFINE RULES ##########
# Resource Update Rule
resource_rule = Cell{Tuple{:state, :npp, :resources}, :resources}() do data, (state, npp, resources), I
    dR_dt = npp - sum(u_i[herbivores] .* state.a[herbivores] .* resources)
    new_R = resources + dR_dt # Assuming dt=0.1
    return new_R
end

# Herbivore Update Rule with Predation
herbivore_rule = Cell{Tuple{:state, :resources, :g_i_grid}, :state}() do data, (state, resources, g_i_grid), I
    # Extract herbivore abundances as a Vector
    herbivores_abundance = Vector{Float64}(state.a[herbivores])
    
    # Extract g_i values for herbivores in the current cell
    g_i_values = Vector{Float64}(g_i_grid[herbivores])
    
    # Extract predator abundances in the current cell
    predators_abundance = Vector{Float64}(state.a[predators])
    
    # Extract interaction coefficients (alpha_ij) for predators and herbivores
    alpha_ij = I_matrix[predators, herbivores]  # Matrix of size (num_predators, num_herbivores)
    
    # Calculate predation: sum_p (alpha_ij * P_p) for each herbivore h
    # This results in a Vector of size (num_herbivores,)
    predation = alpha_ij' * predators_abundance  # Vector{Float64}(num_herbivores)
    
    # Calculate consumption: g_i * H_i * R
    consumption = g_i_values .* herbivores_abundance .* resources
    
    # Update herbivore abundances: H_i + consumption - predation - mortality
    new_abundance = herbivores_abundance .+ consumption .- predation .- (mortality_rate .* herbivores_abundance)
    
    # Check for Inf values
    if any(isinf, new_abundance)
        println("Warning: Herbivore abundance at cell $(data.position) reached infinity.")
    end
    
    # Ensure no negative abundances
    new_abundance = max.(new_abundance, 0.0)
    
    # Update the full abundance vector
    all_abundances = Vector{Float64}(state.a)
    all_abundances[herbivores] = new_abundance
    
    # Return the updated state
    return MyStructs256(SVector{num_species, Float64}(all_abundances))
end



# Predator Update Rule
predator_rule = Cell{Tuple{:state, :g_i_grid}, :state}() do data, (state, g_i_grid_cell), I
    # Extract predator abundances as a Vector
    predators_abundance = Vector{Float64}(state.a[predators])
    
    # Extract g_i values for predators in the current cell
    g_i_values = Vector{Float64}(g_i_grid_cell[predators])
    
    # Extract prey abundances (herbivores) as a Vector
    prey_abundance = Vector{Float64}(state.a[herbivores])
    
    # Extract interaction coefficients (alpha_ij) for predators and prey
    alpha_ij = I_matrix[predators, herbivores]  # Matrix of size (num_predators, num_herbivores)
    
    # Calculate the sum over j of alpha_ij * X_j for each predator i
    sum_alphaX = alpha_ij * prey_abundance  # Vector of size (num_predators,)
    
    # Calculate growth term: g_i * (e_i * sum_alphaX) * P_i
    growth = g_i_values .* (e_i[predators] .* sum_alphaX) .* predators_abundance
    
    # Update predator abundances: P_i + growth - mortality
    new_predators_abundance = predators_abundance .+ growth .- (mortality_rate .* predators_abundance)
    
    # Check for Inf values
    if any(isinf, new_predators_abundance)
        println("Warning: Predator abundance at cell $(data.position) reached infinity.")
    end
    
    # Ensure no negative abundances
    new_predators_abundance = max.(new_predators_abundance, 0.0)
    
    # Update the full abundance vector
    all_abundances = Vector{Float64}(state.a)
    all_abundances[predators] = new_predators_abundance
    
    # Return the updated state
    return MyStructs256(SVector{num_species, Float64}(all_abundances))
end

# Define biotic rules
biotic = Cell{Tuple{:state}, :state}() do data, (state,), _
    # Update Resource
    state = Resource_Update(state, npp_grid[data.position[1], data.position[2]], u_i)
    # Update Herbivores
    state = Herbivore_Update(state, raster_g_i[data.position[1], data.position[2]], c_ij=zeros(Float64, num_species, num_species), m_i=@SVector(fill(mortality_rate, num_species)))
    # Update Predators
    state = Predator_Update(state, I_matrix, e_i, @SVector(fill(mortality_rate, num_species)))
    return MyStructs256(state.a, state.R)
end

disp = OutwardsDispersal{:state, :state}(
    formulation = CustomKernel(alfa),
    distancemethod = AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges(),
)

# Initial state
pepe = (
    state = Matrix(raster_with_abundances),
    npp = Matrix(npp_grid),
    resources = Matrix(resource_grid),
    g_i_grid = Matrix(raster_g_i),
)

# Define mask
mask = Matrix(raster_sum)

# Define output
array_output = ResultOutput(
    pepe, tspan = 1:10;
    mask = mask,
)
p = sim!(array_output, Ruleset(resource_rule, herbivore_rule, predator_rule))

function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end

# Run Simulation and Visualization
MakieOutput(pepe, tspan = 1:100;
    fps = 10, ruleset = Ruleset(resource_rule, herbivore_rule, predator_rule),
    mask = mask) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:Biomass, :NPP, :Resources]
    titles = ["Biomass", "NPP", "Resources"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[j]) for j in 1:3, i in 1:1]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :Biomass
            heatmap!(ax, frame[:state])
        elseif key == :NPP
            heatmap!(ax, frame[:npp])
        elseif key == :Resources
            heatmap!(ax, frame[:resources])
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.titlegap[] = 5  # Adjust the title gap to make it smaller
        ax.titlesize[] = 12 # Set the title font size
        ax.yreversed[] = true
    end
end
