PC = "nicol"
num_species = 256
include("HerpsVsBirmmals.jl")
include("kernels.jl")
include("One-click code.jl")
include("human_footprint.jl")
include("Implicit competition for herbivores.jl")
# include("2010-2100 RasterSeries.jl")

# pepe = (
#     birmmals = Matrix(DA_birmmals_with_abundances),
#     herps = Matrix(DA_herps_with_abundances),
#     k_DA = Matrix(k_DA_hf_additive),
#     birmmals_richness = Matrix(DA_richness_birmmals),
#     herps_richness = Matrix(DA_richness_herps)
# )
DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)
pepe_state = (
    state = Matrix(DA_with_abundances)[50:60, 50:60],
    k_DA = Matrix(k_DA_hf_additive)[50:60, 50:60],
    npp_DA = Matrix(npp_DA)[50:60, 50:60],
    state_richness = Matrix(DA_richness)[50:60, 50:60]
)

##### Parameters #####
self_regulation = 1.0
sigma = 0.1
sigma_comp = 0.0
epsilon = 1.0
beta = 1.0  
assymetry = 1.0
remove_variable(:alfa)
alfa = 0.9
##### Matrices #####
# Trophic
caca = deepcopy(iberian_interact_NA)
full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon, self_regulation, beta))
count = 0
for i in 1:256
    if all(full_IM[:, i] == 0.0)
        count += 1
    end
end
# Competition
competition_NA = deepcopy(iberian_interact_NA)
competition_NA .= 0.0
for i in names(competition_NA, 1), j in names(competition_NA, 2)
    if i in herbivore_names && j in herbivore_names
        competition_NA[i, j] = 1.0
    end
end
full_comp = turn_comp_into_inter(competition_NA, sigma_comp, assymetry)
##########################################################################
# exp(-(1^2) / (2*(alfa^2)))
m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

##### RULES #####
function GLV(state::MyStructs256, k_DA::MyStructs256)
    return MyStructs256(
        SVector{256, Float64}(
            state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a) + ((full_comp * state.a) .* state.a)) 
        )
    )
end

biotic_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return MyStructs256(SVector{256, Float64}(max.(0.0, GLV(state, k_DA).a)))
end

catastrophic_rule = let prob_chaos = 0.1
    Neighbors{:cata_layer, :cata_layer}(Moore(1)) do data, neighborhood, cata_layer, I
        event = rand() <= prob_chaos ? true : false
        if event & isone(cata_layer)
            return false   
        elseif event & !iszero(cata_layer)
            return true
        else
            return cata_layer
        end
    end
end

neighbors_rule = let prob_chaos = 0.01
    Neighbors{:state, :state}(Moore(1)) do data, neighborhood, state, I
        # println(typeof(DG.mask(data)[I[1], I[2]]))
        # println(I)
        event = rand() <= prob_chaos ? true : false
        if event
            return MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        else
            return state
        # elseif any(isone, neighborhood)
            # return state - MyStructs256(SVector{256, Float64}(fill(0.1, 256)))
        end
    end
end

######## DISPERSAL ########
scaling_vector = fill(1.0, 49)
scaling_vector = append!(scaling_vector, fill(1.0, 207))
remix_outdisp = OutwardsDispersalRemix{:state, :state}(
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);

indisp = InwardsDispersal{:state, :state}(;
    formulation=ExponentialKernel(),
    distancemethod=AreaToArea(30)
);

##### MAKIE STATE #####
array_output = ArrayOutput(
    pepe_state; tspan = 1:1000,
    # mask = Matrix(DA_sum)
)

@time c = sim!(array_output, Ruleset(biotic_GLV; boundary = Reflect(), proc = ThreadedCPU()))

num_timesteps = length(c)
# Choose any time step and cell to get the number of species
i = 10  # Specify the row index of the cell of interest
j = 10  # Specify the column index of the cell of interest
num_species = length(c[1].state[i, j].a)

# Rows: Time steps
# Columns: Species abundances
abundances_over_time = Array{Float64}(undef, num_timesteps, num_species)

for t in 1:num_timesteps
    # Extract the abundance vector from the cell at time t
    abundances_over_time[t, :] = c[t].state[i, j].a
end
time = collect(1:num_timesteps)

# fig = Figure()
# ax = Axis(fig[1, 1], xlabel="Time", ylabel="Abundance", title="Species Abundances Over Time at Cell ($i, $j)")

# num_species = size(abundances_over_time, 2)

# for s in 1:num_species
#     # Extract the abundances for species s
#     abundances = abundances_over_time[:, s]
#     # Plot the line for species s
#     lines!(ax, time, abundances, label=species_names[s])
# end

##################################################################################
# Number of time steps in your simulation
num_timesteps = length(c)

# Get the number of species from any cell
i_sample, j_sample = 1, 1  # Sample indices
num_species = length(c[1].state[i_sample, j_sample].a)

# Assuming you have species names
species_names = ["Species $(s)" for s in 1:num_species]  # Replace with actual names if available

# Generate all cell indices in the 11×11 grid
cells = [(i, j) for i in 1:11, j in 1:11]

using Random
# Set a seed for reproducibility (optional)
Random.seed!(123)

# Sample 20 unique cells without replacement
sampled_cells = sample(cells, 20, replace=false)

# Determine the layout for 20 plots (e.g., 4 rows × 5 columns)
nrows = 4
ncols = 5
fig = Figure(resolution=(1000, 800))

# Convert time to a vector (Makie expects arrays)
time = collect(1:num_timesteps)

# Loop over each sampled cell
for idx in 1:length(sampled_cells)
    (i, j) = sampled_cells[idx]
    
    # Extract abundances over time for cell (i, j)
    abundances_over_time = Array{Float64}(undef, num_timesteps, num_species)
    for t in 1:num_timesteps
        abundances_over_time[t, :] = c[t].state[i, j].a
    end
    
    # Determine subplot position
    row = div(idx - 1, ncols) + 1
    col = mod(idx - 1, ncols) + 1
    
    # Create an axis in the figure at the appropriate position
    ax = Axis(fig[row, col],
              xlabel="Time",
              ylabel="Abundance",
              title="Cell ($i, $j)")
    
    # Plot abundances for all species
    for s in 1:num_species
        abundances = abundances_over_time[:, s]
        lines!(ax, time, abundances)
    end
    
end

# Display the figure
fig
