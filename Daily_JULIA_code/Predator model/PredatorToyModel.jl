# Import necessary packages
using DifferentialEquations
using CairoMakie  # For plotting
using Distributions  # For drawing from distributions

#######################################################
################## Parameters #########################
begin
# Set parameters for distributions
num_predators = 1
num_herbivores = 10

# Means and standard deviations for masses
pred_mass_mean = 1.0
pred_mass_sd = 0.1
herb_mass_mean = 0.1
herb_mass_sd = 0.01

# Means and standard deviations for abundances
pred_abundance_mean = 10.0
pred_abundance_sd = 1.0
herb_abundance_mean = 100.0
herb_abundance_sd = 10.0

# Means and standard deviations for attack rates and handling times
pred_a0_mean = 0.1
pred_a0_sd = 0.01
pred_h0_mean = 0.1
pred_h0_sd = 0.01

# Means and standard deviations for mortality rate coefficients
m0_mean = 0.1
m0_sd = 0.01

# Connectance parameter
connectance = 0.5  # Probability that a predator consumes a prey

# Other parameters
e_predator = 0.5  # Assimilation efficiency for predators
d_predator = 100000.0  # Self-limitation coefficient for predators
#####################################################################
#####################################################################

# Define the Species struct with additional parameters
mutable struct Species
    id::Int
    name::String
    mass::Float64
    abundance::Float64
    mortality::Float64
    is_predator::Bool
    prey_ids::Vector{Int}  # List of prey species IDs
    a0::Float64  # Attack rate coefficient (for predators)
    h0::Float64  # Handling time coefficient (for predators)
    m0::Float64  # Mortality rate coefficient
end

# Define the Parameters struct for global parameters
struct Parameters
    alpha::Float64  # Attack rate exponent
    beta::Float64   # Handling time exponent
    e::Float64      # Assimilation efficiency
    d::Float64      # Self-limitation coefficient
end
# Constructor for Parameters
Parameters(;alpha, beta, e, d) = Parameters(alpha, beta, e, d)

# Function to create the species list with parameters drawn from distributions
function create_species_list(num_predators, num_herbivores;
                             pred_mass_mean, pred_mass_sd,
                             herb_mass_mean, herb_mass_sd,
                             pred_abundance_mean, pred_abundance_sd,
                             herb_abundance_mean, herb_abundance_sd,
                             pred_a0_mean, pred_a0_sd,
                             pred_h0_mean, pred_h0_sd,
                             m0_mean, m0_sd)
    species_list = Vector{Species}()
    species_id = 1  # Unique species ID counter

    # Create herbivores first
    for i in 1:num_herbivores
        mass = rand(Normal(herb_mass_mean, herb_mass_sd))
        abundance = rand(Normal(herb_abundance_mean, herb_abundance_sd))
        m0 = rand(Normal(m0_mean, m0_sd))
        sp = Species(
            species_id,
            "Herbivore $i",
            mass,
            abundance,
            0.0,      # Mortality (will be calculated later)
            false,    # is_predator
            [],       # Prey IDs (none for herbivores)
            0.0,      # a0 (not used for herbivores)
            0.0,      # h0 (not used for herbivores)
            m0        # Mortality rate coefficient
        )
        push!(species_list, sp)
        species_id += 1
    end

    # Then create predators
    for i in 1:num_predators
        mass = rand(Normal(pred_mass_mean, pred_mass_sd))
        abundance = rand(Normal(pred_abundance_mean, pred_abundance_sd))
        m0 = rand(Normal(m0_mean, m0_sd))
        a0 = rand(Normal(pred_a0_mean, pred_a0_sd))
        h0 = rand(Normal(pred_h0_mean, pred_h0_sd))
        sp = Species(
            species_id,
            "Predator $i",
            mass,
            abundance,
            0.0,      # Mortality (will be calculated later)
            true,     # is_predator
            [],       # Prey IDs (to be assigned later)
            a0,       # Attack rate coefficient
            h0,       # Handling time coefficient
            m0        # Mortality rate coefficient
        )
        push!(species_list, sp)
        species_id += 1
    end

    return species_list
end

# Function to calculate mortality using species-specific m0
function calculate_mortality(mass, m0)
    return m0 * mass^(-0.25)
end

# Function to create interaction matrix based on connectance
function create_interaction_matrix(species_list, connectance, id_to_index)
    num_species = length(species_list)
    prey_ids = [sp.id for sp in species_list if !sp.is_predator]
    predator_ids = [sp.id for sp in species_list if sp.is_predator]

    # Initialize interaction matrix
    interaction_matrix = zeros(Bool, num_species, num_species)

    # For each predator, decide which prey it eats
    for predator_id in predator_ids
        predator_index = id_to_index[predator_id]
        predator = species_list[predator_index]
        prey_ids_eaten = Int[]
        for prey_id in prey_ids
            if rand() < connectance
                # Predator eats prey
                prey_ids_eaten = vcat(prey_ids_eaten, prey_id)
                prey_index = id_to_index[prey_id]
                interaction_matrix[predator_index, prey_index] = true
            end
        end
        # Assign prey_ids to the predator
        predator.prey_ids = prey_ids_eaten
    end

    return interaction_matrix
end

# Predator dynamics function using species-specific parameters
function predator_dynamics!(du, u, p, t)
    species_list, params, id_to_index = p
    num_species = length(species_list)
    du .= 0.0

    for i in 1:num_species
        sp = species_list[i]
        if sp.is_predator
            P = u[i]
            consumption = 0.0
            for prey_id in sp.prey_ids
                prey_index = id_to_index[prey_id]
                prey_sp = species_list[prey_index]
                H = u[prey_index]
                a = sp.a0 * sp.mass^(params.alpha)
                h = sp.h0 * prey_sp.mass^(params.beta)
                f = (a * H) / (1 + a * h * H)
                consumption += f
            end
            du[i] = P * (params.e * consumption - params.d * P - sp.mortality)
        else
            du[i] = 0.0  # Herbivores are static
        end
    end
end

# Create species list with parameters drawn from distributions
species_list = create_species_list(
    num_predators, num_herbivores;
    pred_mass_mean=pred_mass_mean, pred_mass_sd=pred_mass_sd,
    herb_mass_mean=herb_mass_mean, herb_mass_sd=herb_mass_sd,
    pred_abundance_mean=pred_abundance_mean, pred_abundance_sd=pred_abundance_sd,
    herb_abundance_mean=herb_abundance_mean, herb_abundance_sd=herb_abundance_sd,
    pred_a0_mean=pred_a0_mean, pred_a0_sd=pred_a0_sd,
    pred_h0_mean=pred_h0_mean, pred_h0_sd=pred_h0_sd,
    m0_mean=m0_mean, m0_sd=m0_sd
)

# Build ID to index mapping
id_to_index = Dict{Int, Int}()
for (idx, sp) in enumerate(species_list)
    id_to_index[sp.id] = idx
end

# Create interaction matrix and assign prey_ids to predators
interaction_matrix = create_interaction_matrix(species_list, connectance, id_to_index)

# Calculate mortality rates for all species
for sp in species_list
    sp.mortality = calculate_mortality(sp.mass, sp.m0)
end

# Define global parameters
params = Parameters(
    alpha = 0.01,  # Attack rate exponent
    beta = 0.1,   # Handling time exponent
    e = e_predator,      # Assimilation efficiency
    d = d_predator      # Self-limitation coefficient
)

# Initial abundances
u0 = [sp.abundance for sp in species_list]

# Time span for the simulation
tspan = (0.0, 100.0)

# Define problem parameters
p = (species_list, params, id_to_index)

# Define the ODE problem
prob = ODEProblem(predator_dynamics!, u0, tspan, p)

# Define extinction callbacks and positive domain
callbacks = []
push!(callbacks, PositiveDomain())
cb = CallbackSet(callbacks...)

# Solve the ODE
sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

# Extract times and solutions
times = sol.t
u_array = hcat(sol.u...)  # Each column is the state at a time point

# Get indices of predators
predator_indices = [i for (i, sp) in enumerate(species_list) if sp.is_predator]

# Extract predator abundances over time
predator_abundances = u_array[predator_indices, :]

# Plot using Makie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Abundance", title = "Predator Dynamics")

# Plot each predator
for (i, idx) in enumerate(predator_indices)
    lines!(ax, times, predator_abundances[i, :], label=species_list[idx].name)
end

axislegend(ax)

display(fig)
end