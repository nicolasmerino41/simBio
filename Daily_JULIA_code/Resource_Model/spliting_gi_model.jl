# Import necessary packages
using DifferentialEquations
using CairoMakie  # For plotting

begin
    
# Define common parameters
H0_common = 100.0   # Common characteristic density
m_common = 0.1      # Common mortality rate
mu_value = 1.0      # Common interaction coefficient
NPP = 1000.0        # Net Primary Production

# Number of species
num_species = 5

# Assign p_i values (proportion of resources)
# Equal allocation
p_values = [1.0 / num_species for _ in 1:num_species]
# Random allocation
p_values = rand(num_species)
p_values ./= sum(p_values)  # Normalize so that they sum to 1

# Define the Species struct
mutable struct Species
    id::Int
    name::String
    m::Float64        # Mortality rate (m_i)
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (g_i)
    p::Float64        # Proportion of resources allocated to species i (p_i)
    mu::Float64       # Interaction coefficient (average mu_ij)
end
# Constructor for Species
Species(;id, name, m, H0, H_init, g=0.0, p=0.0, mu=0.0) = Species(id, name, m, H0, H_init, g, p, mu)

# Function to create the interaction matrix
function create_interaction_matrix(num_species::Int, mu_value::Float64)
    mu_matrix = fill(mu_value, num_species, num_species)
    for i in 1:num_species
        mu_matrix[i, i] = 0.0  # No self-interaction
    end
    return mu_matrix
end

# Function to calculate g_i
function calculate_g_i(p_i, NPP, F_i)
    m_i = m_common
    return m_i * ( (1 + sqrt(1 + (4 * p_i * NPP) / F_i)) / 2 )
end

# Initialize species list
species_list = Species[]

# Create species with parameters
for i in 1:num_species
    p_i = p_values[i]
    H0_i = H0_common
    F_i = m_common * H0_i  # F_i = m_i * H_i^0
    g_i = calculate_g_i(p_i, NPP, F_i)
    sp = Species(
        id = i,
        name = "Species $i",
        m = m_common,
        H0 = H0_i,
        H_init = H0_i,  # Initial abundance equal to characteristic density
        g = g_i,
        p = p_i,
        mu = mu_value
    )
    push!(species_list, sp)
end

# Create the interaction matrix
mu_matrix = create_interaction_matrix(num_species, mu_value)

# Initial abundances
u0 = [sp.H_init for sp in species_list]

# Time span for the simulation
tspan = (0.0, 100.0)

# Parameters for the ODEProblem
p = (species_list, mu_matrix)

# Define the ecosystem dynamics function
function ecosystem_dynamics!(du, u, p, t)
    species_list, mu_matrix = p
    num_species = length(species_list)
    du .= 0.0  # Initialize derivatives

    # Calculate total growth rate G
    G = sum([sp.g for sp in species_list])

    for i in 1:num_species
        sp = species_list[i]
        H_i = u[i]
        m_i = sp.m
        g_i = sp.g
        H_i0 = sp.H0

        # Compute interaction term
        interaction = 0.0
        for j in 1:num_species
            H_j = u[j]
            mu_ij = mu_matrix[i, j]
            interaction += mu_ij * H_j
        end

        # Compute the derivative
        du[i] = H_i * m_i * ( (g_i / m_i - 1) - (H_i + interaction) / H_i0 )
    end
end

# Define the ODE problem
prob = ODEProblem(ecosystem_dynamics!, u0, tspan, p)

# Solve the ODE
sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

# Extract times and solutions
times = sol.t
u_array = hcat(sol.u...)  # Each column is the state at a time point

# Plot the abundances over time
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Abundance", title = "Species Abundances")

# Plot each species
for i in 1:num_species
    lines!(ax, times, u_array[i, :], label = species_list[i].name)
end

axislegend(ax)

display(fig)
end