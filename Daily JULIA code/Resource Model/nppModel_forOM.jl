# Required Packages
using DifferentialEquations
using Distributions
using Random

# OpenMOLE-provided variables (assumed to be defined)
# mu                 # Competition coefficient (0 to 1)
# NPP                # Net Primary Productivity (10 to 10000)
# num_predators      # Number of predator species
# num_herbivores     # Number of herbivore species
# H0_mean            # Mean characteristic density of herbivores
# Optionally, you can define 'seed' for reproducibility

# Set random seed for reproducibility (optional)
# Random.seed!(seed)

# Derived parameters
H0_mean_aprox = H0_mean  # H0_mean provided by OpenMOLE
connectivity = 1.0       # You can also make this an input if needed

# Define the Herbivores struct
mutable struct Herbivore
    m::Float64        # Mortality rate
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (to be calculated)
end

# Constructor for Herbivore
Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(m, H0, H_init, g)

# Function to create herbivore_list
function create_herbivore_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                               H0_mean::Float64, H0_sd::Float64)
    herbivore_list = []
    for i in 1:num_herbivores
        m = max(0.01, rand(Normal(m_mean, m_sd)))          # Mortality rate
        H0 = max(1.0, rand(Normal(H0_mean, H0_sd)))        # Characteristic density
        H_init = H0  # Initial abundance set to characteristic density
        push!(herbivore_list, Herbivore(m=m, H0=H0, H_init=H_init))
    end
    return herbivore_list
end

# Function to calculate growth rates based on NPP
function calculate_growth_rates(herbivore_list, NPP, mu)
    S_star = length(herbivore_list)
    # Calculate Fi for each herbivore
    F_list = [sp.H0 * sp.m for sp in herbivore_list]
    # Calculate the numerator of the competition term
    competition_numerator = 1 + mu * (S_star - 1)
    # Calculate gi for each herbivore
    for (i, sp) in enumerate(herbivore_list)
        Fi = F_list[i]
        sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
    end
end

# Predator struct definition
mutable struct Predator
    m::Float64        # Mortality rate
    a::Float64        # Attack rate
    h::Float64        # Handling time
    e::Float64        # Conversion efficiency
    P_init::Float64   # Initial abundance
end

# Predator constructor
Predator(; m::Float64, a::Float64, h::Float64, e::Float64, P_init::Float64) = Predator(m, a, h, e, P_init)

# Function to create a list of predators
function create_predator_list(num_predators::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                              a_mean::Float64=0.001, a_sd::Float64=0.0001,
                              h_mean::Float64=0.1, h_sd::Float64=0.01,
                              e_mean::Float64=0.1, e_sd::Float64=0.01)
    predator_list = []
    for _ in 1:num_predators
        m = max(0.01, rand(Normal(m_mean, m_sd)))              # Mortality rate
        a = max(0.0001, rand(Normal(a_mean, a_sd)))            # Attack rate
        h = max(0.01, rand(Normal(h_mean, h_sd)))              # Handling time
        e = max(0.01, rand(Normal(e_mean, e_sd)))              # Conversion efficiency
        P_init = max(0.1, rand(Normal(5.0, 1.0)))              # Initial abundance
        push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init))
    end
    return predator_list
end

# Function to generate the interaction matrix
function generate_interaction_matrix(num_predators::Int, num_prey::Int, connectivity::Float64)
    IM = zeros(Bool, num_predators, num_prey)
    for k in 1:num_predators
        for i in 1:num_prey
            IM[k, i] = rand() < connectivity  # Assign 1 with probability equal to connectivity
        end
    end
    return IM
end

# Ecosystem dynamics function with IM
function ecosystem_dynamics!(du, u, p, t)
    herbivore_list, beta_matrix, predator_list, IM = p
    S_star = length(herbivore_list)
    num_predators = length(predator_list)
    H = u[1:S_star]  # Herbivore densities
    P = u[S_star+1:end]  # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        sp = herbivore_list[i]
        # Compute competition term
        competition = sum(beta_matrix[i, :] .* H) / sp.H0
        # Compute predation term
        predation = 0.0
        for k in 1:num_predators
            if IM[k, i]
                pred = predator_list[k]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                predation += P[k] * f_ki
            end
        end
        # Compute derivative for herbivores
        du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
    end

    # Predator dynamics
    for k in 1:num_predators
        pred = predator_list[k]
        ingestion = 0.0
        for i in 1:S_star
            if IM[k, i]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                ingestion += f_ki
            end
        end
        # Compute derivative for predators
        du_P[k] = P[k] * (pred.e * ingestion - pred.m)
    end

    # Assign derivatives
    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end

# Main simulation function
function run_simulation()
    # Create herbivore list and calculate growth rates
    herbivore_list = create_herbivore_list(
        num_herbivores,
        H0_mean=H0_mean_aprox,
        H0_sd=H0_mean_aprox / 10
    )
    calculate_growth_rates(herbivore_list, NPP, mu)

    # Create beta_matrix
    S_star = length(herbivore_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0  # Self-competition
    end

    # Create predator list
    predator_list = create_predator_list(num_predators)

    # Generate the interaction matrix
    IM = generate_interaction_matrix(num_predators, S_star, connectivity)

    # Initial conditions for herbivores and predators
    H_init_values = [sp.H_init for sp in herbivore_list]
    P_init_values = [pred.P_init for pred in predator_list]

    # Combined initial conditions for the system
    u_init = vcat(H_init_values, P_init_values)

    # Define the time span for simulation
    tspan = (0.0, 200.0)

    # Define the ODE problem
    p = (herbivore_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Solve the ODE
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract herbivore data
    S_star = length(herbivore_list)
    H_array = sol[1:S_star, :]

    # Calculate total biomass at the end
    total_biomass = sum(H_array[:, end])

    # Return the total biomass
    return total_biomass/NPP
end

# Run the simulation and obtain the output
total_biomass = run_simulation()

# Output the result (OpenMOLE will capture the variable 'total_biomass')
println("total_biomass = ", total_biomass)
