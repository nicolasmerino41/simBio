PC = "nicol"
using Pkg
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
using GLM
using CSV, DataFrames
using DifferentialEquations
using Distributions
using Random

# Define the Herbivore struct
mutable struct Herbivore
    m::Float64        # Mortality rate
    H0::Float64       # Characteristic density (H_i^0)
    H_init::Float64   # Initial abundance (H_i(0))
    g::Float64        # Growth rate (to be calculated)
end

# Constructor for Herbivore
Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(m, H0, H_init, g)

# Function to create herbivores list
function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
    H0_mean::Float64, H0_sd::Float64, H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
    herbivores_list = Herbivore[]
    for i in 1:num_herbivores
        m = abs(rand(Normal(m_mean, m_sd)))           # Mortality rate
        H0 = abs(rand(Normal(H0_mean, H0_sd)))        # Characteristic density
        H_init = H0  # Initial abundance set to characteristic density
        push!(herbivores_list, Herbivore(m=m, H0=H0, H_init=H_init))
    end
    return herbivores_list
end

# Function to calculate growth rates based on NPP
function calculate_growth_rates(herbivore_list, NPP, mu)
    S_star = length(herbivore_list)
    F_list = [sp.H0 * sp.m for sp in herbivore_list]
    competition_numerator = 1 + mu * (S_star - 1)
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
    c::Float64        # Self-regulation coefficient
end

# Predator constructor
Predator(; m::Float64, a::Float64, h::Float64, e::Float64, P_init::Float64, c::Float64) = Predator(m, a, h, e, P_init, c)

# Function to create a list of predators
function create_predator_list(num_predators::Int; m_mean::Float64, m_sd::Float64, a_mean::Float64, a_sd::Float64,
    h_mean::Float64, h_sd::Float64, e_mean::Float64, e_sd::Float64, P_init_mean::Float64, P_init_sd::Float64,
    c_mean::Float64, c_sd::Float64)
    predator_list = Predator[]
    for _ in 1:num_predators
        m = abs(rand(Normal(m_mean, m_sd)))             # Mortality rate
        a = abs(rand(Normal(a_mean, a_sd)))             # Attack rate
        h = abs(rand(Normal(h_mean, h_sd)))             # Handling time
        e = abs(rand(Normal(e_mean, e_sd)))             # Conversion efficiency
        P_init = abs(rand(Normal(P_init_mean, P_init_sd)))  # Initial abundance
        c = abs(rand(Normal(c_mean, c_sd)))             # Self-regulation coefficient
        push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init, c=c))
    end
    return predator_list
end

# Function to generate the interaction matrix
function generate_interaction_matrix(num_predators::Int, num_prey::Int, connectivity::Float64)
    IM = zeros(Bool, num_predators, num_prey)
    for k in 1:num_predators
        for i in 1:num_prey
            IM[k, i] = rand() < connectivity
        end
    end
    return IM
end

# Ecosystem dynamics function
function ecosystem_dynamics!(du, u, p, t)
    herbivore_list, beta_matrix, predator_list, IM = p
    S_star = length(herbivore_list)
    num_predators = length(predator_list)
    H = u[1:S_star]         # Herbivore densities
    P = u[S_star+1:end]     # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        sp = herbivore_list[i]
        competition = sum(beta_matrix[i, :] .* H) / sp.H0
        predation = 0.0
        for k in 1:num_predators
            if IM[k, i]
                pred = predator_list[k]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])
                predation += P[k] * f_ki
            end
        end
        du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
    end

    # Predator dynamics with self-regulation
    for k in 1:num_predators
        pred = predator_list[k]
        ingestion = 0.0
        for i in 1:S_star
            if IM[k, i]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])
                ingestion += f_ki
            end
        end
        du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
    end

    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end
# Import SpinLock for thread-safety
using Base.Threads: SpinLock

# Initialize a SpinLock for file access
const file_lock = SpinLock()

# Function to write results safely
function write_results_safely(results_row, csv_filename)
    lock(file_lock)  # Lock the file for thread-safe access
    try
        if isfile(csv_filename)
            CSV.write(csv_filename, results_row, append=true)
        else
            CSV.write(csv_filename, results_row)
        end
    finally
        unlock(file_lock)  # Unlock the file
    end
end

# Updated run_simulation function with thread-safe CSV writing
function run_simulation(num_herbivores, num_predators, mu, NPP, m_mean_h, m_mean_p, H0_mean, connectivity, c_mean_p)
    # Your simulation logic remains the same
    herbivore_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean, H0_sd=H0_mean/10)
    calculate_growth_rates(herbivore_list, NPP, mu)

    S_star = length(herbivore_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0
    end

    predator_list = create_predator_list(num_predators; m_mean=m_mean_p, m_sd=0.02, a_mean=0.01, a_sd=0.0001,
        h_mean=0.1, h_sd=0.01, e_mean=0.1, e_sd=0.01, P_init_mean=5.0, P_init_sd=1.0, c_mean=c_mean_p, c_sd=0.01)

    IM = generate_interaction_matrix(num_predators, S_star, connectivity)
    H_init_values = [sp.H_init for sp in herbivore_list]
    P_init_values = [pred.P_init for pred in predator_list]
    u_init = vcat(H_init_values, P_init_values)
    tspan = (0.0, 200.0)

    p = (herbivore_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    herbivore_biomass = sum(sol[1:S_star, end])
    predator_biomass = sum(sol[S_star+1:end, end])
    total_biomass = herbivore_biomass + predator_biomass
    total_biomass_vs_npp = total_biomass / NPP
    num_extinct_herbivores = count(sol[1:S_star, end] .<= 1.0)
    num_extinct_predators = count(sol[S_star+1:end, end] .<= 1.0)
    are_there_extinctions = num_extinct_herbivores > 0 || num_extinct_predators > 0
    prop_of_sp_extinct = (num_extinct_herbivores + num_extinct_predators) / (num_herbivores + num_predators)

    # Create results row
    results_row = DataFrame(
        num_herbivores=num_herbivores,
        num_predators=num_predators,
        mu=mu,
        NPP=NPP,
        m_mean_h=m_mean_h,
        m_mean_p=m_mean_p,
        H0_mean=H0_mean,
        connectivity=connectivity,
        c_mean_p=c_mean_p,
        total_biomass=total_biomass,
        num_extinct_herbivores=num_extinct_herbivores,
        num_extinct_predators=num_extinct_predators,
        herbivore_biomass=herbivore_biomass,
        predator_biomass=predator_biomass,
        are_there_extinctions=are_there_extinctions,
        prop_of_sp_extinct=prop_of_sp_extinct,
        total_biomass_vs_npp=total_biomass_vs_npp
    )

    # Write results using the thread-safe function
    csv_filename = "Daily JULIA code/Resource Model/simulation_results.csv"
    write_results_safely(results_row, csv_filename)

    println("go")
    return 1.0
end

# Simulation parameters
num_herbivores_range = 1:16
num_predators_range = 0:10
mu_values = range(0.0, 1.0, length=10)
NPP_values = range(10.0, 10000.0, length=10)
m_mean_h_values = range(0.001, 0.6, length=10)
m_mean_p_values = range(0.001, 0.6, length=10)
H0_mean_values = range(1.0, 1000.0, length=10)
connectivity_values = range(0.0, 1.0, length=10)
c_mean_p_values = range(0.0, 1.0, length=10)

# Use Threads.@threads to parallelize the loop
Threads.@threads for num_herbivores in num_herbivores_range
    for num_predators in num_predators_range
        for mu in mu_values
            for NPP in NPP_values
                for m_mean_h in m_mean_h_values
                    for m_mean_p in m_mean_p_values
                        for H0_mean in H0_mean_values
                            for connectivity in connectivity_values
                                for c_mean_p in c_mean_p_values
                                    run_simulation(num_herbivores, num_predators, mu, NPP, m_mean_h, m_mean_p, H0_mean, connectivity, c_mean_p)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
    

