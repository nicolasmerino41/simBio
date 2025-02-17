using Evolutionary, DifferentialEquations, Random, Logging

const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0
const art_pi              = false

cell = 1
local_i, local_j = idx[cell][1], idx[cell][2]

# Gather cell data
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
local_S, local_R      = identify_n_of_herbs_and_preds(sp_nm)
predator_has_prey     = check_predator_has_prey(sp_nm)

if !predator_has_prey[1]
    local_R -= predator_has_prey[2]
    filter!(name -> !(name in predator_has_prey[3]), sp_nm)
    @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    species_names = sp_nm
else
    species_names = sp_nm
end

localNPP              = Float64(npp_DA_relative_to_1000[local_i, local_j]) 
localH0_vector        = Vector{Float64}(H0_DA[local_i, local_j].a)

# ===========================
# 1️⃣ FITNESS FUNCTION (Survival Rate)
# ===========================
function fitness(params)
    mu, mu_predation, epsilon = params  # Unpack optimization parameters
    
    # Run the ecosystem model with these parameters
    results = attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true;
        localNPP      = localNPP,
        localH0_vector= localH0_vector,
        species_names = sp_nm,
        artificial_pi = art_pi
    )

    # If setup fails, return a bad fitness score
    if results === nothing
        return -0.01  # Small negative value to penalize failure
    end

    # Destructure results
    (S2, R2, species_names, herbivore_list, predator_list, H_i0, m_i, p_vec, x_final, g_i, localHatH, G, M_modified, a_matrix, A, epsilon_vector, m_alpha) = results

    # Solve ODE
    H_init = H_i0
    # println("H_init[1:R2]: ", H_init[1:R2])
    # println("Type of H_init[1:R2]: ", typeof(H_init[1:R2]))
    # println("Element types: ", [typeof(x) for x in H_init[1:R2]])

    P_init = H_init[1:R2] ./ 10.0
    u0     = vcat(H_init, P_init)
    
    params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
    prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)

    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end

    # Check if simulation failed
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return -0.01
    end

    # Evaluate survival rate
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species  # This is what we maximize

    return -survival_rate  # GA minimizes, so return negative
end

# ===========================
# 2️⃣ SETTING UP GA PARAMETERS
# ===========================

# Define parameter bounds
bounds = [
    (0.1, 0.9),  # mu range
    (0.0, 0.2),  # mu_predation range
    (0.0, 1.0)   # epsilon range
]

# Initialize an initial guess for parameters
initial_params = [rand(low:high) for (low, high) in bounds]

# Configure the Genetic Algorithm (GA)
ga_algorithm = GA(
    populationSize = 100,      # Number of individuals in the population
    selection = tournament(3), # Tournament selection with size 3
    mutationRate = 0.5,        # Mutation probability
    crossoverRate = 0.9        # Crossover probability
)

# Set optimization options (number of generations, etc.)
options = Evolutionary.Options(
    iterations = 1000,  # Number of generations
    show_trace = true  # Display progress
)

# ===========================
# 3️⃣ RUNNING THE OPTIMIZATION
# ===========================
result = Evolutionary.optimize(
    fitness,        # Function to minimize
    initial_params, # Initial parameters
    ga_algorithm,   # GA optimization settings
    options         # Options for max iterations and logging
)

# ===========================
# 4️⃣ DISPLAY BEST SOLUTION
# ===========================
best_params = Evolutionary.minimizer(result)
best_survival_rate = -Evolutionary.minimum(result)  # Convert back to positive

println("Best Parameters: mu=$(best_params[1]), mu_predation=$(best_params[2]), epsilon=$(best_params[3])")
println("Best Survival Rate: $best_survival_rate")
