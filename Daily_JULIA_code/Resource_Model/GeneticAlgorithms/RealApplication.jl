using Evolutionary, DifferentialEquations, Random, Logging

# ===========================
# 1️⃣ GLOBAL PARAMETERS
# ===========================
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0
const art_pi              = false

global cell = 2
global local_i, local_j = idx[cell][1], idx[cell][2]

# Gather cell data
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
global local_S, local_R      = identify_n_of_herbs_and_preds(sp_nm)
global predator_has_prey     = check_predator_has_prey(sp_nm)

if !predator_has_prey[1]
    local_R -= predator_has_prey[2]
    filter!(name -> !(name in predator_has_prey[3]), sp_nm)
    @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    global species_names = sp_nm
else
    global species_names = sp_nm
end

global localNPP              = Float64(npp_DA_relative_to_1000[local_i, local_j]) 
global localH0_vector        = Vector{Float64}(H0_DA[local_i, local_j].a)

# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function fitness(params)
    # Ensure parameters stay within bounds (in case of mutation errors)
    params = clamp.(params, [0.1, 0.0, 0.0], [0.9, 0.2, 1.0])
    mu, mu_predation, epsilon = params

    # Debug: Print current evaluation
    # println("🔍 Evaluating: mu=$(mu), mu_predation=$(mu_predation), epsilon=$(epsilon)")

    # Run the ecosystem model with these parameters
    results = attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true;
        localNPP = localNPP,
        localH0_vector = localH0_vector,
        species_names = sp_nm,
        artificial_pi = art_pi
    )

    if results === nothing
        return 0.0  # Harsh penalty for failure
    end

    # Destructure results
    (
        S2, R2, species_names, herbivore_list, predator_list,
        H_i0, m_i, p_vec, x_final, g_i,
        localHatH, G, M_modified, a_matrix, A, epsilon_vector, m_alpha
    ) = results

    # Solve ODE
    H_init = H_i0
    P_init = H_init[1:R2] ./ 10.0
    u0     = vcat(H_init, P_init)

    params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
    prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)

    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end

    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return 0.0  # Harsh penalty for failed simulations
    end

    # Evaluate survival rate
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species  # Maximize this

    # Debug: Print high survival rate cases
    if survival_rate = 1.0
        println("🎯 High SR found! SR = $(survival_rate)")
    end

    return -survival_rate  # GA minimizes, so return negative SR
end

# ===========================
# 3️⃣ SPECIFY PARAMETER BOUNDS (BoxConstraints)
# ===========================
lower_bounds = [0.1, 0.0, 0.0]  # Lower bounds for mu, mu_predation, epsilon
upper_bounds = [0.9, 0.2, 1.0]  # Upper bounds

bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

# ===========================
# 4️⃣ CONFIGURE GENETIC ALGORITHM
# ===========================
# Configure the Genetic Algorithm (GA)
ga_algorithm = GA(
    populationSize = 300,      # Large population size for better diversity
    selection = tournament(3), # Higher tournament size for better selection
    mutationRate = 0.2,        # 🔥 Drastically increase mutation rate
    crossoverRate = 0.5,       # Reduce crossover rate to favor mutations
    ε = 0.1                    # Allow some randomness in selection
)

# ===========================
# 5️⃣ RUN OPTIMIZATION
# ===========================
options = Evolutionary.Options(
    iterations = 100,  # Number of generations
    show_trace = true  # Display progress
)

result = Evolutionary.optimize(
    fitness,        # Fitness function
    bounds,         # Box constraints (no need for manual initial params)
    ga_algorithm,   # GA optimization settings
    options         # Options for iterations and logging
)

# ===========================
# 6️⃣ DISPLAY BEST SOLUTION
# ===========================
best_params = Evolutionary.minimizer(result)
best_survival_rate = -Evolutionary.minimum(result)  # Convert back to positive

println("🎯 Best Parameters: mu=$(best_params[1]), mu_predation=$(best_params[2]), epsilon=$(best_params[3])")
println("🔥 Best Survival Rate: $best_survival_rate")

# Print top 5 best solutions GA found
println("🔥 GA Best Solutions Found:")
for (i, sol) in enumerate(result.history[1:min(5, length(result.history))])
    println("Solution $i: params=", sol.minimizer, " → SR=", -sol.minimum)
end
