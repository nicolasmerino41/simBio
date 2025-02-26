using Evolutionary, DifferentialEquations, Random, Logging, CSV, DataFrames, Base.Threads
const global art_pi = true
# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function fitness(params, local_i, local_j, sp_nm, localH0_vector)
    # Ensure parameters stay within bounds
    params = clamp.(params, [0.1, 0.0, 0.0, 0.0], [0.9, 0.2, 1.0, 1000.0])
    mu, mu_predation, epsilon, localNPP = params

    # Debug: Print current evaluation
    # println("🔍 Evaluating (cell $local_i, $local_j): mu=$(mu), mu_predation=$(mu_predation), epsilon=$(epsilon)")

    # Run the ecosystem model with these parameters
    results = attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true;
        localNPP = localNPP,
        localH0_vector = localH0_vector,
        species_names = sp_nm,
        artificial_pi = true
    )

    if results === nothing
        return 10.0  # Harsh penalty for failure
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
        return 10.0  # Harsh penalty for failed simulations
    end

    # Evaluate survival rate
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species  # Maximize this
    biomass_end = vcat(H_end, P_end)
    # Checking empirical correlation
    real_abundances = filter(x -> !iszero(x), DA_birmmals_with_pi_corrected[local_i, local_j].a)
    cor_real_pred = cor(real_abundances, biomass_end)

    return -cor_real_pred  # GA minimizes, so return negative SR
end

# ===========================
# 3️⃣ SPECIFY PARAMETER BOUNDS (BoxConstraints)
# ===========================
lower_bounds = [0.1, 0.0, 0.0, 0.0]  # Lower bounds for mu, mu_predation, epsilon, NPP
upper_bounds = [0.9, 0.5, 1.0, 1000.0]  # Upper bounds

bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

# ===========================
# 4️⃣ CONFIGURE GENETIC ALGORITHM
# ===========================
# Configure the Genetic Algorithm (GA)
ga_algorithm = GA(
    populationSize = 10,      # Large population size for better diversity
    selection = tournament(3), # Higher tournament size for better selection
    mutationRate = 0.2,        # 🔥 Drastically increase mutation rate
    crossoverRate = 0.5,       # Reduce crossover rate to favor mutations
    ε = 0.1                    # Allow some randomness in selection
)

# ===========================
# 5️⃣ RUN OPTIMIZATION
# ===========================
options = Evolutionary.Options(
    iterations = 1000,  # Number of generations
    show_trace = true  # Display progress
)

# result = Evolutionary.optimize(
#     fitness,        # Fitness function
#     bounds,         # Box constraints (no need for manual initial params)
#     ga_algorithm,   # GA optimization settings
#     options         # Options for iterations and logging
# )

# ===========================
# 1️⃣ PREPARE RESULTS & THREAD LOCK
# ===========================
results = DataFrame(
    cell = Int[],
    best_similarity = Float64[],
    reached_sim_1 = Bool[],  # Whether sim = 1.0 was reached
    # survival_rate = Float64[],
    mu = Float64[],
    mu_predation = Float64[],
    epsilon = Float64[],
    localNPP = Float64[],
    species_richness = Int[],
    local_i = Int[],
    local_j = Int[],
    localH0_vector = Vector{Float64}[],
    # H_vector = Vector{Float64}[],
    # P_vector = Vector{Float64}[]
)

num_cells = 2  # Change to 5950 later
results_lock = Threads.SpinLock()  # Ensures thread-safe access

# ===========================
# 2️⃣ RUN MULTI-THREADED GA OPTIMIZATION
# ===========================
Threads.@threads for cell in 1:1
    println("\n🟢 Processing cell $cell / $num_cells...\n")

    # Local variables to avoid global conflicts
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey = check_predator_has_prey(sp_nm)

    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        species_names = sp_nm
    else
        species_names = sp_nm
    end

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j]) 
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # ===========================
    # Pass thread-local parameters via a closure
    # ===========================
    fitness_closure = function(params)
        fitness(params, local_i, local_j, sp_nm, localH0_vector)
    end

    # ===========================
    # Run GA Optimization for Current Cell
    # ===========================
    result = Evolutionary.optimize(
        fitness_closure,  # Fitness function (now thread-local)
        bounds,           # Box constraints
        ga_algorithm,     # GA settings
        options           # Optimization options
    )

    # Extract best results
    best_params = Evolutionary.minimizer(result)
    best_similarity = -Evolutionary.minimum(result)  # Convert back to positive

    # model = single_run(
    #     cell, best_params[1], best_params[2], best_params[3], true;
    #     plot=false, sp_removed_name=nothing,
    #     artificial_pi=art_pi, NPP=best_params[4],
    #     herbivore_m = 0.1,
    #     predator_m = 0.1
    # )

    # Accumulate results locally to prevent race conditions
    row = (
        cell,
        best_similarity,
        best_similarity == 1.0,  # Check if it reached SR = 1.0,
        # model.survival_rate[1],
        best_params[1],
        best_params[2],
        best_params[3],
        best_params[4],
        local_S+local_R,
        local_i,
        local_j,
        # localNPP,
        localH0_vector,
        # model.H_vector[1],
        # model.P_vector[1]
    )

    # Lock before pushing to shared `results` DataFrame
    lock(results_lock) do
        push!(results, row)
    end

    println("✅ Finished Cell $cell: Best sim = $best_similarity\n")  #, SR = $(model.survival_rate[1])\n")
end

# ===========================
# 3️⃣ CHECK THE OUTPUT
# ===========================
cell = 1
model = single_run(
    cell, results[cell, :mu], results[cell, :mu_predation], results[cell, :epsilon], true;
    plot=true, sp_removed_name=nothing,
    artificial_pi=true, NPP=results[cell, :localNPP],
    herbivore_m = 0.1,
    predator_m = 0.1
)

real_abundances = filter(x -> !iszero(x), DA_birmmals_with_pi_corrected[local_i, local_j].a)
cor_real_pred = cor(real_abundances, biomass_end)

# ===========================
# 3️⃣ SAVE RESULTS SAFELY
# ===========================
CSV.write("ga_cell_results.csv", results)
println("\n📂 Results saved to 'ga_cell_results.csv'! 🚀")
