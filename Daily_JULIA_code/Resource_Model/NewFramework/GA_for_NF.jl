using Evolutionary, DifferentialEquations, Random, Logging

# Import SpinLock for thread-safety
using Base.Threads: SpinLock
using Base.Threads: @threads
# ===========================
# 2️⃣ PARSE ARGUMENTS & INIT
# ===========================
start_val = parse(Int, ARGS[1])  # Starting cell index (for HPC batch processing)
end_val   = parse(Int, ARGS[2])  # Ending cell index
start_val = 1
end_val = 8

const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000      # Up to 2000 generations
const SURVIVAL_THRESHOLD  = 0.0       # Only save configs with SR >= threshold
const global art_pi = false

const file_lock = SpinLock()  # Ensures only one thread writes to CSV at a time

# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function fitness(params, local_i, local_j, sp_nm, localNPP, localH0_vector)
    # Ensure parameters stay within bounds
    params = clamp.(params, [0.1, 0.0, 0.0], [0.9, 0.2, 1.0])
    mu, mu_predation, epsilon = params

    # Debug: Print current evaluation
    # println("🔍 Evaluating (cell $local_i, $local_j): mu=$(mu), mu_predation=$(mu_predation), epsilon=$(epsilon)")

    # Run the ecosystem model with these parameters
    # Attempt setup using the updated parametrisation
    results = new_attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true;
        localNPP       = localNPP,
        localH0_vector = localH0_vector,
        species_names  = sp_nm,
        artificial_pi  = art_pi
    )

    if isnothing(results)
        return 0.0  # Harsh penalty for failure
    end

    # Destructure the returned NamedTuple.
    S2          = results.S
    R2          = results.R
    H_i0        = results.H_i0
    m_i         = results.m_i
    g_i         = results.g_i
    x_final     = results.x_final   # scaling parameter (if needed)
    beta        = results.beta      # niche parameters
    G           = results.G
    M_modified  = results.M_modified
    a_matrix    = results.a_matrix
    A           = results.A
    epsilon = results.epsilon   # predator conversion vector
    m_alpha     = results.m_alpha
    
    # Solve ODE
    H_init = H_i0
    P_init = H_init[1:R2] ./ 10.0
    u0     = vcat(H_init, P_init)

    params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
    prob   = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params)

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

    return -survival_rate  # GA minimizes, so return negative SR
end

# ===========================
# 4️⃣ GA CONFIGURATION
# ===========================
lower_bounds = [0.1, 0.0, 0.0]
upper_bounds = [0.9, 0.5, 1.0]
bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

ga_algorithm = GA(
    populationSize = 100,
    selection = tournament(3),
    mutationRate = 0.25,
    crossoverRate = 0.6,
    epsilon = 0.1
)

options = Evolutionary.Options(
    iterations = 100,
    show_trace = true  # Suppress verbose output for HPC efficiency
)

@time Threads.@threads for cell in start_val:end_val

    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]

    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey = check_predator_has_prey(sp_nm)
    
    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    end

    localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # Pass thread-local parameters via a closure
    fitness_closure = params -> fitness(params, local_i, local_j, sp_nm, localNPP, localH0_vector)

    # Run GA Optimization for Current Cell
    result = Evolutionary.optimize(fitness_closure, bounds, ga_algorithm, options)

    # Extract best results
    best_params = Evolutionary.minimizer(result)
    best_survival_rate = -Evolutionary.minimum(result)  # Convert back to positive

    # Build CSV row
    row = string(cell) * "," * string(best_survival_rate) * "," * string(best_survival_rate == 1.0) * "," *
          string(best_params[1]) * "," * string(best_params[2]) * "," * string(best_params[3]) * "," *
          string(local_S + local_R) * "," * string(local_i) * "," * string(local_j) * "," * string(localNPP) * "\n"

    # Write row to CSV immediately
    lock(file_lock) do
        open("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF.csv", "a") do file
            write(file, row)
        end
    end

    println("✅ Finished Cell $cell: Best SR = $best_survival_rate\n")
end

println("\n📂 Results saved to 'ga_cell_results.csv'! 🚀")