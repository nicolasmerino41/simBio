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
end_val = 1

const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000      # Up to 2000 generations
const SURVIVAL_THRESHOLD  = 0.0       # Only save configs with SR >= threshold
const global art_pi = false

const file_lock = SpinLock()  # Ensures only one thread writes to CSV at a time

# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function fitness(params, local_i, local_j, sp_nm, localNPP)
    # Ensure parameters stay within bounds
    params = clamp.(params, [0.1, 0.0, 0.0, 0.0], [0.9, 0.2, 1.0, 1.0])
    mu, mu_predation, epsilon, m_alpha = params
    # mu, mu_predation, epsilon, m_alpha, alpha = 0.0, 0.0, 0.0, 0.1, 0.25
    # Debug: Print current evaluation
    # println("🔍 Evaluating (cell $local_i, $local_j): mu=$(mu), mu_predation=$(mu_predation), epsilon=$(epsilon)")

    # Run the ecosystem model with these parameters
    # Attempt setup using the updated parametrisation
    results = new_attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true, m_alpha;
        localNPP       = localNPP,
        # localH0_vector = localH0_vector,
        species_names  = sp_nm,
        artificial_pi  = art_pi,
        alpha          = 0.25,
        hollingII      = true, h = 0.1
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
    beta        = results.beta
    M_mod       = results.M_modified    # Modified competition matrix (S×S)
    A_star      = results.A_star        # Nondimensional predation rates (S×R)
    a_matrix    = results.a_matrix      # Herbivore–predator interaction matrix (S×R)
    A           = results.A             # Predator interaction matrix (R×R)
    m_alpha     = results.m_alpha
    h           = results.h
    
    # Build predator attack matrix (R×S) by transposing a_matrix.
    A_pred = transpose(a_matrix)

    # Set predator overpopulation thresholds.
    # (Here we assume P0 = m_alpha, which implies a self-regulation coefficient of 1.)
    P0 = m_alpha
    params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha, h)
    # SOLVING EQUILIBRIUM DENSITIES
    S = S2
    R = R2
    H0 = H_i0
    # Use H0 as the initial guess for herbivores, and a small constant (e.g., 0.1) for predators.
    x0 = vcat(H0, fill(0.1, R))
    # Now call the root finder.
    sol = nlsolve((F, x) -> equilibrium_system!(F, x, params), x0)

    if sol.f_converged && any(sol.zero[S+1:end] .< 0.0)
        H_eq = sol.zero[1:S]
        P_eq = sol.zero[S+1:end]
    elseif sol.f_converged && !any(sol.zero[S+1:end] .< 0.0)
        H_eq = sol.zero[1:S]
        P_eq = sol.zero[S+1:end]
    else
        @error "Equilibrium solver did not converge."
        return 0.0
    end
    # Solve ODE
    # H_init = H_i0 # THIS MAKES MORE SENSE BUT THE FOLLOWING IS A STANDARD
    H_eq[H_eq .< 0.0] .= 0.0
    H_init = H_eq
    # P_init = H_init[1:R2] ./ 10.0 # THIS MAKES MORE SENSE BUT THE FOLLOWING IS A STANDARD
    P_eq[P_eq .< 0.0] .= 0.0
    P_init = P_eq
    u0     = vcat(H_init, P_init)

    params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha, h)
    prob = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params)

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
    H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
    P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
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
# lower_bounds = [0.1, 0.0, 0.0, 0.0]
# upper_bounds = [0.9, 0.5, 1.0, 1.0]
# bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

ga_algorithm = GA(
    populationSize = 10,
    selection = tournament(3),
    mutationRate = 0.25,
    crossoverRate = 0.6,
    epsilon = 0.1
)

options = Evolutionary.Options(
    iterations = 100,
    show_trace = false  # Suppress verbose output for HPC efficiency
)

cuts_for_mu = range(0.0, 0.99, length = 9)
cuts_for_mu_predation = range(0.0, 1.0, length = 9)
cuts_for_epsilon = range(0.0, 1.0, length = 9)
cuts_for_m_alpha = range(0.0, 1.0, length = 9)
# cuts_for_alpha = range(0.0, 1.0, length = 9)

csv_filepath = "Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF.csv"
open(csv_filepath, "w") do file
    header = "cell_id,survival_rate,flag,mu,mu_predation,epsilon_val,m_alpha,total_species,i,j,NPP\n"
    write(file, header)
end
function compute_equilibrium(
    mu, mu_predation, epsilon, m_alpha,
    local_i, local_j, sp_nm, localNPP;
    alpha = 0.25,
    hollingII = true,
    h = 0.1
)
    # 1) Set up community parameters
    results = new_attempt_setup_community(
        local_i, local_j,
        mu, mu_predation, epsilon, true, m_alpha;
        localNPP = localNPP,
        species_names = sp_nm,
        artificial_pi = art_pi,
        alpha = alpha,
        hollingII = hollingII, h = h
    )
    if isnothing(results)
        return nothing, nothing
    end
    
    # 2) Destructure
    S2          = results.S
    R2          = results.R
    H_i0        = results.H_i0
    m_i         = results.m_i
    g_i         = results.g_i
    beta        = results.beta
    M_mod       = results.M_modified
    A_star      = results.A_star
    a_matrix    = results.a_matrix
    A           = results.A
    m_alpha_vec = results.m_alpha
    h_val       = results.h
    
    # 3) Build predator matrix and define P0
    A_pred = transpose(a_matrix)
    P0 = m_alpha_vec
    
    # 4) Set up param tuple
    params_tuple = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha_vec, h_val)
    
    # 5) Solve equilibrium system using nlsolve
    S = S2
    R = R2
    x0 = vcat(H_i0, fill(0.1, R))
    sol = nlsolve((F, x) -> equilibrium_system!(F, x, params_tuple), x0)
    
    if sol.f_converged
        # Extract solutions
        H_eq = copy(sol.zero[1:S])
        P_eq = copy(sol.zero[S+1:S+R])
        # zero out any negative solutions
        H_eq[H_eq .< 0.0] .= 0.0
        P_eq[P_eq .< 0.0] .= 0.0
        return H_eq, P_eq
    else
        return nothing, nothing
    end
end

# Initialize a DataFrame to store GA + equilibrium results.
df = DataFrame(
    cell_id = Int[],
    survival_rate = Float64[],
    mu = Float64[],
    mu_predation = Float64[],
    epsilon_val = Float64[],
    m_alpha = Float64[],
    i = Int[],
    j = Int[],
    localNPP = Float64[],
    H_eq = Vector{Vector{Float64}}[],  # store as a vector-of-vectors
    P_eq = Vector{Vector{Float64}}[]
)


@time Threads.@threads for portion in 1:8
    lower_bounds = [cuts_for_mu[portion], cuts_for_mu_predation[portion], cuts_for_epsilon[portion], cuts_for_m_alpha[portion]]
    upper_bounds = [cuts_for_mu[portion+1], cuts_for_mu_predation[portion+1], cuts_for_epsilon[portion+1], cuts_for_m_alpha[portion+1]]
    bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

    cell = 1    
    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]

    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey = check_predator_has_prey(sp_nm)
    
    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    end

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])

    fitness_closure = params -> fitness(params, local_i, local_j, sp_nm, localNPP)

    result = Evolutionary.optimize(fitness_closure, bounds, ga_algorithm, options)

    best_params = Evolutionary.minimizer(result)
    best_survival_rate = -Evolutionary.minimum(result)

    mu_opt, mu_pred_opt, eps_opt, m_alpha_opt = best_params

    # 1) Compute equilibrium with these best parameters
    H_eq, P_eq = compute_equilibrium(
        mu_opt, mu_pred_opt, eps_opt, m_alpha_opt,
        local_i, local_j, sp_nm, localNPP;
        alpha = 0.25, hollingII = true, h = 0.1
    )

    # 2) Append to DataFrame
    push!(df, (
        cell_id         = cell,
        survival_rate   = best_survival_rate,
        mu              = mu_opt,
        mu_predation    = mu_pred_opt,
        epsilon_val     = eps_opt,
        m_alpha         = m_alpha_opt,
        i               = local_i,
        j               = local_j,
        localNPP        = localNPP,
        H_eq            = isnothing(H_eq) ? Vector{Float64}() : [H_eq],
        P_eq            = isnothing(P_eq) ? Vector{Float64}() : [P_eq]
    ))

    # 3) Also write a CSV row if needed
    H_eq_str = isnothing(H_eq) ? "None" : join(H_eq, ";")
    P_eq_str = isnothing(P_eq) ? "None" : join(P_eq, ";")
    row = string(cell) * "," *
          string(best_survival_rate) * "," *
          string(best_survival_rate == 1.0) * "," *
          string(mu_opt) * "," *
          string(mu_pred_opt) * "," *
          string(eps_opt) * "," *
          string(m_alpha_opt) * "," *
          string(local_S + local_R) * "," *
          string(local_i) * "," *
          string(local_j) * "," *
          string(localNPP) * "," *
          H_eq_str * "," *
          P_eq_str * "\n"

    lock(file_lock) do
        open(csv_filepath, "a") do file
            write(file, row)
        end
    end

    println("✅ Finished Cell $cell: Best SR = $best_survival_rate")
end

serialize(df, "Results/3-3/ga_results_NF_hollingII.jls")

row = df[8, :]
mu_val = row.mu
mu_pred_val = row.mu_predation
eps_val = row.epsilon_val
m_alpha_val = row.m_alpha
H_init = copy(row.H_eq[1])
P_init = copy(row.P_eq[1])

sim_result = herbivore_run(
    row.cell_id, mu_val, mu_pred_val, eps_val, true, m_alpha_val;
    include_predators = true,
    plot = true,
    H_init = H_init,
    P_init = P_init,
    alpha = 0.25,
    hollingII = true,
    h = 0.1
)
