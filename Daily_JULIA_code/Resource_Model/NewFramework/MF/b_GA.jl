using Evolutionary, DifferentialEquations, Random, Logging
using Distributions
using Base.Threads: SpinLock, @threads

# ===========================
# 1️⃣ GLOBAL CONSTANTS & SETTINGS
# ===========================
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000      # Up to 2000 generations
const SURVIVAL_THRESHOLD  = 0.0       # Only save configs with SR >= threshold
const art_pi = false

# Global parameters for fitness (can be tuned):
const connectance = 1.0
const d_i_global = 1.0         # Can be a scalar or a vector
const d_alpha_global = 1.0     # For predators; note: lower d_alpha means weaker self-regulation → higher K_alpha.
const delta_nu_global = 0.05

# (Assume DA_birmmals_with_pi_corrected, npp_DA_relative_to_1000, local_i, local_j, cell are defined globally for the cell.)

# File lock for thread-safety
const file_lock = SpinLock()

# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function fitness(params, S, R)
    # Ensure parameters are within bounds (3 parameters: mu, epsilon, m_alpha).
    params = clamp.(params, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    mu, epsilon, m_alpha = params

    # Generate random herbivore names and predator names.
    herbivore_names = ["H$(i)" for i in 1:S]
    predator_names  = ["P$(i)" for i in 1:R]
    species_names = vcat(herbivore_names, predator_names)
    
    # Build species dictionary (herbivores first, then predators).
    species_dict = Dict{String, Int}()
    for (i, sp) in enumerate(species_names)
        species_dict[sp] = i
    end
    
    Random.seed!(123)
    # Generate random herbivore abundances (observed data).
    # Here using uniform between 5.0 and 10.0.
    H_i0 = [rand(Uniform(5.0, 10.0)) for _ in 1:S]
       
    # Generate a random iberian_interact matrix.
    # Dimensions: (S+R) x (S+R). For predator–herbivore interactions,
    # for rows S+1:S+R and columns 1:S, set entry = 1 with probability connectance.
    N = S + R
    iberian_interact = zeros(Float64, N, N)
    for j in 1:R  # predators are in rows S+1 to S+R
        global_pred_idx = S + j
        for i in 1:S  # herbivores are in columns 1:S
            if rand() < connectance
                iberian_interact[global_pred_idx, i] = 1.0
            end
        end
    end
    
    # Call the parameterisation function using our global settings.
    results = random_parametrise_the_community(
            species_names;
            mu = mu,
            epsilon_val = epsilon,
            mean_m_alpha = m_alpha,
            d_i = d_i_global,
            d_alpha = d_alpha_global,
            iberian_interact = iberian_interact,
            species_dict = species_dict,
            cell_abundance = H_i0,
            delta_nu = delta_nu_global
        )

    # Destructure returned results.
    S = results.S
    R = results.R
    H_star = results.H_star   # equilibrium herbivore abundances
    P_star = results.P_star   # computed predator equilibria
    
    # Build ODE parameter tuple.
    # Order: (S, R, K_i, r_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha)
    params_ode = (results.S, results.R, results.K_i, results.r_i, results.mu, results.nu,
                  results.P_matrix, results.epsilon, results.m_alpha, results.K_alpha)
    
    # Set initial conditions as the computed equilibria.
    u0 = vcat(H_star, P_star)
    prob = ODEProblem(bipartite_dynamics!, u0, (0.0, 500.0), params_ode)
    
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
    end

    # Penalize if simulation fails.
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return 0.0
    end

    # Extract final state from solution.
    H_end = sol[1:S, end]
    P_end = sol[S+1:S+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    survival_rate = total_surv / total_species  # not used in fitness but for debugging
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    pred_herb_ratio = P_biomass / H_biomass

    # We want to maximize pred_herb_ratio, but GA minimizes, so return negative.
    return -pred_herb_ratio
end

function fitness_combined(params, S, R; w_surv=0.5, w_pred=0.5)
    # Clamp parameters (3 parameters: mu, epsilon, m_alpha)
    params = clamp.(params, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    mu, epsilon, m_alpha = params

    # Generate species names.
    herbivore_names = ["H$(i)" for i in 1:S]
    predator_names  = ["P$(i)" for i in 1:R]
    species_names = vcat(herbivore_names, predator_names)
    
    # Build species dictionary (herbivores first, then predators).
    species_dict = Dict{String, Int}()
    for (i, sp) in enumerate(species_names)
        species_dict[sp] = i
    end
    
    Random.seed!(123)
    # Generate random herbivore abundances (observed data), e.g. Uniform between 5 and 10.
    H_i0 = [rand(Uniform(5.0, 10.0)) for _ in 1:S]
       
    # Generate a random interaction matrix.
    N = S + R
    iberian_interact = zeros(Float64, N, N)
    for j in 1:R  # predators in rows S+1:S+R
        global_pred_idx = S + j
        for i in 1:S  # herbivores in columns 1:S
            if rand() < connectance
                iberian_interact[global_pred_idx, i] = 1.0
            end
        end
    end

    # Call our parameterisation function (using global settings for d_i, d_alpha, and delta_nu).
    results = random_parametrise_the_community(
        species_names;
        mu = mu,
        epsilon_val = epsilon,
        mean_m_alpha = m_alpha,
        d_i = d_i_global,
        d_alpha = d_alpha_global,
        iberian_interact = iberian_interact,
        species_dict = species_dict,
        cell_abundance = H_i0,
        delta_nu = delta_nu_global
    )
    
    # Extract returned values.
    S = results.S
    R = results.R
    H_star = results.H_star   # equilibrium herbivore abundances (proxy from H_i0)
    P_star = results.P_star   # computed predator equilibria

    # Build ODE parameters tuple (order must match bipartite_dynamics!).
    ode_params = (results.S, results.R, results.K_i, results.r_i, results.mu, results.nu,
                  results.P_matrix, results.epsilon, results.m_alpha, results.K_alpha)
    
    # Set initial conditions to the computed equilibrium values.
    u0 = vcat(H_star, P_star)
    prob = ODEProblem(bipartite_dynamics!, u0, (0.0, 500.0), ode_params)
    
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
    end

    # Penalize if simulation fails.
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return 0.0
    end

    # Extract final state from the solution.
    H_end = sol[1:S, end]
    P_end = sol[S+1:S+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    survival_rate = total_surv / total_species
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    pred_herb_ratio = iszero(H_biomass) ? NaN : P_biomass / H_biomass

    # Combine objectives (we wish to maximize both survival_rate and pred_herb_ratio).
    # Since GA minimizes, we return the negative weighted sum.
    if iszero(pred_herb_ratio)
        combined_obj = 0.0
    else
        combined_obj = w_surv * survival_rate + w_pred * pred_herb_ratio
    end
    return -combined_obj
end

# ===========================
# 3️⃣ GA CONFIGURATION
# ===========================
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

# Divide parameter space among threads.
cuts_for_mu = range(0.0, 1.0, length = Threads.nthreads() + 1)
cuts_for_epsilon = range(0.0, 1.0, length = Threads.nthreads() + 1)
cuts_for_m_alpha = range(0.0, 1.0, length = Threads.nthreads() + 1)

# Initialize a DataFrame to store GA + equilibrium results.
A_df = DataFrame(
    pred_herb_ratio = Float64[],
    mu = Float64[],
    epsilon_val = Float64[],
    m_alpha = Float64[]
)

S, R = 2, 1

cb_no_trigger, cb_trigger = build_callbacks(S, R, EXTINCTION_THRESHOLD, T_ext, 1)

@time Threads.@threads for portion in 1:Threads.nthreads()
    lower_bounds = [cuts_for_mu[portion], cuts_for_epsilon[portion], cuts_for_m_alpha[portion]]
    upper_bounds = [cuts_for_mu[portion+1], cuts_for_epsilon[portion+1], cuts_for_m_alpha[portion+1]]
    bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)
    
    fitness_closure = params -> fitness(params, S, R)

    result = Evolutionary.optimize(fitness_closure, bounds, ga_algorithm, options)

    best_params = Evolutionary.minimizer(result)
    best_pred_herb_ratio = -Evolutionary.minimum(result)

    mu_opt, eps_opt, m_alpha_opt = best_params

    # Append to DataFrame (wrap vector outputs so each column is a single entry).
    @lock file_lock begin
        push!(A_df, (
            pred_herb_ratio = best_pred_herb_ratio,
            mu              = mau_opt,
            epsilon_val     = eps_opt,
            m_alpha         = m_alpha_opt
        ))
    end

    println("✅ Finished Cell: Best pred/herb ratio = $best_pred_herb_ratio")
end

# serialize(A_df, "Results/3-3/random_GA_NF_hollingII.jls")
