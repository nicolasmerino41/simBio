using Evolutionary, DifferentialEquations, Random, Logging, NLsolve, ForwardDiff, DataFrames, Serialization
using Base.Threads: SpinLock, @threads

# ===========================
# 1️⃣ CONSTANTS & INIT
# ===========================
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000       # Up to 2000 generations
const SURVIVAL_THRESHOLD  = 0.0        # Only save configs with SR >= threshold
const art_pi              = false      # artificial_pi flag

const file_lock = SpinLock()           # Ensures only one thread writes to CSV/DF at a time

# (Set starting and ending indices for testing or HPC)
start_val = 1
end_val = 1

# ===========================
# 2️⃣ FITNESS FUNCTION
# ===========================
function abstract_fitness(params, S::Int, R::Int, localNPP::Float64)
    # Clamp parameters to bounds.
    params = clamp.(params, [0.1, 0.0, 0.0, 0.0], [0.9, 0.2, 1.0, 1.0])
    mu, mu_predation, epsilon, m_alpha = params
    
    # Generate abstract community parameters.
    results = abstract_parametrise_community(S, R;
                        NPP = localNPP, mu = mu, mu_pred = mu_predation,
                        epsilon = epsilon, c = 0.4, f = (R > 0 ? 0.3 : 0.0),
                        M_mean = m_alpha, alpha = 0.25)
    if isnothing(results)
        return 0.0
    end

    # Destructure parameters.
    S2       = results.S
    R2       = results.R
    H0_eff   = results.H0_eff
    m_i      = results.m
    g_i      = results.g
    beta     = results.beta
    comp_matrix = results.comp_matrix
    A_star   = results.A_star
    A_pred   = (R2 > 0 ? results.trophic_matrix : zeros(0, S2))
    A        = results.A
    m_alpha_vec = results.m_alpha
    h_val    = 0.1   # fixed handling time

    # Predator overpopulation thresholds.
    P0 = m_alpha_vec

    # Build parameter tuple (order must match abstract_equilibrium_system!).
    params_tuple = (S2, R2, H0_eff, m_i, g_i, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha_vec, h_val)
    x0 = vcat(H0_eff, fill(0.1, R2))
    
    sol = nlsolve((F, x) -> abstract_equilibrium_system!(F, x, params_tuple), x0)
    if !sol.f_converged
        @error "Equilibrium solver did not converge."
        return 0.0
    end

    H_eq = sol.zero[1:S2]
    P_eq = (R2 > 0 ? sol.zero[S2+1:end] : [])
    H_eq[H_eq .< 0.0] .= 0.0
    if R2 > 0
        P_eq[P_eq .< 0.0] .= 0.0
    end

    # Run ODE simulation.
    sim_results, sol_sim = abstract_herbivore_run(S, R, mu, mu_predation, epsilon, 0.1;
                                                  time_end = 500.0, include_predators = (R2 > 0),
                                                  NPP = localNPP, alpha = 0.25, hollingII = true, h = h_val,
                                                  do_you_want_sol = true)
    if sol_sim.t[end] < 500.0 || any(isnan, sol_sim.u[end]) || any(isinf, sol_sim.u[end])
        return 0.0
    end

    H_end = (R2 > 0 ? sol_sim.u[end][1:S2] : sol_sim.u[end])
    P_end = (R2 > 0 ? sol_sim.u[end][S2+1:end] : zeros(0))
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species

    return -survival_rate  # GA minimizes, so return negative SR
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
    show_trace = false
)

# Define parameter bounds.
cuts_for_mu = range(0.1, 0.9, length = 50)
cuts_for_mu_predation = range(0.0, 0.2, length = 50)
cuts_for_epsilon = range(0.0, 1.0, length = 50)
cuts_for_m_alpha = range(0.0, 1.0, length = 50)

csv_filepath = "abstract_ga_results.csv"
open(csv_filepath, "w") do file
    header = "cell_id,survival_rate,mu,mu_predation,epsilon_val,m_alpha,localNPP,S,R\n"
    write(file, header)
end

df = DataFrame(
    cell_id = Int[],
    survival_rate = Float64[],
    mu = Float64[],
    mu_predation = Float64[],
    epsilon_val = Float64[],
    m_alpha = Float64[],
    localNPP = Float64[],
    S = Int[],
    R = Int[]
)

# ===========================
# 4️⃣ MULTITHREADED GA OPTIMIZATION
# ===========================
@time Threads.@threads for portion in 1:8
    lower_bounds = [cuts_for_mu[portion], cuts_for_mu_predation[portion], cuts_for_epsilon[portion], cuts_for_m_alpha[portion]]
    upper_bounds = [cuts_for_mu[portion+1], cuts_for_mu_predation[portion+1], cuts_for_epsilon[portion+1], cuts_for_m_alpha[portion+1]]
    bounds = Evolutionary.BoxConstraints(lower_bounds, upper_bounds)

    # For the abstract framework, choose arbitrary values:
    cell = 1
    S_abstract = 5
    R_abstract = 2
    localNPP = 1000.0

    fitness_closure = params -> abstract_fitness(params, S_abstract, R_abstract, localNPP)
    result = Evolutionary.optimize(fitness_closure, bounds, ga_algorithm, options)
    best_params = Evolutionary.minimizer(result)
    best_survival_rate = -Evolutionary.minimum(result)

    # Compute equilibrium with these best parameters.
    results = abstract_parametrise_community(S_abstract, R_abstract;
                        NPP = localNPP, mu = best_params[1], mu_pred = best_params[2],
                        epsilon = best_params[3], c = 0.4, f = (R_abstract > 0 ? 0.3 : 0.0),
                        M_mean = 0.1, alpha = 0.25)
    S2 = results.S
    R2 = results.R
    H0_eff = results.H0_eff
    m_i = results.m
    g_i = results.g
    beta = results.beta
    comp_matrix = results.comp_matrix
    A_star = results.A_star
    A_pred = (R2 > 0 ? results.trophic_matrix : zeros(0, S2))
    A = results.A
    m_alpha_vec = results.m_alpha
    h_val = 0.1
    P0 = m_alpha_vec
    params_tuple = (S2, R2, H0_eff, m_i, g_i, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha_vec, h_val)
    x0 = vcat(H0_eff, fill(0.1, R2))
    sol = nlsolve((F, x) -> abstract_equilibrium_system!(F, x, params_tuple), x0)
    if sol.f_converged
        H_eq = sol.zero[1:S2]
        P_eq = (R2 > 0 ? sol.zero[S2+1:end] : [])
        H_eq[H_eq .< 0.0] .= 0.0
        if R2 > 0
            P_eq[P_eq .< 0.0] .= 0.0
        end
    else
        H_eq = nothing
        P_eq = nothing
    end

    # Use a lock to ensure thread-safe push! to the shared DataFrame.
    lock(file_lock) do
        push!(df, (
            cell,
            best_survival_rate,
            best_params[1],
            best_params[2],
            best_params[3],
            best_params[4],
            localNPP,
            S_abstract,
            R_abstract
        ))
    end

    # Write CSV row.
    H_eq_str = isnothing(H_eq) ? "None" : join(H_eq, ";")
    P_eq_str = isnothing(P_eq) ? "None" : join(P_eq, ";")
    row = string(cell) * "," *
          string(best_survival_rate) * "," *
          string(best_params[1]) * "," *
          string(best_params[2]) * "," *
          string(best_params[3]) * "," *
          string(best_params[4]) * "," *
          string(localNPP) * "," *
          string(S_abstract) * "," *
          string(R_abstract) * "\n"
    
    lock(file_lock) do
        open(csv_filepath, "a") do file
            write(file, row)
        end
    end

    @info "Finished abstract cell $cell: Best SR = $best_survival_rate"
end

serialize(df, "abstract_ga_results.jls")
