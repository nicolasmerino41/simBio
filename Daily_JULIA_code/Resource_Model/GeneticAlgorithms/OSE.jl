using Evolutionary, Random, DifferentialEquations, Logging

# ===========================
# 1️⃣ FITNESS FUNCTION (Survival Rate)
# ===========================
function fitness(params)
    mu, mu_predation, epsilon = params  # Unpack parameters

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
        return -1.0  # Penalize failure
    end

    # Destructure results
    (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha) = results

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

    # Check if simulation failed
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return -1.0
    end

    # Evaluate survival rate
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species  # This is what we maximize

    return survival_rate  # Directly maximize survival rate
end

# ===========================
# 2️⃣ NOVELTY FUNCTION (Encouraging Exploration)
# ===========================
function novelty_metric(candidate, archive)
    if isempty(archive)
        return 0.0  # If archive is empty, no novelty yet
    end
    distances = [norm(candidate .- archived) for archived in archive]
    return mean(sort(distances)[1:min(5, length(distances))])  # Distance from 5 nearest neighbors
end

# ===========================
# 3️⃣ OSE SEARCH ALGORITHM
# ===========================
function ose_search(; pop_size=50, generations=100, archive_size=200)
    
    # Define parameter bounds
    bounds = [(0.1, 1.0),  # mu range
              (0.1, 1.0),  # mu_predation range
              (0.1, 1.0)]  # epsilon range

    # Initialize population and archive
    population = [rand(low:high) for (low, high) in bounds for _ in 1:pop_size]
    archive = []

    for gen in 1:generations
        println("Generation $gen")

        # Evaluate fitness & novelty
        fitness_scores = [fitness(p) for p in population]
        novelty_scores = [novelty_metric(p, archive) for p in population]

        # Store successful solutions in archive
        for (i, score) in enumerate(fitness_scores)
            if score == 1.0  # Only keep perfect survival cases
                push!(archive, population[i])
        end
        if length(archive) > archive_size
            archive = archive[end-archive_size+1:end]  # Trim archive

        # Select based on **both** novelty and fitness
        combined_scores = [(0.7 * fitness_scores[i] + 0.3 * novelty_scores[i]) for i in 1:pop_size]
        selected = population[sortperm(combined_scores, rev=true)][1:pop_size]

        # Apply genetic operators (mutation, crossover)
        new_population = evolve(selected)

        population = new_population
    end

    return archive  # Return all diverse solutions with survival_rate = 1
end

# ===========================
# 4️⃣ RUN OSE EXPLORATION
# ===========================
best_solutions = ose_search(pop_size=100, generations=200, archive_size=300)

# ===========================
# 5️⃣ DISPLAY FINAL SOLUTIONS
# ===========================
println("Found $(length(best_solutions)) viable parameter sets where survival_rate = 1")

for (i, params) in enumerate(best_solutions)
    println("Solution $i: mu=$(params[1]), mu_predation=$(params[2]), epsilon=$(params[3])")
end