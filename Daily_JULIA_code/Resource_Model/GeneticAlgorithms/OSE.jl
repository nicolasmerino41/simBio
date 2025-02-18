using Evolutionary, Random, DifferentialEquations, Logging

# ===========================
# 1️⃣ GLOBAL PARAMETERS
# ===========================
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0
global const art_pi              = true

global cell = 1
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
# 1️⃣ FITNESS FUNCTION (Survival Rate)
# ===========================
function fitness(params)
    mu, mu_predation, epsilon = params  # Unpack parameters

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
        return -1.0  # Penalize failure
    end

    # Destructure results
    (S2, R2, _, _, _, H_i0, m_i, _, _, g_i, _, G, M_modified, a_matrix, A, epsilon_vector, m_alpha) = results

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
        return -1.0  # Penalize failed simulations
    end

    # Evaluate survival rate
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species

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
# 3️⃣ EVOLUTIONARY OPERATORS
# ===========================
function mutate(params, bounds, mutation_rate=0.2)
    mutated = copy(params)
    for i in eachindex(mutated)
        if rand() < mutation_rate
            low, high = bounds[i]
            mutated[i] = clamp(mutated[i] + randn() * 0.1, low, high)
        end
    end
    return mutated
end

function crossover(parent1, parent2)
    α = rand()
    return α * parent1 + (1 - α) * parent2
end

function evolve_population(population, bounds)
    new_population = []
    for _ in 1:length(population) ÷ 2
        p1, p2 = rand(population, 2)
        child1 = mutate(crossover(p1, p2), bounds)
        child2 = mutate(crossover(p2, p1), bounds)
        push!(new_population, child1, child2)
    end
    return new_population
end

# ===========================
# 4️⃣ OSE SEARCH ALGORITHM
# ===========================
function ose_search(; pop_size=50, generations=100, archive_size=200)
    
    # Define parameter bounds
    bounds = [(0.1, 1.0),  # mu range
              (0.1, 1.0),  # mu_predation range
              (0.1, 1.0)]  # epsilon range

    # Initialize population & archive
    population = [[rand(low:high) for (low, high) in bounds] for _ in 1:pop_size]
    archive = []

    for gen in 1:generations
        println("📌 Generation $gen")

        # Evaluate fitness & novelty
        fitness_scores = [fitness(p) for p in population]
        novelty_scores = [novelty_metric(p, archive) for p in population]

        # Store successful solutions in archive
        for (i, score) in enumerate(fitness_scores)
            if score == 1.0  # Only keep perfect survival cases
                push!(archive, population[i])
            end
        end
        if length(archive) > archive_size
            archive = archive[end-archive_size+1:end]  # Trim archive
        end
        # Select based on **both** novelty & fitness
        combined_scores = [(0.7 * fitness_scores[i] + 0.3 * novelty_scores[i]) for i in 1:pop_size]
        ranked_indices = sortperm(combined_scores, rev=true)
        selected = [population[i] for i in ranked_indices[1:pop_size]]

        # Apply genetic operators (mutation, crossover)
        population = evolve_population(selected, bounds)

    end

    return archive  # Return all diverse solutions with survival_rate = 1
end

# ===========================
# 5️⃣ RUN OSE EXPLORATION
# ===========================
best_solutions = ose_search(pop_size=10, generations=10, archive_size=50)

# ===========================
# 6️⃣ DISPLAY FINAL SOLUTIONS
# ===========================
println("✅ Found $(length(best_solutions)) viable parameter sets where survival_rate = 1")

for (i, params) in enumerate(best_solutions)
    println("🟢 Solution $i: mu=$(params[1]), mu_predation=$(params[2]), epsilon=$(params[3])")
end

single_run(1, 0.1, 0.1, 0.1172, true; artificial_pi = true)