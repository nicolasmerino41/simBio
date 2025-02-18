using Evolutionary, DifferentialEquations, Random, Logging, CSV, DataFrames, Base.Threads

# ===========================
# 1️⃣ PREPARE RESULTS & THREAD LOCK
# ===========================
results = DataFrame(
    cell = Int[],
    sr_max = Float64[],  # Best survival rate found
    reached_sr_1 = Bool[],  # Whether SR = 1.0 was reached
    mu = Float64[],
    mu_predation = Float64[],
    epsilon = Float64[],
    species_richness = Int[],
    local_i = Int[],
    local_j = Int[],
    localNPP = Float64[],
    localH0_vector = Vector{Float64}[]
)

num_cells = 8  # Change to 5950 later
results_lock = Threads.SpinLock()  # Ensures thread-safe access

# ===========================
# 2️⃣ RUN MULTI-THREADED GA OPTIMIZATION
# ===========================
Threads.@threads for cell in 1:num_cells
    # println("\n🟢 Processing cell $cell / $num_cells...\n")

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
        fitness(params, local_i, local_j, sp_nm, localNPP, localH0_vector)
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
    best_survival_rate = -Evolutionary.minimum(result)  # Convert back to positive

    # Accumulate results locally to prevent race conditions
    row = (
        cell,
        best_survival_rate,
        best_survival_rate == 1.0,  # Check if it reached SR = 1.0
        best_params[1],
        best_params[2],
        best_params[3],
        local_S+local_R,
        local_i,
        local_j,
        localNPP,
        localH0_vector
    )

    # Lock before pushing to shared `results` DataFrame
    lock(results_lock) do
        push!(results, row)
    end

    println("✅ Finished Cell $cell: Best SR = $best_survival_rate\n")
end

# ===========================
# 3️⃣ SAVE RESULTS SAFELY
# ===========================
CSV.write("ga_cell_results.csv", results)
println("\n📂 Results saved to 'ga_cell_results.csv'! 🚀")
