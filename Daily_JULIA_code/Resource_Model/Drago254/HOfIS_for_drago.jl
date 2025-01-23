################################################################################
# Instead of storing the best configurations (max survivors) in a global DataFrame,
# we write each cell's best result (i.e. the line of best configuration) directly
# to a CSV file immediately after processing that cell.
################################################################################
include("prior.jl")
include("prior2.jl")
include("DA_birmmals_with_pi.jl")
include("generate_competition_matrix.jl")
include("species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
include("attempt_setup_community.jl")
include("Callbacks_function.jl")
include("npp_DA_relative_to_1000.jl")

# (If not already declared in your module, you might declare the following globals as const)
# const iberian_interact_NA = iberian_interact_NA
# const species_dict = species_dict
# const predator_names = predator_names
# const herbivore_names = herbivore_names
# const spain_names = spain_names

# Define your parameter ranges
mu_vals               = range(0.1, 0.9, length=10)
mu_predation_vals     = range(0.0, 0.5, length=130)
epsilon_vals          = range(0.1, 1.0, length=20)
sym_competition_vals  = [true]

# Optionally, for testing:
# mu_vals = [0.811111111]
# mu_predation_vals = [0.012121212]
# epsilon_vals = [0.8897959]
# sym_competition_vals = [true]

# 2) Build the parameter combinations and shuffle them 
param_combinations = [
    (mu, mu_predation, epsilon_val, sym_comp) 
    for mu in mu_vals
    for mu_predation in mu_predation_vals
    for epsilon_val in epsilon_vals
    for sym_comp in sym_competition_vals
]
param_combinations = Random.shuffle!(param_combinations)
param_combinations = param_combinations[1:end]

# Define constants
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 50000      # Up to 2000 combos
const SURVIVAL_THRESHOLD  = 0.0       # For example, store best if survival rate >= threshold

# Import SpinLock for thread-safety
using Base.Threads: SpinLock

# Initialize a SpinLock for file access
const file_lock = SpinLock()

function write_result_row(result::NamedTuple, filename::String)
    lock(file_lock)  # Lock the file for thread-safe access
    # Convert the result NamedTuple to a DataFrame with a single row.
    df = DataFrame([result])  # Wrap result in an array for a one-row DataFrame.
    try
        if !isfile(filename)
            # If the file does not exist, write with a header.
            CSV.write(filename, df)
        else
        # Otherwise, append without a header.
            CSV.write(filename, df; append=true)
        end
    finally
        unlock(file_lock)  # Unlock the file
    end
end

# Set the output file name (make sure the directory exists)
output_filename = "Results/best_params_5950_cells_not_rounded_254_26000iter_newherbivores.csv"

# Process cells (each cell will write its result immediately upon finding a good configuration)
# (Assuming that `idx` is defined and gives the (i,j) for each cell.)
@time Threads.@threads for cell in 1:5950  # or however many cells you want to process
    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]
    
    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey     = check_predator_has_prey(sp_nm)

    if !predator_has_prey[1]
       local_R -= predator_has_prey[2]
       filter!(name -> !(name in predator_has_prey[3]), sp_nm)
       @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
        species_names = sp_nm
    else
        species_names = sp_nm
    end
    
    localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j]) #1000.0
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)
    # cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)

    # We'll track the best result in local variables.
    best_survival_rate  = 0.0
    best_result         = nothing
    found_full_survival = false
    
    # 4) Up to 2000 combos (or fewer, if param_combinations < 2000)
    for (p_idx, combo) in enumerate(param_combinations[1:min(end,MAX_ITERS)])
        mu_val, mu_pred_val, eps_val, sym_competition = combo
        
        # Attempt setup
        results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_pred_val, eps_val, sym_competition;
            localNPP      = localNPP,
            localH0_vector= localH0_vector,
            species_names = species_names
        )
        if results === nothing
            continue
        end

        # Destructure the results:
        (S2, R2, H_i0, m_i, g_i, G, M_modified,
         a_matrix, A, epsilon_vector, m_alpha) = (
            results.S, results.R, results.H_i0, results.m_i,
            results.g_i, results.G, results.M_modified,
            results.a_matrix, results.A, results.epsilon_vector,
            results.m_alpha
        )

        if (S2 + R2) == 0 || R2 > length(H_i0)
            @error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
            continue
        end

        # Build initial conditions.
        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0     = vcat(H_init, P_init)

        params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
        prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
        
        #### THIS APPROACH ELIMINATES THE WARNING ####
        # Create a logger that only shows errors (i.e. skips warnings)
        logger = SimpleLogger(stderr, Logging.Error)

        # Wrap the solver call in the logging context:
        sol = with_logger(logger) do
            solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
        #######################################################

        # Skip if integration did not complete to 500 or if nan/inf occurred.
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
            continue
        end

        # Evaluate
        H_end          = sol[1:S2, end]
        P_end          = sol[S2+1:S2+R2, end]
        survived_herb  = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred  = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv     = survived_herb + survived_pred
        total_species  = S2 + R2
        survival_rate  = total_surv / total_species
        giHi           = sum(g_i .* H_end)
        ratio_ok       = (giHi / localNPP > 0.5) && (giHi / localNPP < 10.0)

        # If better than what's stored, update best result
        if survival_rate > best_survival_rate && ratio_ok
            # compute the rest:
            herbivore_survival_rate = survived_herb / S2
            predator_survival_rate  = survived_pred / R2
            H_biomass               = sum(H_end)
            P_biomass               = sum(P_end)
            biomass_at_the_end      = H_biomass + P_biomass
            ratio                   = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)

            best_survival_rate = survival_rate
            best_result = (
                cell_id                 = cell,
                i                       = local_i,
                j                       = local_j,
                mu                      = mu_val,
                mu_predation            = mu_pred_val,
                epsilon_val             = eps_val,
                symmetrical_competition = sym_competition,
                NPP                     = localNPP,
                g_iH_i                  = giHi,
                g_iH_i_over_NPP         = round(giHi / localNPP, digits=4), 
                survived_herbivores     = survived_herb,
                survived_predators      = survived_pred,
                total_survivors         = total_surv,
                total_species           = total_species,
                survival_rate           = survival_rate,
                herbivore_survival_rate = herbivore_survival_rate,
                predator_survival_rate  = predator_survival_rate,
                H_biomass               = H_biomass,
                P_biomass               = P_biomass,
                biomass_at_the_end      = biomass_at_the_end,
                herb_pred_ratio         = ratio
            )
        end

        # Early stop if full survival + the ratio check
        if isapprox(survival_rate, 1.0; atol=1e-10)
            ratio_ok = (giHi / localNPP > 0.5) && (giHi / localNPP < 10.0)
            if ratio_ok
                @info "Cell $cell => full survival with (mu=$mu_val, mu_pred=$mu_pred_val, eps=$eps_val). Stopping early."
                found_full_survival = true
                break
            end
        end
    end

    # Write the result immediately if conditions are met.
    if found_full_survival || (best_survival_rate > SURVIVAL_THRESHOLD)
         # Write best_result to CSV immediately.
         write_result_row(best_result, output_filename)
    else
         @info "Cell $cell: best survival was $(round(best_survival_rate, digits=2)) < $SURVIVAL_THRESHOLD => Not storing."
    end
end  # end Threads.@threads

@info "Finished. Stored combos with survival >= $SURVIVAL_THRESHOLD (or full survival)."
# error("Done")
