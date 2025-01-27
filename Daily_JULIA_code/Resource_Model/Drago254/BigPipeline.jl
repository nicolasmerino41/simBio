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


# Import SpinLock for thread-safety
using Base.Threads: SpinLock

# Initialize a SpinLock for file access
const file_lock = SpinLock()

# Define a helper function to write one result row to CSV.
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

# Define a helper function to write problematic species to a separate CSV.
function write_problematic_species(cell::Int, species::Vector{String}, filename::String)
    df = DataFrame(cell_id = [cell], problematic_species = [join(species, ";")])
    if !isfile(filename)
        CSV.write(filename, df)
    else
        CSV.write(filename, df; append=true, writeheader=false)
    end
end

# Define your parameter ranges
mu_vals               = range(0.1, 0.9, length=10)
mu_predation_vals     = range(0.0, 0.5, length=130)
epsilon_vals          = range(0.1, 1.0, length=20)
sym_competition_vals  = [true]

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
const MAX_ITERS           = 2000      # Up to 2000 combos
const SURVIVAL_THRESHOLD  = 0.0       # For example, store best if survival rate >= threshold

# Set the output file names (ensure directories exist)
output_filename = "Results/Big_pipeline_results.csv"
problematic_species_filename = "Results/problematic_species.csv"

@threads for cell in 1:8  # Adjust the range as needed
    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]
    
    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    
    # Check and remove predators without any prey.
    predator_has_prey = check_predator_has_prey(sp_nm)
    # Assuming check_predator_has_prey returns a tuple like (all_have_prey, num_without, names_without)
    
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

    # Track the best result and problematic species
    best_survival_rate  = 0.0
    best_result         = nothing
    found_full_survival = false
    problematic_species = String[]  # To store species that, when removed, enable feasibility
    
    # Function to attempt feasibility with current species list
    function attempt_feasibility(current_sp_nm)
        # Update local_S and local_R based on current species
        current_S, current_R = identify_n_of_herbs_and_preds(current_sp_nm)
        # Iterate over parameter combinations
        for combo in param_combinations[1:min(end, MAX_ITERS)]
            mu_val, mu_pred_val, eps_val, sym_competition = combo

            # Attempt to set up the community
            results = attempt_setup_community(
                local_i, local_j,
                mu_val, mu_pred_val, eps_val, sym_competition;
                localNPP       = localNPP,
                localH0_vector = localH0_vector,
                species_names  = current_sp_nm
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
                continue
            end

            # Build initial conditions.
            H_init = H_i0
            P_init = H_init[1:R2] ./ 10.0
            u0     = vcat(H_init, P_init)

            params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
            prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
            
            #### THIS APPROACH ELIMINATES THE SOLVER WARNING ####
            logger = SimpleLogger(stderr, Logging.Error)
            sol = with_logger(logger) do
                solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
            end
            #######################################################

            # Skip if integration did not complete to 500 or if nan/inf occurred.
            if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
                continue
            end

            # Evaluate survival
            H_end      = sol[1:S2, end]
            P_end      = sol[S2+1:S2+R2, end]
            survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
            survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
            total_surv    = survived_herb + survived_pred
            total_species = S2 + R2
            survival_rate = total_surv / total_species

            # Update the best result if survival_rate is improved.
            if survival_rate > best_survival_rate
                herbivore_survival_rate = survived_herb / S2
                predator_survival_rate  = survived_pred / R2
                H_biomass               = sum(H_end)
                P_biomass               = sum(P_end)
                biomass_at_the_end      = H_biomass + P_biomass
                ratio                   = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
                giHi                    = sum(g_i .* H_end)

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

            # Stop early if full survival is achieved.
            if isapprox(survival_rate, 1.0; atol=1e-10)
                @info "Cell $cell => full survival with (mu=$mu_val, mu_pred=$mu_pred_val, eps=$eps_val). Stopping early."
                return true  # Feasibility achieved
            end
        end
        # @info "No feasible combination found for cell $cell."
        return false  # Feasibility not achieved
    end

    # Step 1: Attempt full feasibility
    feasible = attempt_feasibility(sp_nm)

    if feasible
        # Write the best_result to CSV
        if best_result !== nothing
            write_result_row(best_result, output_filename)
        end
    else
        # Step 2: Iteratively remove species to achieve feasibility
        problematic_species_list = String[]
        for species in sp_nm
            modified_sp_nm = filter(s -> s != species, sp_nm)
            if attempt_feasibility(modified_sp_nm)
                push!(problematic_species_list, species)
                @info "Removing species '$species' makes cell $cell feasible."
                break  # Optionally, remove multiple species iteratively
            end
        end

        # If feasibility achieved after removal, record the problematic species
        if !isempty(problematic_species_list)
            write_problematic_species(cell, problematic_species_list, problematic_species_filename)
            # Optionally, write the best_result to CSV if desired
            if best_result !== nothing && best_survival_rate >= SURVIVAL_THRESHOLD
                write_result_row(best_result, output_filename)
            end
        else
            # If still not feasible, flag for further investigation
            @info "Cell $cell remains infeasible after attempting species removal."
            # Optionally, log this cell for manual review
        end
    end
end

@info "Finished Step 1: Identified problematic species for cells where full survival wasn't initially feasible."
