################################################################################
# Instead of storing the best configurations (max survivors) in a global DataFrame,
# we write each cell's best result (i.e. the line of best configuration) directly
# to a CSV file immediately after processing that cell.
################################################################################
include("Scripts/prior.jl")
include("Scripts/prior2.jl")
include("Scripts/DA_birmmals_with_pi.jl")
include("Scripts/generate_competition_matrix.jl")
include("Scripts/species_dict.jl")

include("Scripts/ecosystem_dynamics!.jl")
include("Scripts/FI_functions.jl")
include("Scripts/extract_H0_DA.jl")
include("Scripts/attempt_setup_community.jl")
include("Scripts/Callbacks_function.jl")
include("Scripts/npp_DA_relative_to_1000.jl")
include("Scripts/attempt_feasibility.jl")

# Import SpinLock for thread-safety
using Base.Threads: SpinLock
using Base.Threads: @threads

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
mu_predation_vals     = range(0.0, 0.2, length=130)
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
global param_combinations = param_combinations[1:end]

# Define constants
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000      # Up to 2000 combos
const SURVIVAL_THRESHOLD  = 0.0       # For example, store best if survival rate >= threshold

# Set the output file names (ensure directories exist)
output_filename = "Results/Big_pipeline_results_drago_fixed.csv"
problematic_species_filename = "Results/problematic_species_drago_fixed.csv"
slightly_problematic_species_filename = "Results/slightly_problematic_species_drago_fixed.csv"

Threads.@threads for cell in 1:5950 # Adjust the range as needed
    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]
    
    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
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

    # Step 1: Attempt full feasibility
    feasible, best_result = attempt_feasibility(species_names, local_i, local_j, localNPP, localH0_vector, param_combinations; many_params = true, sp_removed = false)

    if feasible
        # Write the best_result to CSV
        if !isnothing(best_result)
            @info "Cell $cell is feasible initially."
            write_result_row(best_result, output_filename)
        end
    else
        if !isnothing(best_result)
            
            parameters = [(best_result.mu, best_result.mu_predation, best_result.epsilon_val, best_result.symmetrical_competition)]
            # println("length(parameters) = ", length(parameters))
            key = false
            @info "Cell $cell did have at least one non-nothing result..."
            
        else
            key = true
            parameters = param_combinations
            @info "Cell $cell did not have any non-nothing results..."
        end

        @info "Cell $cell is not feasible initially. Attempting species removal..."
        # Step 2: Iteratively remove species to achieve feasibility
        problematic_species_list = String[]
	    slightly_problematic_species_list = String[]
        best_survival_rate_from_sp_removed = isnothing(best_result) ? 0.0 : best_result.survival_rate
        for species in sp_nm
            modified_sp_nm = filter(s -> s != species, sp_nm)
            feasible, best_result = attempt_feasibility(modified_sp_nm, local_i, local_j, localNPP, localH0_vector, parameters; many_params = key, sp_removed = true, sp_removed_name = species)
            if feasible
                push!(problematic_species_list, species)
                @info "Removing species '$species' makes cell $cell feasible."
                write_result_row(best_result, output_filename)

            elseif !feasible && !isnothing(best_result) && best_result.survival_rate > best_survival_rate_from_sp_removed
		        push!(slightly_problematic_species_list, species)
                write_result_row(best_result, output_filename)
                best_survival_rate_from_sp_removed = best_result.survival_rate
            end
        end

        # If feasibility achieved after removal, record the problematic species
        if !isempty(problematic_species_list)
            write_problematic_species(cell, problematic_species_list, problematic_species_filename)
            # Optionally, write the best_result to CSV if desired
        elseif !isempty(slightly_problematic_species_list)
            write_problematic_species(cell, slightly_problematic_species_list, slightly_problematic_species_filename)
            # Optionally, write the best_result to CSV if desired
        else
            # If still not feasible, flag for further investigation
            @info "Cell $cell remains infeasible after attempting species removal."
            # Optionally, log this cell for manual review
        end
    end
end

@info "Finished Step 1: Identified problematic species for cells where full survival wasn't initially feasible."
