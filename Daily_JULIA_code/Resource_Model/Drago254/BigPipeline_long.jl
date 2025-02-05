################################################################################
# Distributed Cell Processing with Serial Parameter Search
#
# We distribute cells across multiple processes via pmap (1 chunk per process),
# then within each process, we use threads to handle cells in parallel.
# For each cell, we do our parameter search (serial) and immediately
# append the best result row to a CSV. We use a SpinLock for thread-safety
# within a single process, but *not* across processes.
################################################################################

using Distributed
using ClusterManagers
using CSV, DataFrames, Logging, Random
using Base.Threads: SpinLock, @threads
num_nodes = 16
num_threads = 48
queue = "long"
# -------------------------------------------------------------------
# 0) Launching workers: Suppose we want 16 total processes
#    (one per node?), each with 48 threads:
#    *** If you only requested 1 process per node in your Slurm job,
#    change SlurmManager(16) to SlurmManager(1) etc. as needed. ***
# -------------------------------------------------------------------
addprocs(SlurmManager(num_nodes), exeflags=["--project", "--threads=$num_threads"])

# --------------------------
# 1) Bring all needed code/files to every worker
# --------------------------
@everywhere begin
    using CSV, DataFrames, DifferentialEquations, Distributions, Random, Logging 
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
end

# --------------------------
# 2) Global constants & parameter combos
# --------------------------
@everywhere const EXTINCTION_THRESHOLD = 1e-6
@everywhere const T_ext               = 250.0
@everywhere const MAX_ITERS           = 2000
@everywhere const SURVIVAL_THRESHOLD  = 0.0
@everywhere const art_pi              = true

# We build 50,000 combos total (10*100*50*1=50k):
@everywhere global param_combinations = Random.shuffle!([
    (mu, mu_predation, epsilon_val, sym_comp)
    for mu in range(0.1, 0.9, length=10),
        mu_predation in range(0.0, 0.1, length=200),
        epsilon_val in range(0.1, 1.0, length=50),
        sym_comp in [true]
])[1:MAX_ITERS]

# --------------------------
# 3) A global spinlock for safe *threaded* appends
#    Note: This does not lock across multiple processes!
# --------------------------
@everyconst file_lock = SpinLock()

# --------------------------
# 4) Define a helper to append a NamedTuple row to CSV
# --------------------------
@everywhere function write_result_row(result::NamedTuple, filename::String)
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
@everywhere function write_problematic_species(cell::Int, species::Vector{String}, filename::String)
    df = DataFrame(cell_id = [cell], problematic_species = [join(species, ";")])
    if !isfile(filename)
        CSV.write(filename, df)
    else
        CSV.write(filename, df; append=true, writeheader=false)
    end
end

# --------------------------
# 5) process_cell: do the param search, return best NamedTuple or nothing
# --------------------------
@everywhere function process_cell(
    cell::Int, output_filename::String,
    problematic_species_filename::String,
    slightly_problematic_species_filename::String
)

    # Suppose idx is a global array: idx[k] = (i, j)
    local_i, local_j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell on worker $(myid()) @ node..."

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

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # Build callbacks as needed
    # cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)
    
    # Step 1: Attempt full feasibility
    feasible, best_result = attempt_feasibility(
        species_names, local_i, local_j,
        localNPP, localH0_vector, param_combinations; 
        many_params = true, sp_removed = false,
        artificial_pi = art_pi
    )
    
    if feasible
        # Write the best_result to CSV
        if !isnothing(best_result)
            @info "Cell $cell is feasible initially."
            write_result_row(best_result, output_filename)
        end
    else
        if !isnothing(best_result)
            parameters = param_combinations
            write_result_row(best_result, output_filename)
            key = true
            @info "Cell $cell did have at least one non-nothing result but we don't care"
            
        else
            key = true
            parameters = param_combinations
            @info "Cell $cell did not have any non-nothing results but we don't care"
        end
        
        @info "Cell $cell is not feasible initially. Attempting species removal..."
        # Step 2: Iteratively remove species to achieve feasibility
        problematic_species_list = String[]
	    slightly_problematic_species_list = String[]
        best_survival_rate_from_sp_removed = isnothing(best_result) ? 0.0 : best_result.survival_rate
        for species in sp_nm
            modified_sp_nm = filter(s -> s != species, sp_nm)
            feasible, best_result = attempt_feasibility(
                modified_sp_nm, local_i, local_j,
                localNPP, localH0_vector, parameters;
                many_params = key,
                sp_removed = true, sp_removed_name = species,
                artificial_pi = art_pi
                )
            if feasible
                push!(problematic_species_list, species)
                @info "Removing species '$species' makes cell $cell feasible."
                write_result_row(best_result, output_filename)

            elseif !feasible && best_result.survival_rate > best_survival_rate_from_sp_removed
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

# --------------------------
# 6) process_cell_range: run multiple cells in one chunk (one process),
#    using @threads for concurrency. Write results immediately.
# --------------------------
@everywhere function process_cell_range(
    cell_range::UnitRange{Int}, output_filename::String,
    problematic_species_filename::String,
    slightly_problematic_species_filename::String
)
    @threads for i in eachindex(cell_range)
        cell = cell_range[i]
        process_cell(cell, output_filename, problematic_species_filename, slightly_problematic_species_filename)
    end
    # We don't return anything, because we've already appended
    return nothing
end

# --------------------------
# 7) Utility to chunk the 5950 cells into 16 blocks (for 16 processes)
# --------------------------
function chunk_indices(total_cells::Int, n_chunks::Int)
    chunk_size = ceil(Int, total_cells / n_chunks)
    chunks = Vector{UnitRange{Int}}(undef, n_chunks)
    start_idx = 1
    for i in 1:n_chunks
        end_idx = min(start_idx + chunk_size - 1, total_cells)
        chunks[i] = start_idx:end_idx
        start_idx = end_idx + 1
        if start_idx > total_cells
            return filter(x -> !isempty(x), chunks)
        end
    end
    return filter(x -> !isempty(x), chunks)
end

# --------------------------
# 8) Main driver: pmap over chunks, each chunk calls `process_cell_range`
# --------------------------
all_cells = 1:5950
chunks = chunk_indices(length(all_cells), num_nodes)

# Set the output file names (ensure directories exist)
negation = art_pi ? "even_pi" : "not_even_pi"
output_filename = "Results/Big_pipeline_results_($queue)_$negation.csv"
problematic_species_filename = "Results/problematic_species_($queue)_$negation.csv"
slightly_problematic_species_filename = "Results/slightly_problematic_species_($queue)_$negation.csv"

# Each process gets a chunk, then threads over cells in that chunk:
pmap(chunk -> process_cell_range(chunk, output_filename, problematic_species_filename, slightly_problematic_species_filename), chunks)

@info "All cells completed. Rows appended to $output_filename."
