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

# -------------------------------------------------------------------
# 0) Launching workers: Suppose we want 16 total processes
#    (one per node?), each with 48 threads:
#    *** If you only requested 1 process per node in your Slurm job,
#    change SlurmManager(16) to SlurmManager(1) etc. as needed. ***
# -------------------------------------------------------------------
addprocs(SlurmManager(16), exeflags=["--project", "--threads=48"])

# --------------------------
# 1) Bring all needed code/files to every worker
# --------------------------
@everywhere begin
    using CSV, DataFrames, DifferentialEquations, Distributions, Random, Logging
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
end

# --------------------------
# 2) Global constants & parameter combos
# --------------------------
@everywhere const EXTINCTION_THRESHOLD = 1e-6
@everywhere const T_ext               = 250.0
@everywhere const MAX_ITERS           = 50000
@everywhere const SURVIVAL_THRESHOLD  = 0.0

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
# const file_lock = SpinLock()

# --------------------------
# 4) Define a helper to append a NamedTuple row to CSV
# --------------------------
@everywhere function write_result_row(result::NamedTuple, filename::String)
    
    # Convert the result NamedTuple to a 1-row DataFrame
    df = DataFrame([result])
    if !isfile(filename)
        # If the file does not exist, write with header
        CSV.write(filename, df)
    else
        # Otherwise, append w/o header
        CSV.write(filename, df; append=true)
    end 
end

# --------------------------
# 5) process_cell: do the param search, return best NamedTuple or nothing
# --------------------------
@everywhere function process_cell(cell::Int, output_file::String)
    # Suppose idx is a global array: idx[k] = (i, j)
    local_i, local_j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell on worker $(myid()) @ node..."

    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

    # Optionally remove predators with no prey
    predator_pp = check_predator_has_prey(sp_nm)
    if !predator_pp[1]
        local_R -= predator_pp[2]
        filter!(name -> !(name in predator_pp[3]), sp_nm)
        @info "Cell $cell: removed $(predator_pp[2]) preds w/o prey."
    end

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # Build callbacks as needed
    cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)

    best_survival_rate  = 0.0
    best_result         = nothing
    found_full_survival = false

    for combo in param_combinations
        mu_val_local, mu_pred_val, eps_val, sym_comp = combo

        results = attempt_setup_community(
            local_i, local_j, mu_val_local, mu_pred_val, eps_val, sym_comp;
            localNPP=localNPP, localH0_vector=localH0_vector
        )
        if results === nothing
            continue
        end

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

        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0 = vcat(H_init, P_init)

        params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
        prob = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)

        logger = SimpleLogger(stderr, Logging.Error)
        sol = with_logger(logger) do
            solve(prob, Tsit5(); callback=cb_no_trigger, reltol=1e-6, abstol=1e-8)
        end

        # Skip if integration incomplete or if NaNs
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            continue
        end

        H_end = sol[1:S2, end]
        P_end = sol[S2+1:S2+R2, end]
        survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv    = survived_herb + survived_pred
        total_species = S2 + R2
        survival_rate_val = total_surv / total_species

        if survival_rate_val > best_survival_rate
            herbivore_survival_rate = (S2 == 0) ? 0.0 : (survived_herb / S2)
            predator_survival_rate  = (R2 == 0) ? 0.0 : (survived_pred / R2)
            H_biomass = sum(H_end)
            P_biomass = sum(P_end)
            biomass_at_the_end = H_biomass + P_biomass
            ratio = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
            giHi  = sum(g_i .* H_end)

            best_survival_rate = survival_rate_val
            best_result = (
                cell_id                 = cell,
                i                       = local_i,
                j                       = local_j,
                mu                      = mu_val_local,
                mu_predation            = mu_pred_val,
                epsilon_val             = eps_val,
                symmetrical_competition = sym_comp,
                NPP                     = localNPP,
                g_iH_i                  = giHi,
                survived_herbivores     = survived_herb,
                survived_predators      = survived_pred,
                total_survivors         = total_surv,
                total_species           = total_species,
                survival_rate           = survival_rate_val,
                herbivore_survival_rate = herbivore_survival_rate,
                predator_survival_rate  = predator_survival_rate,
                H_biomass               = H_biomass,
                P_biomass               = P_biomass,
                biomass_at_the_end      = biomass_at_the_end,
                herb_pred_ratio         = ratio
            )
        end

        # Early stop if full survival + the ratio check
        if isapprox(survival_rate_val, 1.0; atol=1e-10)
            ratio_ok = (giHi / localNPP > 0.5) && (giHi / localNPP < 10.0)
            if ratio_ok
                @info "Cell $cell: Full survival & ratio check => stopping early."
                found_full_survival = true
                break
            end
        end
    end

    # If there's a valid best_result, write it out now
    if found_full_survival || (best_survival_rate >= SURVIVAL_THRESHOLD && best_result !== nothing)
        write_result_row(best_result, output_file)
    end

    return nothing  # or return best_result if you want
end

# --------------------------
# 6) process_cell_range: run multiple cells in one chunk (one process),
#    using @threads for concurrency. Write results immediately.
# --------------------------
@everywhere function process_cell_range(cell_range::UnitRange{Int}, output_file::String)
    @threads for i in eachindex(cell_range)
        cell = cell_range[i]
        process_cell(cell, output_file)
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
chunks = chunk_indices(length(all_cells), 16)

output_filename = "Results/best_params_5950_cells_distributed.csv"

# Each process gets a chunk, then threads over cells in that chunk:
pmap(chunk -> process_cell_range(chunk, output_filename), chunks)

@info "All cells completed. Rows appended to $output_filename."
