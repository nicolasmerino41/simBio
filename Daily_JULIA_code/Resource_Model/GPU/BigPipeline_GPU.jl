using CUDA, CSV, DataFrames, Random

# --------------------------
# 1) Load Required Files
# --------------------------
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

# --------------------------
# 2) Global Constants & Parameters
# --------------------------
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0
const art_pi              = true

# Parameter Combinations
param_combinations = CuArray(Random.shuffle!([
    (mu, mu_predation, epsilon_val, sym_comp)
    for mu in range(0.1, 0.9, length=3),
        mu_predation in range(0.0, 0.1, length=3),
        epsilon_val in range(0.1, 1.0, length=3),
        sym_comp in [true]
])[1:MAX_ITERS])

# --------------------------
# 3) GPU Kernel for Cell Processing
# --------------------------
function process_cell_kernel(cell_idx, output_gpu)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i > length(cell_idx)  # Ensure within bounds
        return
    end
    
    cell = cell_idx[i]

    # Extract species and set up parameters
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

    # Ensure predators have prey
    predator_has_prey = check_predator_has_prey(sp_nm)
    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
    end

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = CuArray(H0_DA[local_i, local_j].a)

    # Attempt feasibility
    feasible, best_result = attempt_feasibility(
        sp_nm, local_i, local_j, localNPP, localH0_vector, param_combinations;
        many_params = true, sp_removed = false, artificial_pi = art_pi
    )

    if feasible && !isnothing(best_result)
        output_gpu[i] = best_result.survival_rate  # Store result in GPU array
    else
        output_gpu[i] = 0.0  # Default value if infeasible
    end
end

# --------------------------
# 4) Run GPU Kernel
# --------------------------
function process_cells_on_gpu(all_cells::Vector{Int})
    cell_idx_gpu = CuArray(all_cells)  # Move cell indices to GPU
    output_gpu = CUDA.zeros(length(all_cells))  # Create result array on GPU

    threads_per_block = 256
    blocks_per_grid = ceil(Int, length(all_cells) / threads_per_block)

    @cuda threads=threads_per_block blocks=blocks_per_grid process_cell_kernel(cell_idx_gpu, output_gpu)

    return Array(output_gpu)  # Move results back to CPU
end

# --------------------------
# 5) Run and Save Results
# --------------------------
all_cells = collect(1:5950)
results = process_cells_on_gpu(all_cells)

# Convert to DataFrame and save results
df = DataFrame(cell=all_cells, survival_rate=results)
CSV.write("Results/GPU_results.csv", df)

@info "All cells processed using GPU. Results saved to 'GPU_results.csv'."
