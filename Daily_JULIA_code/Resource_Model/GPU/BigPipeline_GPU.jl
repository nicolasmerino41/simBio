using CUDA, CSV, DataFrames, Random

# --------------------------
# 1) Load Required Files (CPU Only)
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

# Convert parameter combinations into separate arrays
param_combinations_cpu = [
    (mu, mu_predation, epsilon_val, sym_comp)
    for mu in range(0.1, 0.9, length=3),
        mu_predation in range(0.0, 0.1, length=3),
        epsilon_val in range(0.1, 1.0, length=3),
        sym_comp in [true]
]

# Flatten into individual CuArrays
mu_vals = CuArray([p[1] for p in param_combinations_cpu])
mu_predation_vals = CuArray([p[2] for p in param_combinations_cpu])
epsilon_vals = CuArray([p[3] for p in param_combinations_cpu])
sym_comp_vals = CuArray([p[4] for p in param_combinations_cpu])

# --------------------------
# 3) Precompute CPU Values (Move Expensive Computations Out of GPU)
# --------------------------
function preprocess_cell(cell::Int)
    local_i, local_j = idx[cell][1], idx[cell][2]
    
    # Run CPU functions *before* sending data to GPU
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

    predator_has_prey = check_predator_has_prey(sp_nm)
    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
    end

    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = H0_DA[local_i, local_j].a  # Keep on CPU

    return (local_i, local_j, localNPP, localH0_vector, sp_nm)
end

# Preprocess all cells before GPU execution
function preprocess_all_cells(all_cells::Vector{Int})
    return map(preprocess_cell, all_cells)
end

# --------------------------
# 4) GPU Kernel for Processing Feasibility
# --------------------------
function process_cell_kernel(cell_idx, output_gpu, mu_vals, mu_predation_vals, epsilon_vals, sym_comp_vals)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i > length(cell_idx)  # Ensure within bounds
        return
    end

    # Compute feasibility (Simplified for GPU Execution)
    feasible = (mu_vals[i] + mu_predation_vals[i] + epsilon_vals[i]) > 0.5  # Placeholder Logic

    if feasible
        output_gpu[i] = mu_vals[i] * 10  # Store survival rate in GPU array (Example Computation)
    else
        output_gpu[i] = 0.0  # Default if infeasible
    end

    return
end

# --------------------------
# 5) Run GPU Kernel
# --------------------------
function process_cells_on_gpu(all_cells::Vector{Int})
    cell_idx_gpu = CuArray(all_cells)  # Move cell indices to GPU
    output_gpu = CUDA.zeros(Float32, length(all_cells))  # Allocate result array on GPU

    threads_per_block = 256
    blocks_per_grid = ceil(Int, length(all_cells) / threads_per_block)

    @cuda threads=threads_per_block blocks=blocks_per_grid process_cell_kernel(
        cell_idx_gpu, output_gpu, mu_vals, mu_predation_vals, epsilon_vals, sym_comp_vals
    )

    return Array(output_gpu)  # Move results back to CPU
end

# --------------------------
# 6) Run and Save Results
# --------------------------
all_cells = collect(1:5)
preprocessed_data = preprocess_all_cells(all_cells)  # Run CPU Preprocessing
results = process_cells_on_gpu(all_cells)

# Convert to DataFrame and save results
df = DataFrame(cell=all_cells, survival_rate=results)
CSV.write("Results/GPU_results.csv", df)

@info "All cells processed using GPU. Results saved to 'GPU_results.csv'."
