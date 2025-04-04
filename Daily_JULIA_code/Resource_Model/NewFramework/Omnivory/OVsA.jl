# -----------------------------------------------------------------------------
# Assumptions:
# - analytical_equilibrium(cell, ...) is defined elsewhere and returns a NamedTuple with:
#     • equilibrium: a NamedTuple containing H_star (herbivore/omnivore equilibrium abundances),
#                    P_star (predator equilibrium abundances),
#                    herbivore_list (a vector of species names in the herbivore compartment),
#                    predator_list (a vector of predator names).
# - parameters: a NamedTuple containing interaction matrices and the computed numbers S and R.
#
# - omnivore_names and herbivore_names are defined (for example, as global sets or arrays)
#   so that one can check membership (e.g.:
#       omnivore_names = Set(["Omnivore 1", "Omnivore 2", ...])
#       herbivore_names = Set(["Herbivore 1", "Herbivore 2", ...])
# )
# -----------------------------------------------------------------------------

# --- 1. Structural Metrics Computation for Real Communities ---
function compute_structural_metrics_real(comm)
    # Extract the herbivore list from the equilibrium.
    all_herb = comm.equilibrium.herbivore_list
    # Identify pure herbivores (not in omnivore_names) and omnivores.
    pure_H = [sp for sp in all_herb if !(sp in omnivore_names)]
    actual_O = [sp for sp in all_herb if sp in omnivore_names]
    # Define proportion omnivores as the fraction of omnivores among all herbivore compartment species.
    prop_omn = isempty(all_herb) ? 0.0 : length(actual_O) / length(all_herb)
    # Total species = herbivores+omnivores (in all_herb) plus predators.
    total_species = length(all_herb) + comm.parameters.R

    # Retrieve interaction matrices.
    S = comm.parameters.S         # Number of herbivore compartment species (herbivores + omnivores)
    R = comm.parameters.R         # Number of predators
    P_matrix = comm.parameters.P_matrix
    O_matrix = comm.parameters.O_matrix
    T_matrix = comm.parameters.T_matrix
    B_matrix = comm.parameters.B_matrix
    D_matrix = comm.parameters.D_matrix

    # Compute connectance for each species.
    conn_vals = zeros(S+R)
    # For each herbivore (or omnivore) species (indices 1:S):
    for i in 1:S
        d = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
        # Maximum possible interactions: attacked by all predators and interactions with all other herbivores (both dicircleions).
        max_possible = R + 2 * (S - 1)
        conn_vals[i] = max_possible > 0 ? d / max_possible : 0.0
    end
    # For each predator (indices S+1:S+R):
    for alpha in 1:R
        d = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
        # Maximum possible interactions: can feed on all herbivores and interact with all other predators.
        max_possible = S + 2 * (R - 1)
        conn_vals[S+alpha] = max_possible > 0 ? d / max_possible : 0.0
    end
    avg_conn = mean(conn_vals)

    # Compute the degree (total number of interactions) for each species.
    n = S + R
    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
    end
    avg_degree = mean(degree)
    
    return (prop_omn = prop_omn, total_species = total_species, avg_conn = avg_conn, avg_degree = avg_degree)
end

# --- 2. Stable Configuration Search for Real Communities ---
"""
    find_stable_configurations_real(cell; mu_range, eps_range, m_alpha_range, time_end)

# Loops over ranges for μ, ε, and m_alpha for the given cell (real community),
using analytical_equilibrium to generate the community. Returns all configurations
with a locally stable equilibrium (i.e. all eigenvalues of the Jacobian have negative real parts).
"""
function find_stable_configurations_real(cell::Int; 
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0,
    m_alpha_range=0.05:0.05:0.2, time_end=500.0,
    plot = false
)
    
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = analytical_equilibrium(
                    cell, 
                    mu_val, eps_val, m_alpha_val;
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = true, pi_size = 10.0,
                    H_init = nothing,
                    P_init = nothing,
                    nu_omni_proportion = 1.0,
                    nu_b_proportion = 1.0,
                    r_omni_proportion = 1.0,
                    callbacks = false, plot = plot
                )
                if comm === nothing
                    continue
                end
                u0 = vcat(comm.equilibrium.H_star, comm.equilibrium.P_star)
                params = (comm.parameters.S, comm.parameters.R, comm.parameters.K_i, comm.parameters.r_i,
                          comm.parameters.mu, comm.parameters.nu, comm.parameters.P_matrix, comm.parameters.O_matrix,
                          comm.parameters.T_matrix, comm.parameters.epsilon, comm.parameters.m_alpha,
                          comm.parameters.K_alpha, comm.parameters.B_matrix, comm.parameters.D_matrix,
                          comm.parameters.nu_omni, comm.parameters.nu_b)
                sol, J = simulate_community(u0, params, time_end)                          
                # sol, J = simulate_community(u0, params, time_end)
                # Check stability: all eigenvalues must have negative real part.
                if all(real.(eigen(J).values) .< 0)
                    # println("Cell $cell: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val)")
                    push!(stable_configs, (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val))
                end
            end
        end
    end
    if isempty(stable_configs)
        error("No stable configuration found for cell $cell")
    else
        return stable_configs
    end
end

# --- 3. Metrics Computation for Real Communities ---
function compute_resilience(J)
    # Resilience approximated as the spectral abscissa (maximum real part of eigenvalues)
    return maximum(real.(eigen(J).values))
end

function compute_reactivity(J)
    S_sym = (J + transpose(J)) / 2
    return maximum(real.(eigen(S_sym).values))
end

function compute_persistence(sol; threshold=1e-6)
    final_state = sol.u[end]
    n = length(final_state)
    survivors = count(x -> x > threshold, final_state)
    return survivors / n
end

# --- 4. Experiment Pipeline for Real Communities ---
"""
    run_experiment_real(cell_range; conn, mu_range, eps_range, m_alpha_range, time_end)

Loops over each cell in cell_range (each representing a real community).
For each cell, uses find_stable_configurations_real to get at least one stable configuration,
computes stability metrics (resilience, reactivity, persistence) and structural metrics
(proportion omnivores, total species, average connectance), and stores the results.
"""
function run_experiment_real(
    cell_range=1:10;
    time_end=500.0,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:0.2,
    plot = false
)
    
    results = DataFrame(cell=Int[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[],
                        prop_omn=Float64[], total_species=Int[], avg_conn=Float64[],
                        avg_degree=Float64[])
    
    lock_obj = ReentrantLock()
    
    Threads.@threads for cell in cell_range
        try
            stable_configs = find_stable_configurations_real(
                cell;
                mu_range=mu_range, eps_range=eps_range,
                m_alpha_range=m_alpha_range, time_end=time_end,
                plot = plot
                )
            # Average metrics over all stable configurations found:
            all_resilience = mean([compute_resilience(stable.J) for stable in stable_configs])
            all_reactivity = mean([compute_reactivity(stable.J) for stable in stable_configs])
            all_persistence = mean([compute_persistence(stable.sol) for stable in stable_configs])
            println("Cell $cell: Resilience=$(all_resilience), Reactivity=$(all_reactivity), Persistence=$(all_persistence)")
            # Compute structural metrics from the first stable configuration.
            struct_metrics = compute_structural_metrics_real(stable_configs[1].comm)
            
            # Thread-safe push!
            lock(lock_obj) do
                push!(results, (
                    cell,
                    all_resilience, all_reactivity, all_persistence,
                    struct_metrics.prop_omn, struct_metrics.total_species, struct_metrics.avg_conn,
                    struct_metrics.avg_degree
                ))
            end
        catch e
            println("Skipping cell $cell: ", e)
        end
    end
    return results
end

# --- 5. Plotting Results for Real Communities ---
function plot_experiment_results_real(results::DataFrame)
    fig1 = Figure(resolution=(800,600))
    ax1 = Axis(fig1[1,1], xlabel="Average Connectance", ylabel="Resilience", title="Resilience vs. Connectance")
    scatter!(ax1, results.avg_conn, results.resilience, markersize=8, color=:blue)
    
    ax2 = Axis(fig1[1,2], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
    scatter!(ax2, results.prop_omn, results.reactivity, markersize=8, color=:green)
    
    ax3 = Axis(fig1[2,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Total Species")
    scatter!(ax3, results.total_species, results.persistence, markersize=8, color=:red)

    ax4 = Axis(fig1[2,2], xlabel="Average Degree", ylabel="Resilience", title="Resilience vs. Average Degree")
    scatter!(ax4, results.avg_degree, results.resilience, markersize=8, color=:blue)
    
    display(fig1)
end

# --- 6. Running the Real-Community Experiment Pipeline ---
# A_real_results = run_experiment_real(
#     1:1000;
#     time_end=500.0,
#     mu_range=0.0:0.25:1.0, eps_range=0.0:0.25:1.0, m_alpha_range=0.05:0.05:0.2,
#     plot = false
# )
# println(real_results)
plot_experiment_results_real(A_real_results)

####################################################################################
####################################################################################
############ THE REAL COMPARISON BETWEEN ABSTRACT AND REAL COMMUNITIES #############
####################################################################################
####################################################################################
# =============================================================================
# 1. Structural Metrics Computation
# =============================================================================
# For real communities:
function compute_structural_metrics_real(comm)
    all_herb = comm.equilibrium.herbivore_list
    # Identify omnivores (assume omnivore_names is defined globally)
    actual_O = [sp for sp in all_herb if sp in omnivore_names]
    prop_omn = isempty(all_herb) ? 0.0 : length(actual_O) / length(all_herb)
    total_species = length(all_herb) + comm.parameters.R

    S = comm.parameters.S   # herbivore compartment size (herbivores+omnivores)
    R = comm.parameters.R
    P_matrix = comm.parameters.P_matrix
    O_matrix = comm.parameters.O_matrix
    T_matrix = comm.parameters.T_matrix
    B_matrix = comm.parameters.B_matrix
    D_matrix = comm.parameters.D_matrix

    conn_vals = zeros(S+R)
    for i in 1:S
        d = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
        max_possible = R + 2 * (S - 1)
        conn_vals[i] = max_possible > 0 ? d / max_possible : 0.0
    end
    for alpha in 1:R
        d = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
        max_possible = S + 2 * (R - 1)
        conn_vals[S+alpha] = max_possible > 0 ? d / max_possible : 0.0
    end
    avg_conn = mean(conn_vals)

    n = S + R
    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
    end
    avg_degree = mean(degree)
    return (prop_omn = prop_omn, total_species = total_species, avg_conn = avg_conn, avg_degree = avg_degree)
end

# For abstract communities:
function compute_structural_metrics_abstract(comm)
    all_herb = comm.herbivore_list
    # Assume abstract omnivores have names starting with "Omnivore"
    actual_O = [sp for sp in all_herb if occursin("Omnivore", sp)]
    prop_omn = isempty(all_herb) ? 0.0 : length(actual_O) / length(all_herb)
    total_species = length(all_herb) + length(comm.predator_list)

    S = comm.S
    R = comm.R
    P_matrix = comm.P_matrix
    O_matrix = comm.O_matrix
    T_matrix = comm.T_matrix
    B_matrix = comm.B_matrix
    D_matrix = comm.D_matrix

    conn_vals = zeros(S+R)
    for i in 1:S
        d = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
        max_possible = R + 2 * (S - 1)
        conn_vals[i] = max_possible > 0 ? d / max_possible : 0.0
    end
    for alpha in 1:R
        d = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
        max_possible = S + 2 * (R - 1)
        conn_vals[S+alpha] = max_possible > 0 ? d / max_possible : 0.0
    end
    avg_conn = mean(conn_vals)

    n = S + R
    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
    end
    avg_degree = mean(degree)
    return (prop_omn = prop_omn, total_species = total_species, avg_conn = avg_conn, avg_degree = avg_degree)
end

# =============================================================================
# 2. Stable Configuration Search Functions
# =============================================================================
"""
    find_stable_configurations_real(cell; mu_range, eps_range, m_alpha_range, time_end)

Searches over μ, ε, and mₐ for a given real community cell (using analytical_equilibrium).
Returns all stable configurations (all eigenvalues of the Jacobian < 0).
"""
function find_stable_configurations_real(cell::Int; 
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:0.2, time_end=500.0, plot=false)
    
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = analytical_equilibrium(cell, mu_val, eps_val, m_alpha_val;
                    delta_nu = 0.05, d_alpha = 1.0, d_i = 1.0,
                    include_predators = true, include_omnivores = true,
                    sp_removed_name = nothing, artificial_pi = true, pi_size = 10.0,
                    H_init = nothing, P_init = nothing,
                    nu_omni_proportion = 1.0, nu_b_proportion = 1.0, r_omni_proportion = 1.0,
                    callbacks = false, plot = plot)
                if comm === nothing continue end
                u0 = vcat(comm.equilibrium.H_star, comm.equilibrium.P_star)
                params = (comm.parameters.S, comm.parameters.R, comm.parameters.K_i, comm.parameters.r_i,
                          comm.parameters.mu, comm.parameters.nu, comm.parameters.P_matrix, comm.parameters.O_matrix,
                          comm.parameters.T_matrix, comm.parameters.epsilon, comm.parameters.m_alpha,
                          comm.parameters.K_alpha, comm.parameters.B_matrix, comm.parameters.D_matrix,
                          comm.parameters.nu_omni, comm.parameters.nu_b)
                sol, J = simulate_community(u0, params, time_end)
                if all(real.(eigen(J).values) .< 0)
                    println("Cell $cell: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val)")
                    push!(stable_configs, (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val))
                end
            end
        end
    end
    if isempty(stable_configs)
        error("No stable configuration found for cell $cell")
    else
        return stable_configs
    end
end

"""
    find_stable_configurations_abstract(S, O, R; conn, mu_range, eps_range, m_alpha_range, time_end)

Searches over μ, ε, and mₐ for an abstract community defined by S, O, and R (with target connectance conn).
Returns all stable configurations.
"""
function find_stable_configurations_abstract(S::Int, O::Int, R::Int; conn=0.2,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = a_parametrise_the_community(S, O, R;
                    conn=conn, mu=mu_val, epsilon_val=eps_val, mean_m_alpha=m_alpha_val,
                    cell_abundance_h = [ i <= S ? rand()*10.0 : rand()*5.0 for i in 1:(S+O) ],
                    cell_abundance_p = [ rand() for i in 1:R ])
                u0 = vcat(comm.H_eq, comm.P_eq)
                params = (comm.S, comm.R, comm.K_i, comm.r_i, mu_val, comm.nu,
                          comm.P_matrix, comm.O_matrix, comm.T_matrix, comm.epsilon,
                          comm.m_alpha, comm.K_alpha, comm.B_matrix, comm.D_matrix,
                          comm.nu*1.0, comm.nu*1.0)
                sol, J = simulate_community(u0, params, time_end)
                if all(real.(eigen(J).values) .< 0)
                    println("Abstract community: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val), conn=$(conn)")
                    push!(stable_configs, (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val))
                end
            end
        end
    end
    if isempty(stable_configs)
        error("No stable configuration found for abstract community with S=$S, O=$O, R=$R, conn=$conn")
    else
        return stable_configs
    end
end

# =============================================================================
# 3. Metrics Computation (shared)
# =============================================================================
function compute_resilience(J)
    return maximum(real.(eigen(J).values))
end

function compute_reactivity(J)
    S_sym = (J + transpose(J)) / 2
    return maximum(real.(eigen(S_sym).values))
end

function compute_persistence(sol; threshold=1e-6)
    final_state = sol.u[end]
    n = length(final_state)
    survivors = count(x -> x > threshold, final_state)
    return survivors / n
end

# =============================================================================
# 4. Experiment Pipelines
# =============================================================================
"""
    run_experiment_real(; cell_range, ...)

For each real community (cell) in cell_range, uses find_stable_configurations_real to obtain all stable configurations,
computes stability metrics and structural metrics, and stores one row per configuration in a DataFrame.
"""
function run_experiment_real(; cell_range=1:10, time_end=500.0,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, plot=false)
    
    results = DataFrame(community_type=String[], cell=Int[],
                        S=Int[], O=Int[], R=Int[], conn=Float64[],
                        mu=Float64[], eps=Float64[], m_alpha=Float64[], nu=Float64[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[],
                        prop_omn=Float64[], total_species=Int[], avg_conn=Float64[], avg_degree=Float64[])
    lock_obj = ReentrantLock()
    
    Threads.@threads for cell in cell_range
        try
            stable_configs = find_stable_configurations_real(cell;
                                    mu_range=mu_range, eps_range=eps_range, 
                                    m_alpha_range=m_alpha_range, time_end=time_end, plot=plot)
            for (config_idx, stable) in enumerate(stable_configs)
                sol, J = stable.sol, stable.J
                resil = compute_resilience(J)
                react = compute_reactivity(J)
                persist = compute_persistence(sol)
                struct_metrics = compute_structural_metrics_real(stable.comm)
                # Extract S and O from the herbivore list:
                total_herb = length(stable.comm.equilibrium.herbivore_list)
                num_omnivores = count(sp -> sp in omnivore_names, stable.comm.equilibrium.herbivore_list)
                S_real = total_herb - num_omnivores
                O_real = num_omnivores
                R_real = stable.comm.parameters.R
                conn_val = struct_metrics.avg_conn
                lock(lock_obj) do
                    push!(results, (
                        "Real", cell, S_real, O_real, R_real, conn_val,
                        stable.mu, stable.eps, stable.m_alpha, stable.comm.parameters.nu, resil, react, persist,
                        struct_metrics.prop_omn, struct_metrics.total_species, struct_metrics.avg_conn, struct_metrics.avg_degree
                    ))
                end
            end
        catch e
            println("Skipping cell $cell: ", e)
        end
    end
    return results
end

"""
    run_experiment_abstract(; S_range, O_range, R_range, conn_range, ...)

For each abstract community configuration (given by S, O, R, and conn), uses find_stable_configurations_abstract
to obtain all stable configurations, computes metrics, and stores one row per configuration in a DataFrame.
"""
function run_experiment_abstract(; S_range=10:10:30, O_range=0:10:30, R_range=5:5:15, conn_range=0.1:0.1:1.0,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    results = DataFrame(community_type=String[], cell = Int[], S=Int[], O=Int[], R=Int[], conn=Float64[],
                        mu=Float64[], eps=Float64[], m_alpha=Float64[], nu =Float64[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[],
                        prop_omn=Float64[], total_species=Int[], avg_conn=Float64[], avg_degree=Float64[])
    for S in S_range
        for O in O_range
            for R in R_range
                for conn in conn_range
                    try
                        stable_configs = find_stable_configurations_abstract(S, O, R; conn=conn, 
                                                mu_range=mu_range, eps_range=eps_range, m_alpha_range=m_alpha_range,
                                                time_end=time_end)
                        for s in stable_configs
                            sol, J = s.sol, s.J
                            resil = compute_resilience(J)
                            react = compute_reactivity(J)
                            persist = compute_persistence(sol)
                            struct_metrics = compute_structural_metrics_abstract(s.comm)
                            push!(results, (
                                "Abstract", 1, S, O, R, conn, s.mu, s.eps, s.m_alpha, s.comm.nu,
                                resil, react, persist,
                                struct_metrics.prop_omn, struct_metrics.total_species, struct_metrics.avg_conn, struct_metrics.avg_degree
                            ))
                        end
                    catch e
                        println("Skipping configuration S=$S, O=$O, R=$R, conn=$conn: ", e)
                    end
                end
            end
        end
    end
    return results
end

# -----------------------------------------------------------------------------
# 5. Combined Experiment Pipeline
# -----------------------------------------------------------------------------
"""
    run_experiment_comparison(; cell_range, S_range, O_range, R_range, conn_range, ...)

Runs both the real-community experiment (over a range of cells) and the abstract-community experiment
(over a grid of S, O, R, and conn). Returns a combined DataFrame with a "community_type" column.
"""
function run_experiment_comparison(; cell_range=1:10,
    S_range=10:10:30, O_range=0:10:30, R_range=5:5:15, conn_range=0.08:0.02:0.18,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    real_results = run_experiment_real(; cell_range=cell_range, time_end=time_end,
                                        mu_range=mu_range, eps_range=eps_range, m_alpha_range=m_alpha_range)
    abstract_results = run_experiment_abstract(; S_range=S_range, O_range=O_range, R_range=R_range,
                                                conn_range=conn_range, mu_range=mu_range, eps_range=eps_range,
                                                m_alpha_range=m_alpha_range, time_end=time_end)
    combined = vcat(real_results, abstract_results)
    return combined
end

# -----------------------------------------------------------------------------
# 6. Plotting the Combined Results
# -----------------------------------------------------------------------------
function plot_combined_results(
    results::DataFrame;
    variable_to_color_by::Symbol = :mu,
    log_resilience::Bool = false,
    log_reactivity::Bool = false
)
    cmap = cgrad(:viridis)  # Use Viridis colormap

    # Compute overall color range based on the chosen variable.
    cmin = minimum(results[!, variable_to_color_by])
    cmax = maximum(results[!, variable_to_color_by])
    
    # Separate indices for Real and Abstract communities.
    real_idx = findall(results.community_type .== "Real")
    abs_idx  = findall(results.community_type .== "Abstract")
    
    fig = Figure(resolution=(1000,600))
    
    if log_resilience
        real_res = log.(results.resilience[real_idx].*-1)
        abs_res = log.(results.resilience[abs_idx].*-1)
    else
        real_res = results.resilience[real_idx]
        abs_res = results.resilience[abs_idx]
    end
    
   if log_reactivity
        real_react = results.reactivity[real_idx]
        abs_react = results.reactivity[abs_idx]

        for (i, val) in enumerate(real_react)
            if val > 0.0
                real_react[i] = log(val)
            elseif val < 0.0
                real_react[i] = -log(-val)
            else
                real_react[i] = 0.0  # handle log(0) case if needed
            end
        end
        for (i, val) in enumerate(abs_react)
            if val > 0.0
                abs_react[i] = log(val)
            elseif val < 0.0
                abs_react[i] = -log(-val)
            else
                abs_react[i] = 0.0  # handle log(0) case if needed
            end
        end
    else
        real_react = results.reactivity[real_idx]
        abs_react = results.reactivity[abs_idx]
    end

    # Plot 1: Resilience vs. Average Connectance
    ax1 = Axis(fig[1,1], xlabel="Average Connectance", ylabel=log_resilience ? "Log 1/Resilience" : "Resilience", title="Resilience vs. Connectance")
    scatter!(ax1, results.avg_conn[real_idx], real_res,
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax1, results.avg_conn[abs_idx], abs_res,
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    
    # Plot 2: Reactivity vs. Proportion Omnivores
    ax2 = Axis(fig[1,2], xlabel="Proportion Omnivores", ylabel=log_reactivity ? "Log Reactivity" : "Reactivity", title="Reactivity vs. Proportion Omnivores")
    scatter!(ax2, results.prop_omn[real_idx], real_react,
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax2, results.prop_omn[abs_idx], abs_react,
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    
    # Plot 3: Persistence vs. Total Species
    ax3 = Axis(fig[2,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Total Species")
    scatter!(ax3, results.total_species[real_idx], results.persistence[real_idx],
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax3, results.total_species[abs_idx], results.persistence[abs_idx],
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    
    # Plot 4: Resilience vs. Average Degree
    ax4 = Axis(fig[2,2], xlabel="Average Degree", ylabel=log_resilience ? "Log 1/Resilience" : "Resilience", title="Resilience vs. Average Degree")
    scatter!(ax4, results.avg_degree[real_idx], real_res,
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax4, results.avg_degree[abs_idx], abs_res,
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    
    # Plot 5: Reactivity vs. Average Connectance
    ax5 = Axis(fig[3,1], xlabel="Average Connectance", ylabel=log_reactivity ? "Log Reactivity" : "Reactivity", title="Reactivity vs. Connectance")
    scatter!(ax5, results.avg_conn[real_idx], real_react,
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax5, results.avg_conn[abs_idx], abs_react,
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))

    # Plot 6: Resilience vs. Total Species
    ax6 = Axis(fig[3,2], xlabel="Total Species", ylabel=log_resilience ? "Log 1/Resilience" : "Resilience", title="Resilience vs. Total Species")
    scatter!(ax6, results.total_species[real_idx], real_res,
             markersize=12, marker=:circle, color=results[real_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    scatter!(ax6, results.total_species[abs_idx], abs_res,
             markersize=8, marker=:utriangle, color=results[abs_idx, variable_to_color_by],
             colormap=cmap, colorrange=(cmin, cmax))
    
    # Plot 6: Colorbar
    Colorbar(fig[3,3], limits=(cmin, cmax), colormap=cmap)
    
    display(fig)
end


# -----------------------------------------------------------------------------
# 7. Running the Comparison Pipeline
# -----------------------------------------------------------------------------
combined_results_with_nu_good = run_experiment_comparison(;
    cell_range=1:10,
    S_range=20:10:50, O_range=0:5:10, R_range=5:5:10, conn_range=0.0:0.05:0.2,
    mu_range=0.0:0.25:1.0, eps_range=0.0:0.5:1.0, m_alpha_range=0.05:0.1:0.25,
    time_end=500.0
)
combined_results = CSV.File("combined_results.csv") |> DataFrame

# println(combined_results[1:10, :])
plot_combined_results(
    combined_results_with_nu_good;
    variable_to_color_by = :eps,
    log_resilience = false,
    log_reactivity = false
)

using GLM

# Fit a linear model where resilience is explained by all other variables
lm_model = lm(@formula(resilience ~ (prop_omn + total_species + avg_conn + avg_degree + reactivity + mu + eps + m_alpha)^2), combined_results)
# Display the summary of the linear model
display(coeftable(lm_model))

# Save the combined results to a CSV file
CSV.write("combined_results.csv", combined_results)
