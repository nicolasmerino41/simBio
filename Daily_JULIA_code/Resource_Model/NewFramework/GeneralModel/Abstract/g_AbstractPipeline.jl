# =============================================================================
# 1. ODE Dynamics and Simulation (using your omnivore_dynamics!)
# =============================================================================
function simulate_community(u0, params, time_end=500.0; initial_deviation=0.0)
    u0 = u0 .+ randn(length(u0)) .* initial_deviation
    prob = ODEProblem(general_dynamics!, u0, (0.0, time_end), params)
    sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    function f_wrapper(u)
        du = similar(u)
        general_dynamics!(du, u, params, 0.0)
        return du
    end
    J = ForwardDiff.jacobian(f_wrapper, sol.u[end])
    return sol, J
end

# =============================================================================
# 2. Metrics Computation
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
# 3. Structural Metrics Computation for Abstract Communities
# =============================================================================
function compute_structural_metrics_abstract(comm)
    # In our abstract community, the herbivore compartment is defined by comm.herbivore_list
    # and the predator compartment by comm.predator_list.
    all_herb = comm.herbivore_list
    # Identify omnivores (assumed to have "Omnivore" in their name)
    actual_O = [sp for sp in all_herb if occursin("Omnivore", sp)]
    prop_omn = isempty(all_herb) ? 0.0 : length(actual_O) / length(all_herb)
    total_species = length(all_herb) + length(comm.predator_list)
    
    S = comm.S    # Number of herbivore compartment species (S_total)
    R = comm.R
    A = comm.A_matrix    # Merged interaction matrix (size (S+R)×(S+R))
    
    # Compute connectance for each species.
    conn_vals = zeros(S + R)
    # For herbivores (rows 1:S), define maximum possible interactions as:
    # max_possible = R + 2*(S-1)
    for i in 1:S
        d = sum(abs.(A[i, :]))
        max_possible = R + 2*(S - 1)
        conn_vals[i] = max_possible > 0 ? d / max_possible : 0.0
    end
    # For predators (rows S+1:S+R), maximum possible is:
    # max_possible = S + 2*(R-1)
    for alpha in 1:R
        idx = S + alpha
        d = sum(abs.(A[idx, :]))
        max_possible = S + 2*(R - 1)
        conn_vals[idx] = max_possible > 0 ? d / max_possible : 0.0
    end
    avg_conn = mean(conn_vals)
    
    # Compute the degree for each species as the sum of the absolute interaction strengths.
    degree = [sum(abs.(A[i, :])) for i in 1:(S + R)]
    avg_degree = mean(degree)
    
    return (prop_omn = prop_omn, total_species = total_species, avg_conn = avg_conn, avg_degree = avg_degree)
end

# =============================================================================
# 4. Stable Configuration Search for Abstract Communities
# =============================================================================
"""
    find_stable_configurations_abstract(S, O, R; conn, mu_range, eps_range, m_alpha_range, time_end)

Searches over μ, ε, and mₐ for an abstract community defined by S, O, and R (with target connectance conn).
Returns all configurations yielding a locally stable equilibrium.
"""
function find_stable_configurations_abstract(S::Int, O::Int, R::Int; conn=0.2,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = g_abstract_parametrise_the_community(S, O, R;
                    conn=conn, mu=mu_val, epsilon_val=eps_val, mean_m_alpha=m_alpha_val,
                    r_omni_proportion = 1.0,
                    cell_abundance_h = [i <= S ? 10.0 : 5.0 for i in 1:(S+O)],
                    cell_abundance_p = [1.0 for i in 1:R])
                u0 = vcat(comm.H_eq, comm.P_eq)
                params = (comm.S, comm.R, comm.K_i, comm.r_i, mu_val, comm.nu,
                          comm.A_matrix, comm.C_matrix, comm.epsilon,
                          comm.m_alpha, comm.K_alpha)
                sol, J = simulate_community(u0, params, time_end)
                if all(real.(eigen(J).values) .< 0)
                    # println("Abstract: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val), conn=$(conn)")
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
# 5. Experiment Pipeline for Abstract Communities
# =============================================================================
"""
    run_experiment_abstract(; S_range, O_range, R_range, conn_range, mu_range, eps_range, m_alpha_range, time_end)

Loops over combinations of S, O, R, and conn for abstract communities.
For each configuration, obtains all stable configurations, computes metrics (resilience, reactivity, persistence),
and stores one row per configuration in a DataFrame.
"""
function run_experiment_abstract(; S_range=10:10:30, O_range=0:10:30, R_range=5:5:15, conn_range=0.1:0.1:1.0,
    mu_range=0.1:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    results = DataFrame(community_type=String[], S=Int[], O=Int[], R=Int[], conn=Float64[],
                        mu=Float64[], eps=Float64[], m_alpha=Float64[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[],
                        prop_omn=Float64[], total_species=Int[], avg_conn=Float64[], avg_degree=Float64[])
    
    
    lock_obj = ReentrantLock()

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
                            lock(lock_obj) do
                                push!(results, (
                                    "Abstract", S, O, R, conn, s.mu, s.eps, s.m_alpha,
                                    resil, react, persist,
                                    struct_metrics.prop_omn, struct_metrics.total_species, struct_metrics.avg_conn, struct_metrics.avg_degree
                                ))
                            end
                        end
                    catch e
                        # println("Skipping configuration S=$S, O=$O, R=$R, conn=$conn: ", e)
                    end
                end
            end
        end
    end
    return results
end

# =============================================================================
# 6. Plotting the Abstract Results
# =============================================================================
"""
    plot_experiment_results_abstract(results; variable_to_color_by)

Plots abstract community results using:
- Triangles for markers,
- Colors determined by the chosen variable (default is :mu) using the Viridis colormap.
Subplots include comparisons of:
  - Resilience vs. Average Connectance,
  - Reactivity vs. Proportion Omnivores,
  - Persistence vs. Total Species,
  - Resilience vs. Average Degree,
  - Reactivity vs. Average Connectance.
"""
function plot_experiment_results_abstract(results::DataFrame; variable_to_color_by::Symbol = :mu, log_resilience=false)
    cmap = cgrad(:viridis)
    # Compute overall color range from the chosen variable.
    cmin = minimum(results[!, variable_to_color_by])
    cmax = maximum(results[!, variable_to_color_by])
    
    # Since all results here are abstract, we use triangles.
    fig = Figure(resolution=(1000,600))

    if log_resilience
        res = log.(results.resilience.*-1)
    else
        res = results.resilience
    end
    
    ax1 = Axis(fig[1,1], xlabel="Average Connectance", ylabel=log_resilience ? "Log 1/Resilience" : "Resilience", title="Resilience vs. Connectance")
    scatter!(ax1, results.avg_conn, res, markersize=10, marker=:utriangle,
             color=results[!, variable_to_color_by], colormap=cmap, colorrange=(cmin, cmax))
    
    ax2 = Axis(fig[1,2], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
    scatter!(ax2, results.prop_omn, results.reactivity, markersize=10, marker=:utriangle,
             color=results[!, variable_to_color_by], colormap=cmap, colorrange=(cmin, cmax))
    
    ax3 = Axis(fig[2,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Total Species")
    scatter!(ax3, results.total_species, results.persistence, markersize=10, marker=:utriangle,
             color=results[!, variable_to_color_by], colormap=cmap, colorrange=(cmin, cmax))
    
    ax4 = Axis(fig[2,2], xlabel="Average Degree", ylabel=log_resilience ? "Log 1/Resilience" : "Resilience", title="Resilience vs. Average Degree")
    scatter!(ax4, results.avg_degree, res, markersize=10, marker=:utriangle,
             color=results[!, variable_to_color_by], colormap=cmap, colorrange=(cmin, cmax))
    
    ax5 = Axis(fig[3,1], xlabel="Average Connectance", ylabel="Reactivity", title="Reactivity vs. Average Connectance")
    scatter!(ax5, results.avg_conn, results.reactivity, markersize=10, marker=:utriangle,
             color=results[!, variable_to_color_by], colormap=cmap, colorrange=(cmin, cmax))
    
    # Add a colorbar.
    # Colorbar(fig[3,2], limits=(cmin, cmax), colormap=cmap)
    display(fig)
end

# =============================================================================
# 7. Running the Abstract-Community Pipeline
# =============================================================================
abstract_results = run_experiment_abstract(; 
    S_range=35, O_range=5, R_range=8, conn_range=0.1,
    mu_range=0.5, eps_range=0.29, m_alpha_range=0.1, time_end=500.0)
println(abstract_results)
plot_experiment_results_abstract(
    abstract_results;
    variable_to_color_by=:mu,
    log_resilience=false
)

#####################################################################################
#####################################################################################
##################### THIS WILL NEED TO BE REMOVED AFTER EXPLORATION ################
#####################################################################################
#####################################################################################
# function find_stable_configurations_abstract(S::Int, O::Int, R::Int; conn=0.2,
#     mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0, time_end=500.0)
    
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = g_abstract_parametrise_the_community(S, O, R;
                    conn=conn, mu=mu_val, epsilon_val=eps_val, mean_m_alpha=m_alpha_val,
                    r_omni_proportion = 1.0,
                    cell_abundance_h = [i <= S ? 10.0 : 5.0 for i in 1:(S+O)],
                    cell_abundance_p = [1.0 for i in 1:R])
                u0 = vcat(comm.H_eq, comm.P_eq)
                params = (comm.S, comm.R, comm.K_i, comm.r_i, mu_val, comm.nu,
                          comm.A_matrix, comm.C_matrix, comm.epsilon,
                          comm.m_alpha, comm.K_alpha)
                sol, J = simulate_community(u0, params, time_end)
                if all(real.(eigen(J).values) .< 0)
                    println("Abstract: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val), conn=$(conn)")
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
# end

find_stable_configurations_abstract(35, 5, 8; conn=0.1,
    mu_range=0.2, eps_range=0.2, m_alpha_range=0.15, time_end=500.0
)