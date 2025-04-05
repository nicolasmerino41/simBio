function find_stable_configurations_real_stop(cell::Int; 
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:0.2, time_end=500.0, plot=false,
    stop_early=false)
    
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
                    # println("Cell $cell: Stable config found: μ=$(mu_val), ε=$(eps_val), mₐ=$(m_alpha_val)")
                    push!(stable_configs, (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val))
                    if stop_early
                        return stable_configs
                    end
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

function map_real_resilience(; cell_range=1:10, time_end=500.0,
mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0)

    results = DataFrame(community_type=String[], cell=Int[],
                        S=Int[], O=Int[], R=Int[], conn=Float64[],
                        mu=Float64[], eps=Float64[], m_alpha=Float64[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[],
                        prop_omn=Float64[], total_species=Int[], avg_conn=Float64[], avg_degree=Float64[])
    lock_obj = ReentrantLock()

    Threads.@threads for cell in cell_range
        try
            stable_configs = find_stable_configurations_real_stop(cell;
                                    mu_range=mu_range, eps_range=eps_range, 
                                    m_alpha_range=m_alpha_range, time_end=time_end, stop_early=true)
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
                        stable.mu, stable.eps, stable.m_alpha, resil, react, persist,
                        struct_metrics.prop_omn, struct_metrics.total_species, struct_metrics.avg_conn, struct_metrics.avg_degree
                    ))
                end
            end
        catch e
            # println("Skipping cell $cell: ", e)
        end
    end
    DA_metric = deepcopy(float.(DA_sum))
    DA_metric[DA_metric .== 0] .= NaN
    DA_metric[DA_metric .== 1] .= NaN
    for cell in results.cell
        resil = results.resilience[results.cell .== cell]
        DA_metric[idx[cell][1], idx[cell][2]] = resil[1]
    end

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "DA_metric Visualization", xlabel = "X-axis", ylabel = "Y-axis")
    heatmap!(ax, DA_metric, colormap = :viridis)
    ax.yreversed[] = true
    display(fig)

end

map_real_resilience(
    ; cell_range=1:5950, time_end=500.0,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:1.0
)
