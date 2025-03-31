##### PLOTING EIGENVALUES FOR MULTIPLE CELLS #####
# --- Merged Function ---
function plot_multiple_cells_eigenvalues(num_cells;
    callbacks=true,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.01:1.0, m_alpha_range=0.0:0.1:1.0,
    cut_value = nothing
)
    
    re_vals_all = Float64[]  # initialize as a concrete vector
    im_vals_all = Float64[]
    largest_re_all = DataFrame(
        Cell = Int[],
        Largest_re = Float64[]
        )

    lock_obj = ReentrantLock()  # create a lock object

    Threads.@threads for cell in 1:num_cells
        println("Processing cell $cell")
        stable_found = false
        stable_result = nothing

        # Search for a stable equilibrium configuration over the parameter ranges.
        for eps_val in eps_range
            for mu_val in mu_range
                for m_alpha in m_alpha_range
                    try
                        result = analytical_equilibrium(cell, mu_val, eps_val, m_alpha;
                        delta_nu = 0.05,
                        d_alpha = 1.0, d_i = 1.0,
                        include_predators = true,
                        include_omnivores = true,
                        sp_removed_name = nothing,
                        artificial_pi = true, pi_size = 10.0,
                        H_init = nothing, P_init = nothing,
                        nu_omni_proportion = 1.0, nu_b_proportion = 1.0,
                        r_omni_proportion = 1.0,
                        callbacks = callbacks, plot = false)
                    
                    eigvals = eigen(result.Jacobian).values
                    if all(real.(eigvals) .< 0)
                        println("Stable equilibrium found for cell $cell with μ = $mu_val, ε = $eps_val, mₐ = $m_alpha")
                        stable_found = true
                        stable_result = result
                        break
                    end
                    catch e
                        println("Error for cell $cell with μ = $mu_val, ε = $eps_val, mₐ = $m_alpha: $e")
                        continue
                    end
                end
                if stable_found
                    break
                end
            end
            if stable_found
                break
            end
        end
        
        if !stable_found
            println("No stable equilibrium found for cell $cell")
            return nothing
        end
    
        if !stable_found
            println("No stable equilibrium found for cell $cell")
            continue  # instead of return nothing, skip this cell
        end

        # Compute eigenvalues
        ev = eigen(stable_result.Jacobian)
        vals = ev.values
        re_vals = real.(vals)
        im_vals = imag.(vals)

        # Lock and update shared arrays
        lock(lock_obj) do
            append!(re_vals_all, re_vals)
            append!(im_vals_all, im_vals)

            largest_re = maximum(re_vals)
            largest_re_all = vcat(largest_re_all, DataFrame(Cell = cell, Largest_re = largest_re))
        end
    end

    # println("re_vals_all: ", re_vals_all)
    if !isnothing(cut_value)
        filter!(x -> x >= -cut_value, re_vals_all)
    end

    # fig = Figure()
    # ax = Axis(
    #     fig[1, 1], xlabel = "Real Part", ylabel = "Imaginary Part",
    #     title = "Eigenvalue Plot by Guild"
    # )
    # MK.scatter!(ax, re_vals_all, im_vals_all[1:length(re_vals_all)], color = :blue, marker = :circle)
    # MK.hlines!(ax, [0], linestyle = :dash, color = :black)
    # MK.vlines!(ax, [0], linestyle = :dash, color = :black)
    # display(fig)
    # println("All eigenvalues: ", vcat(re_vals_all, im_vals_all))
    # return vcat(re_vals_all, im_vals_all)

    return largest_re_all

end

largest_re_all = plot_multiple_cells_eigenvalues(
    200;
    callbacks=false,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.01:1.0, m_alpha_range=0.0:0.1:1.0,
    cut_value = 5
)


#################### EIGENVALUES VS METRICS ####################
#######################################################################################
#######################################################################################
# Community-level metrics function
function compute_community_metrics(cell; 
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.01:1.0, m_alpha_range=0.0:0.1:1.0,
    callbacks=false)

    println("Processing cell $cell")
    stable_found = false
    stable_result = nothing

    # Search for a stable equilibrium configuration over the parameter ranges.
    for eps_val in eps_range
        for mu_val in mu_range
            for m_alpha in m_alpha_range
                try
                    result = analytical_equilibrium(cell, mu_val, eps_val, m_alpha;
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = false, pi_size = 10.0,
                    H_init = nothing, P_init = nothing,
                    nu_omni_proportion = 1.0, nu_b_proportion = 1.0,
                    r_omni_proportion = 1.0,
                    callbacks = callbacks, plot = false)
                
                eigvals = eigen(result.Jacobian).values
                if all(real.(eigvals) .< 0)
                    println("Stable equilibrium found for cell $cell with μ = $mu_val, ε = $eps_val, mₐ = $m_alpha")
                    stable_found = true
                    stable_result = result
                    break
                end
                catch e
                    println("Error for cell $cell with μ = $mu_val, ε = $eps_val, mₐ = $m_alpha: $e")
                    continue
                end
            end
            if stable_found
                break
            end
        end
        if stable_found
            break
        end
    end
    
    # Compute eigenvalues
   
    # Extract equilibrium state, parameters, and Jacobian.
    A_eq = stable_result.equilibrium   # Contains u0, herbivore_list, predator_list, etc.
    A_p  = stable_result.parameters    # Contains S, R, interaction matrices, etc.
    J    = stable_result.Jacobian

    # 1. Compute largest real eigenvalue (largest_re):
    eigvals = eigen(J).values
    largest_re = maximum(real.(eigvals))

    # 2. Compute total community biomass:
    total_biomass = sum(A_eq.u0)

    # 3. Compute average degree:
    S = A_p.S
    R = A_p.R
    n = S + R
    P_matrix = A_p.P_matrix
    O_matrix = A_p.O_matrix
    T_matrix = A_p.T_matrix
    B_matrix = A_p.B_matrix
    D_matrix = A_p.D_matrix

    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
    end
    avg_degree = sum(degree) / n

    # 4. Compute average connectance:
    # For herbivores:
    connectance_arr = zeros(n)
    for i in 1:S
        d = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
        max_possible = R + 2 * (S - 1)
        connectance_arr[i] = d / max_possible
    end
    # For predators:
    for alpha in 1:R
        d = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
        max_possible = S + 2 * (R - 1)
        connectance_arr[S+alpha] = d / max_possible
    end
    avg_connectance = sum(connectance_arr) / n

    return (largest_re = largest_re,
            biomass = total_biomass,
            avg_degree = avg_degree,
            avg_connectance = avg_connectance)
end
begin
    
    # Now, assume you have a list of cell indices for which you want to compute these metrics:
    cells = 1:50   # for example, cells 1 to 10

    # We'll store metrics in arrays
    largest_re_array = Float64[]
    biomass_array = Float64[]
    degree_array = Float64[]
    connectance_array = Float64[]

    for cell in cells
        # You can adjust parameter values as needed.
        metrics = compute_community_metrics(cell; callbacks=false)
        push!(largest_re_array, metrics.largest_re)
        push!(biomass_array, metrics.biomass)
        push!(degree_array, metrics.avg_degree)
        push!(connectance_array, metrics.avg_connectance)
    end

    # Now, produce three scatter plots:
    fig = Figure(resolution = (1200, 400))

    ax1 = Axis(fig[1,1], xlabel = "Largest Real Eigenvalue", ylabel = "Total Biomass",
            title = "Biomass vs. Largest Real Eigenvalue")
    scatter!(ax1, largest_re_array, biomass_array, markersize=10)

    ax2 = Axis(fig[1,2], xlabel = "Largest Real Eigenvalue", ylabel = "Avg. Degree",
            title = "Degree vs. Largest Real Eigenvalue")
    scatter!(ax2, largest_re_array, degree_array, markersize=10)

    ax3 = Axis(fig[1,3], xlabel = "Largest Real Eigenvalue", ylabel = "Avg. Connectance",
            title = "Connectance vs. Largest Real Eigenvalue")
    scatter!(ax3, largest_re_array, connectance_array, markersize=10)

    display(fig)
end