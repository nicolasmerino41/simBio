function total_omnivore_run(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 500.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    include_omnivores = true,
    plot = false,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 0,
    include_n_omnivores = 0,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 0.01,
    force_nu_to = nothing,
    use_cb = false,
    removal_fraction = 1.0
)
    AAAA = DataFrame()
    local_i, local_j = idx[cell][1], idx[cell][2]
    # Extract species present in the cell.
    all_sp = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Start with the full species list.
    sp_nm = copy(all_sp)
    
    # Remove predators if not included.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if !(sp in carnivore_names)]
    end

    # Process omnivores:
    if include_omnivores
        # If include_omnivores is true, then all omnivores remain in sp_nm.
        if include_n_omnivores > 0
            error("Cannot specify include_n_omnivores > 0 when include_omnivores is true")
        end
    else
        # Exclude all omnivores.
        sp_nm = [sp for sp in sp_nm if !(sp in omnivore_names)]
        # But if a positive number of omnivores is requested, add them from the cell.
        if include_n_omnivores > 0
            omnivores_in_cell = [sp for sp in all_sp if sp in omnivore_names]
            if length(omnivores_in_cell) < include_n_omnivores
                error("Not enough omnivores in the cell to include the requested number: requested $(include_n_omnivores), found $(length(omnivores_in_cell))")
            end
            sp_nm = union(sp_nm, omnivores_in_cell[1:include_n_omnivores])
        end
    end

    # Also remove species if sp_removed_name is given.
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Get community parameters using the new setup function.
    params_setup = o_attempt_setup_community(
        local_i, local_j,
        mu_val,
        eps_val,
        mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure the returned parameters.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq        # Baseline herbivore (and omnivore) abundances
    P_eq       = params_setup.P_eq        # Baseline predator abundances
    r_i        = params_setup.r_i          # Herbivore intrinsic growth rates
    K_i        = params_setup.K_i          # Herbivore carrying capacities
    nu_val     = params_setup.nu           # Effective predation rate
    P_matrix   = params_setup.P_matrix     # S x R incidence matrix (carnivores preying on herbivores)
    O_matrix   = params_setup.O_matrix     # S x S omnivory interaction matrix (consumer to prey)
    T_matrix   = params_setup.T_matrix     # S x S herbivore-herbivore consumption matrix (prey to consumer)
    epsilon    = params_setup.epsilon      # Predator conversion efficiencies
    m_alpha    = params_setup.m_alpha      # Predator mortality rates
    K_alpha    = params_setup.K_alpha      # Predator carrying capacities
    B_matrix   = params_setup.B_matrix     # R x R beneficial predator-predator interaction matrix
    D_matrix   = params_setup.D_matrix     # R x R detrimental predator-predator interaction matrix
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    println("Herbivore list: ", length(herbivore_list))

    # Initial conditions: use provided or baseline abundances.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        if iszero(extinguishing_species)
            u0 = abs.(vcat(H_init, P_init) .+ randn(S+R) .* initial_deviation)
        else
            if extinguishing_species > S+R
                @error "Error: extinguishing_species > S+R"                
            end
            u0 = vcat(H_init, P_init) .* [i == extinguishing_species ? 0.0 : 1.0 for i in 1:(S+R)]
        end
    else
        u0 = H_init
    end
    
    # Build the parameter tuple for the ODE.
    # Order: (S, R, K_i, r_i, mu, nu, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix)
    if !isnothing(force_nu_to)
        nu_val = force_nu_to
    end
    params = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_val*nu_omni_proportion, nu_val*nu_b_proportion)

    # Define and solve the ODE using omnivore_dynamics!
    prob = ODEProblem(omnivore_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if use_cb
            solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
        else
            solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    end

    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: solution did not finish properly"
            return nothing
        end
    end

    # --- PERTURBATION HALFWAY ---
    fig = Figure(; size = (900, 700))
    ax = Axis(
        fig[1, 1],
        xlabel = "Time",
        ylabel = "Biomass",
        title = "Dynamics for cell $cell",
        yscale = log ? log10 : identity
    )
    for i in 1:S+R
        new_u0 = deepcopy(u0)
        new_u0[i] -= new_u0[i] * removal_fraction
        println("new_u0: ", new_u0)
        prob = ODEProblem(omnivore_dynamics!, new_u0, (time_end, time_end+time_end), params)
        if use_cb
            new_sol = solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
            else
            new_sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
        times_combined = vcat(sol.t, new_sol.t)
        if i <= S
            lines!(ax, times_combined, sum(eachrow(hcat(sol[:, :], new_sol[:, :]))), color = :blue)
            # println("For species $i, the largest difference is ", maximum(abs.(sum(sol[:, end]) .- sum(eachrow(new_sol[:, :])))))
            println("For species $i the biomass before perturbation is ", sum(sol[:, end]), " and the new is ", sum(new_sol[:, end]))
        else
            lines!(ax, times_combined, sum(eachrow(hcat(sol[:, :], new_sol[:, :]))), color = :red)
            # println("For species $i, the largest difference is ", maximum(abs.(sum(sol[:, end]) .- sum(eachrow(new_sol[:, :])))))
            println("For species $i the biomass before perturbation is ", sum(sol[:, end]), " and the new is ", sum(new_sol[:, end]))
        end
    end

    display(fig)
    return params
end

@time A_run = total_omnivore_run(
    1, # cell
    0.5, 0.29, 0.1; # mu, epsilon, m_alpha
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 1000.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    include_omnivores = true,
    plot = true,
    sp_removed_name = nothing,
    artificial_pi = true, pi_size = 10.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = true,
    log = false,
    initial_deviation = 0.0, # leave at 0.0 for no effect
    extinguishing_species = 0, # leave at 0 for no effect
    include_n_omnivores = 0, # leave at 0 for no effect
    nu_omni_proportion = 1.0, # leave at 1.0 for no effect
    nu_b_proportion = 1.0, # leave at 1.0 for no effect
    r_omni_proportion = 1.0, # leave at 1.0 for no effect
    force_nu_to = nothing, # leave at nothing for no effect
    use_cb = false,
    removal_fraction = 0.5
)