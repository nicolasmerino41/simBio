function omnivore_run(
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
    nu_omni_proportion = 1.0
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
        d_i = d_i
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
    # println("Herbivore list: ", herbivore_list)

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
    params = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_val*nu_omni_proportion)

    # Define and solve the ODE using omnivore_dynamics!
    prob = ODEProblem(omnivore_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
    end

    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: solution did not finish properly"
            return nothing
        end
    end

    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (600, 500))
        ax = Axis(
            fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for cell $cell",
            yscale = log ? log10 : identity
        )
        times_combined = sol.t
        # Plot herbivore/omnivore dynamics: species that are in herbivore_names are blue; those in omnivore_names are green.
        for i in 1:S
            if herbivore_list[i] in omnivore_names
                lines!(ax, times_combined, sol[i, :], label = "O$(i)", color = :green)
                # println("O$(i)")
            else
                lines!(ax, times_combined, sol[i, :], label = "H$(i)", color = :blue)
                # println("H$(i)")
            end
        end
        # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
        for alpha in 1:R
            lines!(ax, times_combined, sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
        end
        display(fig)
    end
    println("S: $S")
    # Post-processing: extract final densities.
    H_end = sol[1:S, end]
    P_end = sol[S+1:S+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    survival_rate = total_surv / total_species
    herbivore_survival_rate = survived_herb / S
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    println("Estimated nu = $nu_val")
    
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
        H_eq                    = [H_eq],
        P_eq                    = [P_eq],
        H_star                  = [H_star],
        P_star                  = [P_star],
        H_vector                = [H_end],
        P_vector                = [P_end],
        survived_herbivores     = survived_herb,
        survived_predators      = survived_pred,
        H_biomass               = H_biomass,
        P_biomass               = P_biomass,
        biomass_at_the_end      = biomass_at_the_end,
        herb_pred_ratio         = ratio,
        total_survivors         = total_surv,
        total_species           = total_species,
        herbivore_survival_rate = herbivore_survival_rate,
        predator_survival_rate  = predator_survival_rate,
        mu                      = mu_val,
        nu                      = nu_val,
        epsilon_val             = eps_val,
        mean_m_alpha            = mean_m_alpha
    )
    
    append!(AAAA, single_run_results)
    
    if do_you_want_params && do_you_want_sol
        return AAAA, params, sol
    elseif do_you_want_params || do_you_want_sol
        return AAAA, do_you_want_params ? params : sol
    else
        return AAAA
    end
end

# cb_no_trigger, cb_trigger = build_callbacks(33, 12, EXTINCTION_THRESHOLD, T_ext, 1)
@time A_run = omnivore_run(
    1, # cell
    0.5, 1.0, 0.01; # mu, epsilon, m_alpha
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 1000.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    include_omnivores = false,
    plot = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 10.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = true,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 0,
    include_n_omnivores = 0,
    nu_omni_proportion = 1.0
)

println("Survival rate: ", A_run.survival_rate[1])

println(sol[:, end])
println(sol.u[:, end])
