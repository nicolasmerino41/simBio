function bipartite_run(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 500.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    plot = false,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 0
)
    AAAA = DataFrame()
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species: if predators are not included or if a species is to be removed.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Get community parameters using the new setup function.
    params_setup = b_attempt_setup_community(
        local_i, local_j,
        mu_val,
        # nu = nu_val,
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
    H_eq       = params_setup.H_eq        # baseline herbivore abundances
    P_eq       = params_setup.P_eq        # baseline predator abundances
    r_i        = params_setup.r_i          # herbivore intrinsic growth rates
    K_i        = params_setup.K_i          # herbivore carrying capacities
    # mu here is the input parameters (should equal mu_val and nu_val)
    nu_val     = params_setup.nu
    P_matrix   = params_setup.P_matrix     # S x R predation incidence matrix
    epsilon    = params_setup.epsilon      # predator conversion efficiencies
    m_alpha    = params_setup.m_alpha      # predator mortality
    K_alpha    = params_setup.K_alpha      # predator carrying capacities
    # herbivore_list = params_setup.herbivore_list
    # predator_list  = params_setup.predator_list
    # species_names = params_setup.species_names,
    H_star = params_setup.H_star
    P_star = params_setup.P_star
    # println("H_star: $H_star")
    # println("P_star: $P_star")

    # Initial conditions: if H_init is not provided, use the baseline abundances.
    if isnothing(H_init)
        H_init = H_star
    end
    # Set up predator initial conditions if needed.
    if R > 0
        if isnothing(P_init)
            # Default: use a fraction of the first herbivore's abundance.
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
        # println("typeof(H_init): ", typeof(H_init))
        # println("typeof(P_init): ", typeof(P_init))
    else
        u0 = H_init
    end
    # println("typeof(u0): ", typeof(u0))
    # Build the parameter tuple for the ODE.
    # Order: (S, R, K_i, r_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha)
    params = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, epsilon, m_alpha, K_alpha)
    # println("K_i: $K_i")
    # println("K_alpha: $K_alpha")
    # println("r_i: $r_i")
    # println("m_alpha: $m_alpha")
    # Define and solve the ODE.
    prob = ODEProblem(bipartite_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
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
        # Plot herbivore dynamics (indices 1:S) in blue.
        for i in 1:S
            lines!(ax, times_combined, sol[i, :], label = "H$(i)", color = :blue)
        end
        # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
        for alpha in 1:R
            lines!(ax, times_combined, sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
        end
        display(fig)
    end
    
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
    # giHi = sum(r_i .* H_end)
    herbivore_survival_rate = survived_herb / S
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    # println("H_end: $H_end")
    # println("P_end: $P_end")
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

# cell = 1
# cb_no_trigger, cb_trigger = build_callbacks(33, 12, EXTINCTION_THRESHOLD, T_ext, 1)
@time A_run = bipartite_run(
    1, # cell
    1.0, 1.0, 0.01; # mu, epsilon, m_alpha
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 1000.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    plot = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 10.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = true,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 17 
)

println("Survival rate: ", A_run.survival_rate[1])

println(sol[:, end])
println(sol.u[:, end])