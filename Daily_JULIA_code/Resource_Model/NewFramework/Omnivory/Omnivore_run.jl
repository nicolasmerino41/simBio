################################################################################
# o_omnivore_run (modified run function for the new framework)
################################################################################
function omnivore_run(
    cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
    time_end=500.0, 
    do_you_want_params=false,
    do_you_want_sol = false,
    include_predators=false,
    plot=false, sp_removed_name=nothing,
    artificial_pi=false,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    o_pred_avg::Float64 = 0.0,
    log = false
)
    AAAA = DataFrame()
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    
    if !include_predators
        local_R = 0
        filter!(name -> !(name in predator_names), sp_nm)
    end
    if !isnothing(sp_removed_name)
        filter!(name -> (name != sp_removed_name), sp_nm)
    end

    localNPP = 0.0   # not used

    # Get the community parameters from the new setup function.
    results = o_new_attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
        species_names = sp_nm,
        o_pred_avg = o_pred_avg,
        artificial_pi = artificial_pi
    )
    if isnothing(results)
        @error "Error: results === nothing"
        return nothing
    end
    
    # Destructure the returned tuple.
    # The new community setup function returns:
    # (S, R, H_i0_eff, m_i, g_i, O_loss, O_gain, A_star, a_matrix, A, epsilon, m_alpha,
    #  herbivore_list, predator_list, species_names)
    S, R, H_i0_eff, m_i, g_i, O_loss, O_gain, A_star, a_matrix, A, epsilon, m_alpha, herbivore_list, predator_list, species_names = results
    
    # If no initial herbivore densities provided, use the effective baseline.
    if isnothing(H_init)
        H_init = copy(H_i0_eff)
    end
    # Set up predator initial conditions if needed.
    if R > 0
        if isnothing(P_init)
            # Default predator initial densities based on herbivore densities.
            P_init = ones(R) .* (H_init[1:min(R, length(H_init))] ./ 10.0)
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    
    # Build the parameter tuple.
    # NOTE: The new omnivore_dynamics! expects parameters in this order:
    # (S, R, H_i0_eff, g_i, mu_val, O_loss, O_gain, A_star, a_matrix, epsilon, B, m_alpha)
    # Here we use A (returned from the setup) as the predator–predator interaction matrix.
    params = (S, R, H_i0_eff, g_i, mu_val, O_loss, O_gain, A_star, a_matrix, epsilon, A, m_alpha)
    
    # Define and solve the ODE.
    prob = ODEProblem(omnivore_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end
    
    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: sol did not finish properly"
            return nothing
        end
    end
    
    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (600, 500))
        ax = Axis(
            fig[1, 1],
            xlabel="Time",
            ylabel="Biomass",
            title="Dynamics for cell $cell",
            yscale = log ? log10 : identity
        )
        times_combined = sol.t
        
        # Plot herbivore dynamics (indices 1:S) as blue solid lines.
        for i in 1:S
            lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
        end
        # Plot predator dynamics (indices S+1:S+R) as red dashed lines.
        for α in 1:R
            lines!(ax, times_combined, sol[S+α, :], label="P$(α)", linestyle=:dash, color=:red)
        end
    
        display(fig)
    end
    
    # Post-processing: extract final densities.
    H_end = sol.u[end][1:S]
    P_end = (R > 0 ? sol.u[end][S+1:end] : zeros(0))
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    survival_rate = total_surv / total_species
    giHi = sum(g_i .* H_end)
    herbivore_survival_rate = survived_herb / S
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
        NPP                     = localNPP,   # not used
        g_iH_i                  = giHi,
        g_iH_i_over_NPP         = round(giHi, digits=4),
        g_iHobs                 = round(sum(g_i .* H_i0_eff), digits=4),
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
        mu_predation            = mu_pred_val,
        epsilon_val             = eps_val,
        symmetrical_competition = sym_competition
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

@time h_run, sol = omnivore_run(
    1, 
    0.5, 0.012, 1.0, true, 0.1;
    time_end=500.0, 
    do_you_want_params=false,
    do_you_want_sol = true,
    include_predators=true,
    plot=true, sp_removed_name=nothing,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = true,
    o_pred_avg = 0.0,
    log = false
);
