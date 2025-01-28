# Function to attempt feasibility with current species list
function attempt_feasibility(current_sp_nm, local_i, local_j, localNPP, localH0_vector, parameters; many_params = true, sp_removed = false, sp_removed_name = missing)

    cell = findfirst(isequal(CartesianIndex(local_i, local_j)), idx)
    # Track the best result and problematic species
    best_survival_rate  = 0.0
    best_result         = nothing
    
    # Update local_S and local_R based on current species
    current_S, current_R = identify_n_of_herbs_and_preds(current_sp_nm)
    
    if !many_params && length(parameters) > 1
        error("original is false but params is length(param_combinations) > 1")
    end
    if sp_removed && ismissing(sp_removed_name)
        error("sp_removed is true but sp_removed_name is nothing")
    end
    
    for combo in parameters[1:min(end, MAX_ITERS)]
        mu_val, mu_pred_val, eps_val, sym_competition = combo

        # Attempt to set up the community
        results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_pred_val, eps_val, sym_competition;
            localNPP       = localNPP,
            localH0_vector = localH0_vector,
            species_names  = current_sp_nm
        )
        if results === nothing
            continue
        end

        # Destructure the results:
        (S2, R2, H_i0, m_i, g_i, G, M_modified,
         a_matrix, A, epsilon_vector, m_alpha) = (
            results.S, results.R, results.H_i0, results.m_i,
            results.g_i, results.G, results.M_modified,
            results.a_matrix, results.A, results.epsilon_vector,
            results.m_alpha
        )

        if (S2 + R2) == 0 || R2 > length(H_i0)
            continue
        end

        # Build initial conditions.
        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0     = vcat(H_init, P_init)

        params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
        prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
        
        #### THIS APPROACH ELIMINATES THE SOLVER WARNING ####
        logger = SimpleLogger(stderr, Logging.Error)
        sol = with_logger(logger) do
            solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
        #######################################################

        # Skip if integration did not complete to 500 or if nan/inf occurred.
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            continue
        end

        # Evaluate survival
        H_end      = sol[1:S2, end]
        P_end      = sol[S2+1:S2+R2, end]
        survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv    = survived_herb + survived_pred
        total_species = S2 + R2
        survival_rate = total_surv / total_species

        # Update the best result if survival_rate is improved.
        if survival_rate > best_survival_rate
            
            herbivore_survival_rate = survived_herb / S2
            predator_survival_rate  = survived_pred / R2
            H_biomass               = sum(H_end)
            P_biomass               = sum(P_end)
            biomass_at_the_end      = H_biomass + P_biomass
            ratio                   = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
            giHi                    = sum(g_i .* H_end)

            best_survival_rate = survival_rate
            best_result = (
                cell_id                 = cell,
                i                       = local_i,
                j                       = local_j,
                mu                      = mu_val,
                mu_predation            = mu_pred_val,
                epsilon_val             = eps_val,
                symmetrical_competition = sym_competition,
                NPP                     = localNPP,
                g_iH_i                  = giHi,
                g_iH_i_over_NPP         = round(giHi / localNPP, digits=4),
                survived_herbivores     = survived_herb,
                survived_predators      = survived_pred,
                total_survivors         = total_surv,
                total_species           = total_species,
                survival_rate           = survival_rate,
                herbivore_survival_rate = herbivore_survival_rate,
                predator_survival_rate  = predator_survival_rate,
                H_biomass               = H_biomass,
                P_biomass               = P_biomass,
                biomass_at_the_end      = biomass_at_the_end,
                herb_pred_ratio         = ratio,
                many_params             = many_params,
                sp_removed              = sp_removed,
                sp_removed_name         = sp_removed_name
            )
        end

        # Stop early if full survival is achieved.
        if isapprox(survival_rate, 1.0; atol=1e-10)
            ratio_ok = (giHi / localNPP > 0.5) && (giHi / localNPP < 5.0)
            if ratio_ok
            @info "Cell $cell => full survival with (mu=$mu_val, mu_pred=$mu_pred_val, eps=$eps_val). Stopping early."
            return true, best_result  # Feasibility achieved
            end
        end
    end
    # @info "No feasible combination found for cell $cell."
    return false, best_result  # Feasibility not achieved
end