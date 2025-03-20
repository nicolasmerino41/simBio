function herbivore_run(
    cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
    time_end=500.0, 
    do_you_want_params=false,
    do_you_want_sol = false,
    include_predators=false,
    plot=false, sp_removed_name=nothing,
    artificial_pi=false, NPP=nothing,
    alpha = 0.25,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    hollingII = false,
    h = 0.1,
    log = false
)
    
    # Placeholder for results DataFrame
    AAAA = DataFrame()

    # @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]

    # Extract species names from the cell data.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
     
    if !include_predators
        local_R  = 0
    
        filter!(name -> !(name in predator_names), sp_nm)
    end
    if !isnothing(sp_removed_name)
        filter!(name -> (name != sp_removed_name), sp_nm)
    end

    # Obtain local NPP and empirical herbivore abundances.
    localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
    # localH0_vector = nothing

    if !isnothing(NPP)
        localNPP = NPP
    end

    # Set up community parameters using the updated parametrisation.
    results = new_attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
        localNPP       = localNPP,
        # localH0_vector = localH0_vector,
        species_names  = sp_nm,
        artificial_pi  = artificial_pi,
        alpha          = alpha,
        hollingII      = hollingII,
        h              = h
    )

    if isnothing(results)
        @error "Error: results === nothing"
        if do_you_want_params && do_you_want_sol
            return nothing, nothing, nothing
        elseif do_you_want_params || do_you_want_sol
            return nothing, nothing
        else
            return nothing
        end
    end

    # Destructure the returned NamedTuple.
    S2          = results.S
    R2          = results.R
    H_i0        = results.H_i0
    m_i         = results.m_i
    g_i         = results.g_i
    beta        = results.beta
    M_mod       = results.M_modified    
    A_star      = results.A_star        
    a_matrix    = results.a_matrix      
    A           = results.A             
    m_alpha     = results.m_alpha
    h           = results.h

    A_pred = transpose(a_matrix)

    # Set predator overpopulation thresholds.
    # (Here we assume P0 = m_alpha, which implies a self-regulation coefficient of 1.)
    P0 = m_alpha

    # println("here S = $S2, R = $R2")
    if (S2 + R2) == 0 || R2 > length(H_i0)
        @error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
        if do_you_want_params && do_you_want_sol
            return nothing, params, nothing
        elseif do_you_want_params || do_you_want_sol
            return nothing, do_you_want_params ? params : nothing
        else
            return nothing
        end
    end

    # Build initial conditions.
    # Herbivore initial conditions are set to the empirical abundances (H_i0).
    # Predator initial conditions are initialized as a fraction of the herbivore abundances.
    if isnothing(H_init)
        H_init = fill(1.0, length(H_i0))#H_i0
    end
    if isnothing(P_init)
        P_init = H_init[1:R2] ./ 10.0
    end
    u0 = vcat(H_init, P_init)

    # Package parameters for new_dynamics!.
    # new_dynamics! expects the following parameter ordering:
    # (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha)
    if !hollingII
        params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
    else
        params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha, h)
    end
    # cb_no_trigger, cb_trigger = build_callbacks(S2, R2, EXTINCTION_THRESHOLD, T_ext, 1)
    prob = ODEProblem(new_dynamics!, u0, (0.0, time_end), params)

    # Solve the ODE problem with strict tolerances.
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); callback=cb_no_trigger, abstol=1e-8, reltol=1e-6)
    end

    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
            if do_you_want_params && do_you_want_sol
                return nothing, params, sol
            elseif do_you_want_params || do_you_want_sol
                return nothing, do_you_want_params ? params : sol
            else
                return nothing
            end
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
        
        # @info "1"
        # Plot herbivore dynamics (indices 1:S2) as blue solid lines.
        for i in 1:S2
            lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
        end
        # @info "2"
        # Plot predator dynamics (indices S2+1:S2+R2) as red dashed lines.
        for α in 1:R2
            lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
        end
        # @info "3"
    
        display(fig)
    end
    # @info "4"
    # --- Evaluate Outputs ---
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
    P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
    Hobs = H_i0
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species
    giHi = sum(g_i .* H_end)
    herbivore_survival_rate = survived_herb / S2
    predator_survival_rate = survived_pred / R2
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
    # println("g_i = ", g_i)
    # Store the results in a DataFrame row.
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
        NPP                     = localNPP,
        g_iH_i                  = giHi,
        g_iH_i_over_NPP         = round(giHi / localNPP, digits=4),
        g_iHobs                 = round(sum(g_i .* Hobs), digits=4),
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
    # @info "5"
    append!(AAAA, single_run_results)
    # @info "The survival rate is $(round(survival_rate, digits=4))"
    
    # Return both the DataFrame with summary results and the Makie figure.
    if do_you_want_params && do_you_want_sol
        return AAAA, params, sol
    elseif do_you_want_params || do_you_want_sol
        return AAAA, do_you_want_params ? params : sol
    else
        return AAAA
    end
end

# function run_right_away(
#     params;
#     time_end=500.0,
#     plot=false,
#     H_init = nothing,
#     P_init = nothing,
#     hollingII = false,
#     h = 0.1
# )
    
#     # Build initial conditions.
#     # Herbivore initial conditions are set to the empirical abundances (H_i0).
#     # Predator initial conditions are initialized as a fraction of the herbivore abundances.
#     if isnothing(H_init)
#         H_init = fill(1.0, length(H_i0))#H_i0
#     end
#     if isnothing(P_init)
#         P_init = H_init[1:R2] ./ 10.0
#     end
#     u0 = vcat(H_init, P_init)

#     # Package parameters for new_dynamics!.
#     # new_dynamics! expects the following parameter ordering:
#     # (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha)
#     if !hollingII
#         S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha = params
#     else
#         S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha, h = params
#     end

#     prob = ODEProblem(new_dynamics!, u0, (0.0, time_end), params)

#     # Solve the ODE problem with strict tolerances.
#     logger = SimpleLogger(stderr, Logging.Error)
#     sol = with_logger(logger) do
#         solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
#     end

#     # if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
#     #     @error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
#     #     return nothing
#     # end

#     # --- Plotting the Dynamics using Makie ---
#     if plot    
#         fig = Figure(; size = (600, 500))
#         ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass")
#         times_combined = sol.t
#         # @info "1"
#         # Plot herbivore dynamics (indices 1:S2) as blue solid lines.
#         for i in 1:S2
#             lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
#         end
#         # @info "2"
#         # Plot predator dynamics (indices S2+1:S2+R2) as red dashed lines.
#         for α in 1:R2
#             lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
#         end
#         # @info "3"
    
#         display(fig)
#     end
#     # @info "4"
#     # --- Evaluate Outputs ---
#     H_end = sol[1:S2, end]
#     P_end = sol[S2+1:S2+R2, end]
#     H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
#     P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
#     Hobs = H_i0
#     survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
#     survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
#     total_surv = survived_herb + survived_pred
#     total_species = S2 + R2
#     survival_rate = total_surv / total_species
#     giHi = sum(g_i .* H_end)
#     herbivore_survival_rate = survived_herb / S2
#     predator_survival_rate = survived_pred / R2
#     H_biomass = sum(H_end)
#     P_biomass = sum(P_end)
#     biomass_at_the_end = H_biomass + P_biomass
#     ratio = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
#     # println("g_i = ", g_i)
#     # Store the results in a DataFrame row.
#     single_run_results = DataFrame(
#         # cell_id                 = cell,
#         # i                       = local_i,
#         # j                       = local_j,
#         survival_rate           = survival_rate,
#         # NPP                     = localNPP,
#         g_iH_i                  = giHi,
#         # g_iH_i_over_NPP         = round(giHi / localNPP, digits=4),
#         g_iHobs                 = round(sum(g_i .* Hobs), digits=4),
#         survived_herbivores     = survived_herb,
#         survived_predators      = survived_pred,
#         H_biomass               = H_biomass,
#         P_biomass               = P_biomass,
#         biomass_at_the_end      = biomass_at_the_end,
#         herb_pred_ratio         = ratio,
#         total_survivors         = total_surv,
#         total_species           = total_species,
#         herbivore_survival_rate = herbivore_survival_rate,
#         predator_survival_rate  = predator_survival_rate,
#         # mu                      = mu_val,
#         # mu_predation            = mu_pred_val,
#         # epsilon_val             = eps_val,
#         # symmetrical_competition = sym_competition,
#         H_vector                = [H_end], 
#         P_vector                = [P_end]
#     )
#     # @info "5"
#     # append!(AAAA, single_run_results)
#     @info "The survival rate is $(round(survival_rate, digits=4))"
    
#     # Return both the summary results and the Makie figure.
#     return single_run_results
# end

# # Example call:
# AAAA = herbivore_run(
#     1, 
#     0.01, 0.01, 1.0, true, 0.;
#     include_predators=true,
#     time_end=500.0,
#     do_you_want_params=false,
#     do_you_want_sol = false,
#     plot=true, sp_removed_name=nothing, 
#     NPP=nothing, artificial_pi=false,
#     ignore_inf_error = true,
#     alpha = 0.25,
#     hollingII = true,
#     h = 0.1
# );

# # S, R, Hi0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha = params
# Threads.@threads for mu_pred in 0.0:0.01:0.2
#         for mu in 0.1:0.1:1.0, eps in 0.01:0.1:1.0, m_alpha in 0.0:0.1:1.0, alpha in 0.0:0.1:1.0
#         AAAA, params = new_single_run_with_plot(1, mu, mu_pred, eps, true, m_alpha; plot=false, sp_removed_name=nothing, NPP=nothing, artificial_pi=true, alpha = alpha);
#         if !isnothing(AAAA) && AAAA.survived_predators[1] > 0
#             println("mu = $mu, mu_pred = $mu_pred, eps = $eps, m_alpha = $m_alpha, alpha = $alpha", )
#             break
#         end
#     end
# end