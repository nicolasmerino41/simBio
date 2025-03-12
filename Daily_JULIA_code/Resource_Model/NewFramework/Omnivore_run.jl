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
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    o_pred_avg::Float64 = 1.0,
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
    results = o_new_attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
        localNPP = localNPP,
        species_names = sp_nm,
        o_pred_avg = o_pred_avg
    )
    if isnothing(results)
        @error "Error: results === nothing"
        return nothing
    end
    
    S2         = results.S
    R2         = results.R
    H_i0       = results.H_i0       # effective herbivore abundance (with net herbivore-herbivore interactions)
    m_i        = results.m_i
    g_i        = results.g_i
    beta       = results.beta      # zeros
    O_loss     = results.O_loss
    O_gain     = results.O_gain
    A_loss     = results.A_star    # predator-induced loss matrix
    a_matrix   = results.a_matrix
    A          = results.A
    m_alpha    = results.m_alpha
    h          = results.h
    
    if isnothing(H_init)
        H_init = copy(H_i0)
    end
    if R2 > 0
        if isnothing(P_init)
            P_init = ones(R2) .* (H_init[1:min(R2, length(H_init))] ./ 10.0)
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    
    if R2 > 0
        params = (S2, R2, H_i0, m_i, g_i, beta, O_loss, O_gain, A_loss, A, m_alpha)
    else
        params = (S2, R2, H_i0, m_i, g_i, beta, zeros(S2,S2), zeros(S2,S2), zeros(S2,0), zeros(0,0), zeros(0))
    end
    
    prob = ODEProblem(o_omnivore_dynamics!, u0, (0.0, time_end), params)
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
    
    H_end = sol.u[end][1:S2]
    P_end = (R2 > 0 ? sol.u[end][S2+1:end] : zeros(0))
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species
    giHi = sum(g_i .* H_end)
    herbivore_survival_rate = survived_herb / S2
    predator_survival_rate = (R2 > 0 ? survived_pred / R2 : 0.0)
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
        g_iHobs                 = round(sum(g_i .* H_i0), digits=4),
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

@time h_run = omnivore_run(
    1, 
    0.0, 0.0, 0.0, true, 1.0;
    plot=true,
    o_pred_avg = 0.0
)