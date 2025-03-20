function new_single_run_with_plot(
    cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
    plot=false, sp_removed_name=nothing,
    artificial_pi=false, NPP=nothing,
    alpha = 0.25
)
    
    # Placeholder for results DataFrame
    AAAA = DataFrame()

    # @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]

    # Extract species names from the cell data.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey = check_predator_has_prey(sp_nm)
    # println("here S = $local_S, R = $local_R")
    if !predator_has_prey[1]
        local_R -= predator_has_prey[2]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
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
        alpha          = alpha
    )

    if isnothing(results)
        @error "Error: results === nothing"
        return nothing, nothing
    end

    # Destructure the returned NamedTuple.
    S2          = results.S
    R2          = results.R
    H_i0        = results.H_i0
    m_i         = results.m_i
    g_i         = results.g_i
    beta        = results.beta
    M_mod       = results.M_modified    # Modified competition matrix (S×S)
    A_star      = results.A_star        # Nondimensional predation rates (S×R)
    a_matrix    = results.a_matrix      # Herbivore–predator interaction matrix (S×R)
    A           = results.A             # Predator interaction matrix (R×R)
    m_alpha     = results.m_alpha

    # Build predator attack matrix (R×S) by transposing a_matrix.
    A_pred = transpose(a_matrix)

    # Set predator overpopulation thresholds.
    # (Here we assume P0 = m_alpha, which implies a self-regulation coefficient of 1.)
    P0 = m_alpha

    # println("here S = $S2, R = $R2")
    if (S2 + R2) == 0 || R2 > length(H_i0)
        @error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
        return nothing, nothing
    end

    # Build initial conditions.
    # Herbivore initial conditions are set to the empirical abundances (H_i0).
    # Predator initial conditions are initialized as a fraction of the herbivore abundances.
    H_init = H_i0
    P_init = H_init[1:R2] ./ 10.0
    u0 = vcat(H_init, P_init)

    # Package parameters for new_dynamics!.
    # new_dynamics! expects the following parameter ordering:
    # (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha)
    params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
    prob = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params)

    # Solve the ODE problem with strict tolerances.
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end

    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        @error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
        return nothing, nothing
    end

    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(resolution=(1200, 600))
        ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass", title="Dynamics for cell $cell")
        times_combined = sol.t
        @info "1"
        # Plot herbivore dynamics (indices 1:S2) as blue solid lines.
        for i in 1:S2
            lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
        end
        @info "2"
        # Plot predator dynamics (indices S2+1:S2+R2) as red dashed lines.
        for α in 1:R2
            lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
        end
        @info "3"
    
        display(fig)
    end
    @info "4"
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
    @info "5"
    append!(AAAA, single_run_results)
    @info "The survival rate is $(round(survival_rate, digits=4))"
    
    # Return both the DataFrame with summary results and the Makie figure.
    return AAAA, params
end

# Example call:
AAAA, params = new_single_run_with_plot(
    1, 0.0001, 0.0, 0.0, true, 0.3;
    plot=false, sp_removed_name=nothing, 
    NPP=nothing, artificial_pi=false,
    alpha = 0.25
)

# S, R, Hi0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha = params
Threads.@threads for mu_pred in 0.0:0.01:0.2
        for mu in 0.1:0.1:1.0, eps in 0.01:0.1:1.0, m_alpha in 0.0:0.1:1.0, alpha in 0.0:0.1:1.0
        AAAA, params = new_single_run_with_plot(1, mu, mu_pred, eps, true, m_alpha; plot=false, sp_removed_name=nothing, NPP=nothing, artificial_pi=true, alpha = alpha);
        if !isnothing(AAAA) && AAAA.survived_predators[1] > 0
            println("mu = $mu, mu_pred = $mu_pred, eps = $eps, m_alpha = $m_alpha, alpha = $alpha", )
            break
        end
    end
end