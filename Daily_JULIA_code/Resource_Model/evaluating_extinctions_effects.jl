begin
    ################################################################################
    # run_removals_fixed_configuration_with_collateral_extinctions.jl
    #
    # This script runs a full community simulation (the baseline) for a chosen cell 
    # with a fixed parameter configuration. Then, for each species (first for herbivores, 
    # then for predators), its initial abundance is set to 0; the ODE is re-integrated; 
    # outputs are computed; and collateral extinctions are defined as those species that were 
    # alive in the full simulation but become extinct in the removal run. The extra outputs 
    # include the delta in total biomass, the final biomass vectors (H and P) as a string, 
    # and a record (by indices and names) of collateral extinctions.
    ################################################################################

    # --------------------------
    # 1) Fixed parameter configuration
    # --------------------------
    mu_val           = 0.367
    mu_predation_val = 0.019
    epsilon_val      = 0.982
    sym_competition  = true

    M_mean               = 0.1
    mean_m_alpha         = 0.1
    m_standard_deviation = 0.0
    h_standard_deviation = 0.0
    artificial_pi        = false
    real_H0              = true

    # --------------------------
    # 2) Choose cell (example: cell 2)
    # --------------------------
    cell = 1
    local_i, local_j = idx[cell][1], idx[cell][2]

    # --------------------------
    # 3) Gather cell data from DA and H0_DA
    # --------------------------
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # --------------------------
    # 4) Run the full community simulation (baseline; no species removed)
    # --------------------------
    full_results = setup_community_from_cell(local_i, local_j;
        NPP = localNPP,
        M_mean = M_mean,
        mu = mu_val,
        symmetrical_competition = sym_competition,
        mean_m_alpha = mean_m_alpha,
        epsilon_val = epsilon_val,
        mu_predation = mu_predation_val,
        iberian_interact_NA = iberian_interact_NA,
        species_dict = species_dict,
        m_standard_deviation = m_standard_deviation,
        h_standard_deviation = h_standard_deviation,
        artificial_pi = artificial_pi,
        real_H0 = real_H0,
        H0_vector = localH0_vector
    )
    if full_results === nothing
        error("Full community setup failed for cell $cell.")
    end
    (S_full, R_full, species_names_full, herbivore_list_full, predator_list_full,
     H_i0_full, m_i_full, p_vec_full, x_final_full, g_i_full, localHatH_full,
     G_full, M_modified_full, a_matrix_full, A_full, epsilon_vector_full, m_alpha_full) = full_results

    # Define full species list (herbivores then predators)
    all_species_names = vcat(herbivore_list_full, predator_list_full)

    # Build initial conditions and run the full simulation
    H_init_full = H_i0_full
    if R_full > length(H_init_full)
        error("In cell $cell: More predators than herbivore data available.")
    end
    P_init_full = H_init_full[1:R_full] ./ 10.0
    u0_full = vcat(H_init_full, P_init_full)
    params_full = (S_full, R_full, H_i0_full, m_i_full, g_i_full, G_full,
                   M_modified_full, a_matrix_full, A_full, epsilon_vector_full, m_alpha_full)
    tspan = (0.0, 500.0)
    prob_full = ODEProblem(ecosystem_dynamics!, u0_full, tspan, params_full)
    sol_full = solve(prob_full, Tsit5(); reltol=1e-6, abstol=1e-6)

    if sol_full.t[end] < 500.0 || any(isnan, sol_full.u[end]) || any(isinf, sol_full.u[end])
        full_survival_rate = 0.0
        full_H_biomass = 0.0
        full_P_biomass = 0.0
        full_surv_herb = 0
        full_surv_pred = 0
    else
        H_end_full = sol_full[1:S_full, end]
        P_end_full = sol_full[S_full+1:S_full+R_full, end]
        H_end_full[H_end_full .< EXTINCTION_THRESHOLD] .= 0.0
        P_end_full[P_end_full .< EXTINCTION_THRESHOLD] .= 0.0
        full_surv_herb = count(H_end_full .> EXTINCTION_THRESHOLD)
        full_surv_pred = count(P_end_full .> EXTINCTION_THRESHOLD)
        full_total_surv = full_surv_herb + full_surv_pred
        full_total_species = S_full + R_full
        full_survival_rate = full_total_species > 0 ? full_total_surv / full_total_species : 0.0
        full_H_biomass = sum(H_end_full[H_end_full .> EXTINCTION_THRESHOLD])
        full_P_biomass = sum(P_end_full[P_end_full .> EXTINCTION_THRESHOLD])
    end
    full_herb_pred_ratio = full_H_biomass > 0 ? full_P_biomass / full_H_biomass : 0.0
    full_total_biomass = full_H_biomass + full_P_biomass

    # Save full community extinction masks (only for species that were alive at baseline)
    full_ext_mask = map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, vcat(H_end_full, P_end_full))
    
    # --------------------------
    # 5) Prepare DataFrame for removal experiments
    # --------------------------
    AAA = DataFrame(cell = Int[], sp_removed = String[],
        sp_id_removed = Int[],
        survival_rate = Float64[], total_biomass = Float64[],
        h_biomass = Float64[], p_biomass = Float64[],
        herb_pred_ratio = Float64[],
        herbivore_survival_rate = Float64[],
        predator_survival_rate = Float64[],
        delta_total_biomass = Float64[],
        H_biomass_vector = String[],
        P_biomass_vector = String[],
        ind_ext_num = String[],
        ind_ext_name = String[],
        survived_herbs = Int[],
        survived_preds = Int[]
    )

    # Save the full simulation results as baseline (sp_removed = "none")
    push!(AAA, (
        cell = cell,
        sp_removed = "none",
        sp_id_removed = 0,
        survival_rate = round(full_survival_rate, digits=3),
        total_biomass = round(full_total_biomass, digits=3),
        h_biomass = round(full_H_biomass, digits=3),
        p_biomass = round(full_P_biomass, digits=3),
        herb_pred_ratio =  round(full_herb_pred_ratio, digits=3),
        herbivore_survival_rate = S_full > 0 ? round(full_surv_herb / S_full, digits=3) : 0.0,
        predator_survival_rate = R_full > 0 ? round(full_surv_pred / R_full, digits=3) : 0.0,
        delta_total_biomass = 0.0,
        H_biomass_vector = join(string.(H_end_full), ","),
        P_biomass_vector = join(string.(P_end_full), ","),
        ind_ext_num = "none",
        ind_ext_name = "none",
        survived_herbs = full_surv_herb,
        survived_preds = full_surv_pred
    ))

    # --------------------------
    # 6) Helper function to run simulation and extract outputs.
    # --------------------------
    function run_simulation(u0, params, tspan)
        prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
        sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
        S_local = params[1]
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            return (0.0, 0.0, 0.0, 0.0, 0, 0, zeros(Int, S_local), zeros(S_local), zeros(S_local))
        else
            H_end = sol[1:S_local, end]
            P_end = sol[S_local+1:end, end]
            H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
            P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
            surv_herb = count(H_end .> EXTINCTION_THRESHOLD)
            surv_pred = count(P_end .> EXTINCTION_THRESHOLD)
            tot_surv = surv_herb + surv_pred
            tot_species = params[1] + params[2]
            sr = tot_species > 0 ? tot_surv / tot_species : 0.0
            H_bio = sum(H_end[H_end .> EXTINCTION_THRESHOLD])
            P_bio = sum(P_end[P_end .> EXTINCTION_THRESHOLD])
            ext_mask = vcat(map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, H_end),
                            map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, P_end))
            return (sr, H_bio + P_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end)
        end
    end

    # --------------------------
    # 7) Removal experiments for herbivores.
    # --------------------------
    for (i_sp, sp) in enumerate(herbivore_list_full)
        modified_H_init = copy(H_init_full)
        modified_H_init[i_sp] = 0.0  # Actively remove species sp.
        modified_u0 = vcat(modified_H_init, P_init_full)
        (sr, tot_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end) =
            run_simulation(modified_u0, params_full, tspan)
        herb_pred_ratio_mod = H_bio > 0 ? (P_bio / H_bio) : 0.0
        herb_surv_rate_mod = S_full > 0 ? surv_herb / S_full : 0.0
        pred_surv_rate_mod = R_full > 0 ? surv_pred / R_full : 0.0
        delta_total_biomass = full_total_biomass - tot_bio
        H_biomass_vector = join(string.(H_end), ",")
        P_biomass_vector = join(string.(P_end), ",")
       
        # Those species that (i) were alive in the full simulation (i.e. full_ext_mask_herb[i] == 0)
        # but (ii) become extinct (ext_mask[i] == 1) in the removal run, excluding the species actively removed.
        collateral_indices = [i for i in 1:S_full+R_full if full_ext_mask[i] == 0 && ext_mask[i] == 1 && i != i_sp]
        collateral_names = isempty(collateral_indices) ? "none" : join([sp_nm[i] for i in collateral_indices], ",")
        ind_ext_num_str = isempty(collateral_indices) ? "none" : join(string.(collateral_indices), ",")

        push!(AAA, (
            cell = cell,
            sp_removed = sp,
            sp_id_removed = i_sp,
            survival_rate = round(sr, digits=3),
            total_biomass = round(tot_bio, digits=3),
            delta_total_biomass = round(delta_total_biomass, digits=3),
            h_biomass = round(H_bio, digits=3),
            p_biomass = round(P_bio, digits=3),
            herb_pred_ratio = round(herb_pred_ratio_mod, digits=3),
            herbivore_survival_rate = round(herb_surv_rate_mod, digits=3),
            predator_survival_rate = round(pred_surv_rate_mod, digits=3),
            H_biomass_vector = H_biomass_vector,
            P_biomass_vector = P_biomass_vector,
            ind_ext_num = ind_ext_num_str,
            ind_ext_name = collateral_names,
            survived_herbs = surv_herb,
            survived_preds = surv_pred
        ))
    end

    # --------------------------
    # 8) Removal experiments for predators.
    # --------------------------
    for (i_sp, sp) in enumerate(predator_list_full)
        sp_number = S_full + i_sp
        modified_P_init = copy(P_init_full)
        modified_P_init[i_sp] = 0.0  # Remove species sp actively.
        modified_u0 = vcat(H_init_full, modified_P_init)
        (sr, tot_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end) =
            run_simulation(modified_u0, params_full, tspan)
        herb_pred_ratio_mod = H_bio > 0 ? (P_bio / H_bio) : 0.0
        herb_surv_rate_mod = S_full > 0 ? surv_herb / S_full : 0.0
        pred_surv_rate_mod = R_full > 0 ? surv_pred / R_full : 0.0
        delta_total_biomass = full_total_biomass - tot_bio
        H_biomass_vector = join(string.(H_end), ",")
        P_biomass_vector = join(string.(P_end), ",")
        # FALSE For predators, we consider collateral extinctions only in the predator compartment.
        # In the full community, the predator indices are S_full+1:S_full+R_full.
        collateral_indices = [i for i in 1:S_full+R_full if full_ext_mask[i] == 0 && ext_mask[i] == 1 && i != i_sp]
        collateral_names = isempty(collateral_indices) ? "none" : join([sp_nm[i] for i in collateral_indices], ",")
        ind_ext_num_str = isempty(collateral_indices) ? "none" : join(string.(collateral_indices), ",")
        
        push!(AAA, (
            cell = cell,
            sp_removed = sp,
            sp_id_removed = sp_number,
            survival_rate = round(sr, digits=3),
            total_biomass = round(tot_bio, digits=3),
            delta_total_biomass = round(delta_total_biomass, digits=3),
            h_biomass = round(H_bio, digits=3),
            p_biomass = round(P_bio, digits=3),
            herb_pred_ratio = round(herb_pred_ratio_mod, digits=3),
            herbivore_survival_rate = round(herb_surv_rate_mod, digits=3),
            predator_survival_rate = round(pred_surv_rate_mod, digits=3),
            H_biomass_vector = H_biomass_vector,
            P_biomass_vector = P_biomass_vector,
            ind_ext_num = ind_ext_num_str,
            ind_ext_name = collateral_names,
            survived_herbs = surv_herb,
            survived_preds = surv_pred
        ))
    end

    # --------------------------
    # 9) Final outputs: the DataFrame results_df now holds all removal experiment results.
    # --------------------------
    @info "Removal experiments finished for cell $cell."
    # @info "Results:"
    # println(AAA)
end
