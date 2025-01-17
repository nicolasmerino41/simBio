begin
    ###############################
    # 1) Fixed Parameter Configuration
    ###############################
    mu_val           = 0.189
    mu_predation_val = 0.016
    epsilon_val      = 0.835
    sym_competition  = true

    M_mean               = 0.1
    mean_m_alpha         = 0.1
    m_standard_deviation = 0.0
    h_standard_deviation = 0.0
    artificial_pi        = false
    real_H0              = true

    ###############################
    # 2) Choose a Cell (example: cell 2)
    ###############################
    cell = 4093
    local_i, local_j = idx[cell][1], idx[cell][2]

    ###############################
    # 3) Gather Cell Data from DA and H0_DA
    ###############################
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    ###############################
    # 4) Run Full Community Simulation (Baseline: no species removed)
    ###############################
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

    # Build the full species list (herbivores then predators)
    all_species_names = vcat(herbivore_list_full, predator_list_full)

    # Build initial conditions
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
        # Set values below threshold to zero
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

    # Save full simulation extinction mask for the entire community.
    full_ext_mask = map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, vcat(H_end_full, P_end_full))

    ###############################
    # 5) Prepare DataFrame to Store Removal Experiment Results
    ###############################
    results_df = DataFrame(cell = Int[], sp_removed = String[],
        sp_id_removed = Int[],
        survival_rate = Float64[], total_biomass = Float64[],
        h_biomass = Float64[], p_biomass = Float64[],
        herb_pred_ratio = Float64[],
        herbivore_survival_rate = Float64[],
        predator_survival_rate = Float64[],
        delta_total_biomass = Float64[],
        H_biomass_vector = String[],
        P_biomass_vector = String[],
        H_full_minus_H = String[],
        P_full_minus_P = String[],
        ind_ext_num = String[],
        ind_ext_name = String[],
        survived_herbs = Int[],
        survived_preds = Int[]
    )

    # Save the baseline (no removal) row.
    push!(results_df, (
        cell = cell,
        sp_removed = "none",
        sp_id_removed = 0,
        survival_rate = round(full_survival_rate, digits=3),
        total_biomass = round(full_total_biomass, digits=3),
        h_biomass = round(full_H_biomass, digits=3),
        p_biomass = round(full_P_biomass, digits=3),
        herb_pred_ratio = round(full_herb_pred_ratio, digits=3),
        herbivore_survival_rate = S_full > 0 ? round(full_surv_herb / S_full, digits=3) : 0.0,
        predator_survival_rate = R_full > 0 ? round(full_surv_pred / R_full, digits=3) : 0.0,
        delta_total_biomass = 0.0,
        H_biomass_vector = join(string.(H_end_full), ","),
        P_biomass_vector = join(string.(P_end_full), ","),
        H_full_minus_H = fill(0.0, S_full),
        P_full_minus_P = fill(0.0, R_full),
        ind_ext_num = "none",
        ind_ext_name = "none",
        survived_herbs = full_surv_herb,
        survived_preds = full_surv_pred
    ))

    ###############################
    # 6) Helper Function: Run Simulation and Extract Outputs
    ###############################
    function run_simulation(u0, params, tspan)
        prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
        sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
        S_local = params[1]
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            return (0.0, 0.0, 0.0, 0.0, 0, 0, zeros(Int, S_local+params[2]), zeros(S_local), zeros(params[2]))
        else
            H_end = sol[1:S_local, end]
            P_end = sol[S_local+1:end, end]
            H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
            P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
            H_end_full_minus_H_end = H_end_full .- H_end
            P_end_full_minus_H_end = P_end_full .- P_end
            surv_herb = count(H_end .> EXTINCTION_THRESHOLD)
            surv_pred = count(P_end .> EXTINCTION_THRESHOLD)
            tot_surv = surv_herb + surv_pred
            tot_species = S_local + params[2]
            sr = tot_species > 0 ? tot_surv / tot_species : 0.0
            H_bio = sum(H_end[H_end .> EXTINCTION_THRESHOLD])
            P_bio = sum(P_end[P_end .> EXTINCTION_THRESHOLD])
            # Build a full extinction mask for the community (herbivores + predators)
            ext_mask = vcat(map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, H_end),
                            map(x -> x > EXTINCTION_THRESHOLD ? 0 : 1, P_end))
            return (sr, H_bio + P_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end, H_end_full_minus_H_end, P_end_full_minus_P_end)
        end
    end

    ###############################
    # 7) Removal Experiments for Herbivores and Predators.
    # For each removal experiment, we modify u0_full by setting the given species' initial
    # abundance to 0. We then run the simulation, extract outputs, and compute the collateral extinction
    # list as those species which were alive in full simulation but are now extinct.
    # This is done over the entire community (indices 1:(S_full+R_full)). The only species excluded from
    # the collateral list is the actively removed species.
    ###############################

    # Helper to compute collateral extinctions:
    function collateral_extinctions(active_removed::Int, full_mask::Vector{Int}, current_mask::Vector{Int}, names::Vector{String})
        # active_removed is the full state index of the species actively removed.
        # We return all indices (and corresponding species names) where full_mask == 0 (i.e. species was alive in full run)
        # and current_mask == 1 (species extinct in removal run), excluding active_removed.
        inds = [i for i in 1:length(full_mask) if full_mask[i] == 0 && current_mask[i] == 1 && i != active_removed]
        names_collateral = isempty(inds) ? "none" : join([names[i] for i in inds], ",")
        inds_str = isempty(inds) ? "none" : join(string.(inds), ",")
        return inds_str, names_collateral
    end

    # Prepare an array for species names corresponding to the full state:
    full_state_species_names = vcat(herbivore_list_full, predator_list_full)

    # --- Removal for Herbivores ---
    for (i_sp, sp) in enumerate(herbivore_list_full)
        active_removed = i_sp  # in the full state, herbivore indices are 1:S_full.
        modified_H_init = copy(H_init_full)
        modified_H_init[i_sp] = 0.0  # Remove species sp actively.
        modified_u0 = vcat(modified_H_init, P_init_full)
        sr, tot_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end, H_end_full_minus_H_end, P_end_full_minus_P_end =
            run_simulation(modified_u0, params_full, tspan)
        herb_pred_ratio_mod = H_bio > 0 ? (P_bio / H_bio) : 0.0
        herb_surv_rate_mod = S_full > 0 ? surv_herb / S_full : 0.0
        pred_surv_rate_mod = R_full > 0 ? surv_pred / R_full : 0.0
        delta_total_biomass = full_total_biomass - tot_bio
        H_biomass_vector = join(string.(H_end), ",")
        P_biomass_vector = join(string.(P_end), ",")
        ind_ext_num_str, ind_ext_name_str = collateral_extinctions(active_removed, full_ext_mask, ext_mask, full_state_species_names)
        println("When removing herbivore $sp, collateral extinctions: $ind_ext_name_str; sr = $sr")
        push!(results_df, (
            cell = cell,
            sp_removed = sp,
            sp_id_removed = active_removed,
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
            H_full_minus_H = H_end_full_minus_H_end,
            P_full_minus_P = P_end_full_minus_P_end,
            ind_ext_num = ind_ext_num_str,
            ind_ext_name = ind_ext_name_str,
            survived_herbs = surv_herb,
            survived_preds = surv_pred
        ))
    end

    # --- Removal for Predators ---
    for (i_sp, sp) in enumerate(predator_list_full)
        active_removed = S_full + i_sp  # in the full state, predator indices start at S_full+1.
        modified_P_init = copy(P_init_full)
        modified_P_init[i_sp] = 0.0  # Remove species sp actively.
        modified_u0 = vcat(H_init_full, modified_P_init)
        sr, tot_bio, H_bio, P_bio, surv_herb, surv_pred, ext_mask, H_end, P_end, H_end_full_minus_H_end, P_end_full_minus_P_end =
            run_simulation(modified_u0, params_full, tspan)
        herb_pred_ratio_mod = H_bio > 0 ? (P_bio / H_bio) : 0.0
        herb_surv_rate_mod = S_full > 0 ? surv_herb / S_full : 0.0
        pred_surv_rate_mod = R_full > 0 ? surv_pred / R_full : 0.0
        delta_total_biomass = full_total_biomass - tot_bio
        H_biomass_vector = join(string.(H_end), ",")
        P_biomass_vector = join(string.(P_end), ",")
        ind_ext_num_str, ind_ext_name_str = collateral_extinctions(active_removed, full_ext_mask, ext_mask, full_state_species_names)
        println("When removing predator $sp, collateral extinctions: $ind_ext_name_str; sr = $sr")
        push!(results_df, (
            cell = cell,
            sp_removed = sp,
            sp_id_removed = active_removed,
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
            ind_ext_name = ind_ext_name_str,
            survived_herbs = surv_herb,
            survived_preds = surv_pred
        ))
    end

    ###############################
    # 8) Final Outputs: the DataFrame results_df now holds all removal experiment results.
    ###############################
    @info "Removal experiments finished for cell $cell."
    @info "Final removal experiment results for cell $cell:"
    println(results_df)
end


begin
    
    # Assume that results_df is from your removal experiments and
    # species_metrics is from your compute_food_web_metrics(cell) call.
    # Here, results_df contains a column "sp_removed" with species names (or "none" for no removal)
    # and species_metrics contains a column "species" with species names.

    species_net_metrics = compute_food_web_metrics(cell).species_metrics

    # Filter out the full community baseline row (if sp_removed == "none")
    removal_df = filter(row -> row.sp_removed != "none", results_df)

    # Perform an inner join on the species name.
    merged_df = innerjoin(removal_df, species_net_metrics, on = [:sp_removed => :species])

    # For the full community baseline (sp_removed == "none"), you may wish to add a row with NaN
    # for the network metrics.
    baseline_df = filter(row -> row.sp_removed == "none", results_df)
    baseline_df[!, :indegree] .= NaN
    baseline_df[!, :outdegree] .= NaN
    baseline_df[!, :total_degree] .= NaN
    baseline_df[!, :betweenness] .= NaN
    baseline_df[!, :closeness] .= NaN
    baseline_df[!, :clustering] .= NaN

    # Now, you can combine the merged removal results with the baseline row.
    merged_df = vcat(baseline_df, merged_df)

    @info "Merged removal and network metrics:"
    display(merged_df)
end

haches = merged_df[end, 12]
typeof(haches)
haches = map(x -> parse(Float64, x), split(haches, ","))
pees = merged_df[end, 13]
pees = map(x -> parse(Float64, x), split(pees, ","))
haches
pees

any(iszero, haches)
any(iszero, pees)
count(haches .> 0)
count(pees .> 0)
println("pees: ", pees)