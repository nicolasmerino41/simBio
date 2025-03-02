# Read best parameters and filter full survival cells
best_params = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/best_params_5950_cells_not_rounded_254_1950iter.csv", DataFrame)
best_params = dropmissing(best_params)
new_best_params = similar(best_params, 0)
unique_cells = unique(best_params.cell_id)
for cell in unique_cells
    cell_df = best_params[best_params.cell_id .== cell, :]
    new_best_params = vcat(new_best_params, cell_df[findmax(cell_df.survival_rate)[2], :])
end

full_survival_df = filter(row -> row.survival_rate == 1.0, best_params)
# best_params = filter(row -> row.cell_id != 5919, best_params)
best_params = filter(row -> row.g_iH_i_over_NPP < 10.0, best_params)
average_survival_rate = mean(best_params.survival_rate)
all_results_list = []
# Iterate over all full survival cells
Threads.@threads for i in 1:nrow(best_params)
    
    cell = best_params[i, :cell_id]
    target_sv = best_params[i, :survival_rate]
    
    begin
        ###############################
        # 1) Fixed Parameter Configuration
        ###############################
        mu_val           = best_params[i, :mu]
        mu_predation_val = best_params[i, :mu_predation]
        epsilon_val      = best_params[i, :epsilon_val]
        sym_competition  = true

        M_mean               = 0.1
        mean_m_alpha         = 0.1
        m_standard_deviation = 0.0
        h_standard_deviation = 0.0
        artificial_pi        = false
        real_H0              = true

        ###############################
        # 2) Choose a Cell
        ###############################
        local_i, local_j = idx[cell][1], idx[cell][2]

        ###############################
        # 3) Gather Cell Data from DA and H0_DA
        ###############################
        sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
        local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

        predator_has_prey = check_predator_has_prey(sp_nm)

        if !predator_has_prey[1]
            local_R -= predator_has_prey[2]
            filter!(name -> !(name in predator_has_prey[3]), sp_nm)
            @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
            species_names = sp_nm
        else
            species_names = sp_nm
        end

        localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
        localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

        # cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)
        
        ###############################
        # 4) Run Full Community Simulation (Baseline: no species removed)
        ###############################
        # full_results = setup_community_from_cell(local_i, local_j;
        #     NPP = localNPP,
        #     M_mean = M_mean,
        #     mu = mu_val,
        #     symmetrical_competition = sym_competition,
        #     mean_m_alpha = mean_m_alpha,
        #     epsilon_val = epsilon_val,
        #     mu_predation = mu_predation_val,
        #     iberian_interact_NA = iberian_interact_NA,
        #     species_dict = species_dict,
        #     m_standard_deviation = m_standard_deviation,
        #     h_standard_deviation = h_standard_deviation,
        #     artificial_pi = artificial_pi,
        #     real_H0 = real_H0,
        #     H0_vector = localH0_vector
        # )

        # Attempt setup
        full_results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_predation_val, epsilon_val, sym_competition;
            localNPP      = localNPP,
            localH0_vector= localH0_vector,
            species_names = species_names
        )
        if full_results === nothing
            continue
        end

        # Unpack the full simulation results
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
        params_full = (
            S_full, R_full, H_i0_full, m_i_full, g_i_full, G_full,
            M_modified_full, a_matrix_full, A_full, epsilon_vector_full, m_alpha_full
        )
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

    if abs(full_survival_rate - target_sv) > 0.05
        @error "Cell $cell: Full survival rate is $(full_survival_rate), but target is $target_sv."
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
        H_biomass_vector = Vector[],
        P_biomass_vector = Vector[],
        H_full_minus_H = Vector[],
        P_full_minus_P = Vector[],
        ind_ext_num = String[],
        ind_ext_name = Vector[],
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
        H_biomass_vector = H_end_full,
        P_biomass_vector = P_end_full,
        H_full_minus_H = fill(0.0, S_full),
        P_full_minus_P = fill(0.0, R_full),
        ind_ext_num = "none",
        ind_ext_name = ["none"],
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
            P_end_full_minus_P_end = P_end_full .- P_end
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
        names_collateral = isempty(inds) ? ["none"] : [names[i] for i in inds]
        inds_str = isempty(inds) ? "none" : join(string.(inds), ",")
        return inds_str, names_collateral
    end

    # Prepare an array for species names corresponding to the full state:
    full_state_species_names = vcat(herbivore_list_full, predator_list_full)

    # --- Removal for Herbivores ---
    for (i_sp, sp) in enumerate(herbivore_list_full)
        # println("Removing herbivore $sp")
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
        H_biomass_vector = H_end
        P_biomass_vector = P_end
        ind_ext_num_str, ind_ext_name_str = collateral_extinctions(active_removed, full_ext_mask, ext_mask, full_state_species_names)
        # println("When removing herbivore $sp, collateral extinctions: $ind_ext_name_str; sr = $sr")
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
        H_biomass_vector = H_end
        P_biomass_vector = P_end
        ind_ext_num_str, ind_ext_name_str = collateral_extinctions(active_removed, full_ext_mask, ext_mask, full_state_species_names)
        # println("When removing predator $sp, collateral extinctions: $ind_ext_name_str; sr = $sr")
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
        
        ###############################
        # 5) Compute Network Metrics and Merge
        ###############################
        species_net_metrics = compute_food_web_metrics(cell).species_metrics

        # Filter out the full community baseline row (if sp_removed == "none")
        removal_df = filter(row -> row.sp_removed != ["none"], results_df)

        # Perform an inner join on the species name
        merged_df = innerjoin(removal_df, species_net_metrics, on = [:sp_removed => :species])

        # For the full community baseline (sp_removed == "none"), add a row with NaN for the network metrics
        baseline_df = filter(row -> row.sp_removed == "none", results_df)
        baseline_df[!, :indegree] .= NaN
        baseline_df[!, :outdegree] .= NaN
        baseline_df[!, :total_degree] .= NaN
        baseline_df[!, :betweenness] .= NaN
        baseline_df[!, :closeness] .= NaN
        baseline_df[!, :clustering] .= NaN

        # Combine the merged removal results with the baseline row
        merged_df = vcat(baseline_df, merged_df)

        # Add the result for this cell to the list
        push!(all_results_list, merged_df)

        @info "Merged removal and network metrics for cell $cell" #. sv was $full_survival_rate"
    end
end

# serialize("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/all_results_list_from_best_params_5950_cells_254_1950iter_newherbivores.jls", all_results_list)
all_results_list = deserialize("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/all_results_list_from_best_params_5950_cells_254_1950iter.jls")

# Checking sr matches best_params and all_results_list with full_species by matching cell_id
comparison_df = DataFrame(best_params_sr = [], cell_from_best = [], cell_from_all = [], all_results_list_sr = [], diff = [])
for bp_row in eachrow(best_params)
    cell_id = bp_row[:cell_id]
    all_results_row = findfirst(x -> x[1, :cell] == cell_id, all_results_list)
    
    if all_results_row !== nothing
        best_params_sr = bp_row[:survival_rate]
        all_results_list_sr = all_results_list[all_results_row][1, :survival_rate]
        diff = best_params_sr - all_results_list_sr
        push!(comparison_df, (best_params_sr, cell_id, all_results_list[all_results_row][1, :cell], all_results_list_sr, diff))
    end
end

# Initialize the DataFrame with proper column types
species_count_df = DataFrame(
    species_name = String[],
    count = Int[],
    cell_id = Vector{Int}[],  # Assuming cell IDs are integers
    ind_ext_names = Vector{String}[]
)

# Loop through all results
for dd in all_results_list
    # Find the entry with the lowest survival_rate
    min_survival = minimum(dd.survival_rate)
    min_survived = minimum(dd.survived_herbs + dd.survived_preds)

    # Skip if the minimum survival condition matches the first row's survivors minus 1
    if min_survived == dd[1, :survived_herbs] + dd[1, :survived_preds] - 1
        continue
    end

    # Get all entries with the minimum survival_rate
    min_survival_entries = findall(x -> x == min_survival, dd.survival_rate)
    
    for min_survival_entry in min_survival_entries
        # Get the species name, individual extinction names, and cell ID
        species_name = dd[min_survival_entry, :sp_removed]
        ind_ext_names = dd[min_survival_entry, :ind_ext_name]
        cell_id = dd[min_survival_entry, :cell]

        # Check if species_name is already in the DataFrame
        if species_name in species_count_df.species_name
            # Find the index of the existing entry
            indexx = findfirst(x -> x == species_name, species_count_df.species_name)

            # Update the ind_ext_names and cell_id if they are not already present
            for ind_ext_name in ind_ext_names
                if !(ind_ext_name in species_count_df[indexx, :ind_ext_names])
                    push!(species_count_df[indexx, :ind_ext_names], ind_ext_name)
                end
            end

            if !(cell_id in species_count_df[indexx, :cell_id])
                push!(species_count_df[indexx, :cell_id], cell_id)
            end

            # Increment the count
            species_count_df[indexx, :count] += 1
        else
            # Add a new entry for this species
            push!(species_count_df, (
                species_name,                  # species_name
                1,                             # count starts at 1
                [cell_id],                     # cell_id as a single-element array
                ind_ext_names                  # ind_ext_names as-is
            ))
        end
    end
end

for i in birmmals_names
    if !(i in species_count_df.species_name)
        push!(species_count_df, (
            i,
            0,
            [],
            []
        ))
    end
end
species_count_df[!, :ind_ext_names_length] = length.(species_count_df.ind_ext_names)

begin
    species_count_df = sort(species_count_df, :count, rev=true)
    fig = Figure(resolution = (1000, 800))
    ax = Axis(fig[1, 1], title="Frequency of most influencial species", xlabel="Species", ylabel="Frequency")
    MK.barplot!(ax, 1:length(species_count_df.species_name), species_count_df.count)
    ax.xticks = (1:length(species_count_df.species_name), species_count_df.species_name)
    ax.xticklabelrotation = π/2.5
    ax.xticklabelsize = 7

    sorted_birmmals_biomass_fixed = sort(birmmals_biomass_fixed, :biomass, rev=true)
    ax2 = Axis(fig[1, 2], title="H0_values Biomass", xlabel="Species", ylabel="Biomass")
    MK.barplot!(ax2, 1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.biomass)
    ax2.xticks = (1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.species)
    ax2.xticklabelrotation = π/2.5
    ax2.xticklabelsize = 6
    
    display(fig)
end

begin
    dis = true
    save_the_plot = true
    log = true
    fig = Figure(resolution = (1000, 800))
    sorted_birmmals_biomass_fixed = sort(birmmals_biomass_fixed, :biomass, rev=true)
    ax2 = Axis(
        fig[1, 1], title="H0_values Biomass", 
        xlabel="Species", ylabel="Biomass",
        yscale=log ? log10 : identity
        )
    MK.barplot!(ax2, 1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.biomass)
    ax2.xticks = (1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.species)
    ax2.xticklabelrotation = π/2.5
    ax2.xticklabelsize = 6
    
    if dis
        display(fig)
    end
    word = log ? "log" : "absolute"
    if save_the_plot
        save("Plots/H0_values_biomass_$word.png", fig) 
    end
end

begin
    dis = true
    log = false
    save_the_plot = true
    fig = Figure(resolution = (1000, 800))
    ax = Axis(
        fig[1, 1],
         title="Correlation between H0_values and frequency of most influencial species", 
         xlabel="Frequency", ylabel="H0_values",
         yscale=log ? log10 : identity,
        #  xscale=log ? log10 : identity
    )
    matched_df = innerjoin(species_count_df, sorted_birmmals_biomass_fixed, on=:species_name => :species)
    MK.scatter!(ax, matched_df.count, matched_df.biomass)
    ax.xlabel = "Frequency of most influencial species"
    ax.ylabel = "Biomass of species"
    
    if dis
        display(fig)
    end
    word = log ? "logscale" : ""
    if save_the_plot
        save("Plots/H0_values_VS_most_influencial_$word.png", fig) 
    end
end
