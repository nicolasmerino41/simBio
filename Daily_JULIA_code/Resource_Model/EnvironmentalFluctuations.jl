# CHOOSE A CELL
# Extract best param_combinations from best_params_per_cell
best_params = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/best_params_5950_cells_not_rounded.csv", DataFrame)
best_params = filter(row -> row.g_iH_i_over_NPP < 10.0, best_params)

i = 1

# SIMULATE ENVIRONMENTAL FLUCTUATIONS
begin

    target_sv = best_params[i, :survival_rate]
    cell = best_params[i, :cell_id]
        
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
        localNPP = Float64(npp_DA_relative_to_1000[local_i, local_j])
        localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

        # cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)

        ###############################
        # 4) Run Full Community Simulation (Baseline: no species removed)
        ###############################
        # Attempt setup
        full_results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_predation_val, epsilon_val, sym_competition;
            localNPP      = localNPP,
            localH0_vector= localH0_vector
        )
        if full_results === nothing
            @error  "full results are nothing"  
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
    params_full = (S_full, R_full, H_i0_full, m_i_full, g_i_full, G_full,
                   M_modified_full, a_matrix_full, A_full, epsilon_vector_full, m_alpha_full)
    tspan = (0.0, 500.0)
    prob_full = ODEProblem(ecosystem_dynamics!, u0_full, tspan, params_full)
    sol_full = solve(prob_full, Tsit5(); callback = cb_no_trigger, reltol=1e-6, abstol=1e-6)

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

    times = sol_full.t
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass", title="Simulation")
    for i in 1:S_full
        lines!(ax, times, sol_full[i, :], label="H$i")
    end
    for α in 1:R_full
        lines!(ax, times, sol_full[S_full + α, :], label="P$α", linestyle=:dash)
    end
    # axislegend(ax)
    display(fig)
end

#### APOCALYPSE ####
begin
    apocalypse_vector = rand(0.5:0.1:1.0, S_full)
    # Build initial conditions
    H_init_full = sol_full[1:S_full, end]
    if R_full > length(H_init_full)
        error("In cell $cell: More predators than herbivore data available.")
    end
    P_init_full = sol_full[1:R_full, end] ./ 10.0
    u0_full = vcat(H_init_full, P_init_full)
    params_full = (S_full, R_full, H_i0_full.*apocalypse_vector, m_i_full, g_i_full, G_full,
               M_modified_full, a_matrix_full, A_full, epsilon_vector_full, m_alpha_full)
    tspan = (0.0, 500.0)
    prob_full = ODEProblem(ecosystem_dynamics!, u0_full, tspan, params_full)
    sol_full_after = solve(prob_full, Tsit5(); callback = cb_no_trigger, reltol=1e-6, abstol=1e-6)

    if sol_full_after.t[end] < 500.0 || any(isnan, sol_full_after.u[end]) || any(isinf, sol_full_after.u[end])
        full_survival_rate = 0.0
        full_H_biomass = 0.0
        full_P_biomass = 0.0
        full_surv_herb = 0
        full_surv_pred = 0
    else
        H_end_full = sol_full_after[1:S_full, end]
        P_end_full = sol_full_after[S_full+1:S_full+R_full, end]
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

    # Combine time vectors and solutions for continuous plotting
    times_combined = vcat(sol_full.t, sol_full_after.t[2:end].+500.0)  # Avoid duplicating the first time point
    H_combined = hcat(sol_full[1:S_full, :], sol_full_after[1:S_full, 2:end])  # Merge herbivores' data
    P_combined = hcat(sol_full[S_full+1:S_full+R_full, :], sol_full_after[S_full+1:S_full+R_full, 2:end])  # Merge predators' data

    # Plot both simulations on the same graph
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass", title="Combined Simulation (Baseline and Apocalypse)")

    # Plot herbivores
    for i in 1:S_full
        lines!(ax, times_combined, H_combined[i, :], label="H$i")
    end

    # Plot predators
    for α in 1:R_full
        lines!(ax, times_combined, P_combined[α, :], label="P$α", linestyle=:dash)
    end

    # Add a legend to the plot
    # axislegend(ax, position=:rt)
    display(fig)    
end

index = findfirst(x -> x == "Oenanthe oenanthe", names(iberian_interact_NA, 1))

println(iberian_interact_NA["Oenanthe oenanthe", :])
cucu = DataFrame(species = names(iberian_interact_NA, 1), value = iberian_interact_NA[:, "Oenanthe oenanthe", ])