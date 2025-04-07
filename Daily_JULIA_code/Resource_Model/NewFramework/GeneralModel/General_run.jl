function general_run(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 500.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    include_omnivores = true,
    plot = false,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    log = false,
    initial_deviation = 0.0,
    sp_affected = 0,
    extinguishing_species = 0,
    include_n_omnivores = 0,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 0.01,
    force_nu_to = nothing,
    use_cb = false,
    perturbation_halfway = false,
    species_to_perturb = 0,
    removal_fraction = 0.0,
    randomise_attack_rates = false,
    attack_rate_sd = 0.2
)
    AAAA = DataFrame()
    local_i, local_j = idx[cell][1], idx[cell][2]
    # Extract species present in the cell.
    all_sp = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Start with the full species list.
    sp_nm = copy(all_sp)
    
    # Remove predators if not included.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if !(sp in carnivore_names)]
    end

    # Process omnivores:
    if include_omnivores
        # If include_omnivores is true, then all omnivores remain in sp_nm.
        if include_n_omnivores > 0
            error("Cannot specify include_n_omnivores > 0 when include_omnivores is true")
        end
    else
        # Exclude all omnivores.
        sp_nm = [sp for sp in sp_nm if !(sp in omnivore_names)]
        # But if a positive number of omnivores is requested, add them from the cell.
        if include_n_omnivores > 0
            omnivores_in_cell = [sp for sp in all_sp if sp in omnivore_names]
            if length(omnivores_in_cell) < include_n_omnivores
                error("Not enough omnivores in the cell to include the requested number: requested $(include_n_omnivores), found $(length(omnivores_in_cell))")
            end
            sp_nm = union(sp_nm, omnivores_in_cell[1:include_n_omnivores])
        end
    end

    # Also remove species if sp_removed_name is given.
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Get community parameters using the new setup function.
    params_setup = g_attempt_setup_community(
        local_i, local_j,
        mu_val,
        eps_val,
        mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion,
        randomise_attack_rates = randomise_attack_rates,
        attack_rate_sd = attack_rate_sd
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure the returned parameters.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq        # Baseline herbivore (and omnivore) abundances
    P_eq       = params_setup.P_eq        # Baseline predator abundances
    r_i        = params_setup.r_i          # Herbivore intrinsic growth rates
    K_i        = params_setup.K_i          # Herbivore carrying capacities
    mu_val     = params_setup.mu           # Herbivore competition
    nu_val     = params_setup.nu           # Effective predation rate
    A_matrix   = params_setup.A_matrix     # R x S omnivore-herbivore interaction matrix
    C_matrix   = params_setup.C_matrix     # R x S omnivore-omnivore interaction matrix
    epsilon    = params_setup.epsilon      # Predator conversion efficiencies
    m_alpha    = params_setup.m_alpha      # Predator mortality rates
    K_alpha    = params_setup.K_alpha      # Predator carrying capacities
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    # println("Herbivore list: ", length(herbivore_list))

    # Initial conditions: use provided or baseline abundances.
    if isnothing(H_init)
        H_init = H_star
    end
    if true
        if isnothing(P_init)
            P_init = P_star
        end
        if iszero(extinguishing_species)
            initial_deviation_vec = zeros(S+R)
            if sp_affected > 0
                initial_deviation_vec[sp_affected] = initial_deviation
            end
            u0 = abs.(vcat(H_init, P_init) .+ initial_deviation_vec)
        else
            if extinguishing_species > S+R
                @error "Error: extinguishing_species > S+R"                
            end
            u0 = vcat(H_init, P_init) .* [i == extinguishing_species ? 0.0 : 1.0 for i in 1:(S+R)]
        end
    else
        u0 = H_init
    end
    
    # Build the parameter tuple for the ODE.
    # Order: (S, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)
    if !isnothing(force_nu_to)
        nu_val = force_nu_to
    end
    params = (S, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)

    # Define and solve the ODE using omnivore_dynamics!
    prob = ODEProblem(general_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if use_cb
            solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
        else
            solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    end

    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: solution did not finish properly"
            return nothing
        end
    end

    # --- PERTURBATION HALFWAY ---
    if perturbation_halfway && species_to_perturb > 0 && species_to_perturb <= S+R
        u0[species_to_perturb] -= u0[species_to_perturb] * removal_fraction
        new_u0 = deepcopy(sol[:, end])
        new_u0[species_to_perturb] -= new_u0[species_to_perturb] * removal_fraction
        prob = ODEProblem(general_dynamics!, new_u0, (time_end, time_end+time_end), params)
        if use_cb
            new_sol = solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
            else
            new_sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    elseif perturbation_halfway && species_to_perturb > S+R
        error("Error: species_to_perturb > S+R")
    elseif perturbation_halfway && iszero(species_to_perturb)
        error("Error: You activated perturbation_halfway but didn't specify species_to_perturb")
    end

    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (1000, 600))
        ax = Axis(
            fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for cell $cell with mu = $mu_val, epsilon = $eps_val, m_alpha = $mean_m_alpha",
            yscale = log ? log10 : identity
        )
        times_combined = sol.t
        # Plot herbivore/omnivore dynamics: species that are in herbivore_names are blue; those in omnivore_names are green.
        for i in 1:S
            if herbivore_list[i] in omnivore_names
                lines!(ax, times_combined, sol[i, :], label = "O$(i)", color = :green)
                # println("O$(i)")
            else
                lines!(ax, times_combined, sol[i, :], label = "H$(i)", color = :blue)
                # println("H$(i)")
            end
        end
        # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
        for alpha in 1:R
            lines!(ax, times_combined, sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
        end
        hlines!(ax, [0.0], label = "", color = :black, linestyle = :dot)
        

        if perturbation_halfway
            ax2 = Axis(
                fig[1, 2],
                xlabel = "Time",
                ylabel = "Biomass",
                title = "Dynamics for cell after perturbation",
            )
            times_combined = new_sol.t
            # Plot herbivore/omnivore dynamics: species that are in herbivore_names are blue; those in omnivore_names are green.
            for i in 1:S
                if herbivore_list[i] in omnivore_names
                    lines!(ax2, times_combined, new_sol[i, :], label = "O$(i)", color = :green)
                    # println("O$(i)")
                else
                    lines!(ax2, times_combined, new_sol[i, :], label = "H$(i)", color = :blue)
                    # println("H$(i)")
                end
            end
            
            # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
            for alpha in 1:R
                lines!(ax2, times_combined, new_sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
            end
            
            ax3 = Axis(
                fig[2, 1],
                xlabel = "Time",
                ylabel = "Total Biomass",
                title = "Dynamics for cell before and after perturbation",
            )
            times_combined = vcat(sol.t, new_sol.t)
            # println("times combined: ", length(times_combined))
            # println("combined sol: ", typeof(hcat(sol[:, :], new_sol[:, :])))
            # println("combined sol: ", size(hcat(sol[:, :], new_sol[:, :])))
            # println("sum is ", sum(eachrow(hcat(sol[:, :], new_sol[:, :]))))
            lines!(ax3, times_combined, sum(eachrow(hcat(sol[:, :], new_sol[:, :]))))
        end
        display(fig)
    end
    # println("S: $S")
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
    herbivore_survival_rate = survived_herb / S
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    # println("Estimated nu = $nu_val")
    
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
        h_survival              = "$(survived_herb)/$S",
        p_survival              = "$(survived_pred)/$R",
        biomass_at_the_end      = biomass_at_the_end,
        H_vector                = [H_end],
        P_vector                = [P_end],
        H_eq                    = [H_eq],
        P_eq                    = [P_eq],
        survived_herbivores     = survived_herb,
        survived_predators      = survived_pred,
        H_biomass               = H_biomass,
        P_biomass               = P_biomass,
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

    if perturbation_halfway
    
        new_H_end = new_sol[1:S, end]
        new_P_end = new_sol[S+1:S+R, end]
        new_H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_H_end)
        new_P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_P_end)
        new_survived_herb = count(x -> x > EXTINCTION_THRESHOLD, new_H_end)
        new_survived_pred = count(x -> x > EXTINCTION_THRESHOLD, new_P_end)
        new_total_surv = new_survived_herb + new_survived_pred
        new_total_species = S + R
        new_survival_rate = new_total_surv / new_total_species
        new_herbivore_survival_rate = new_survived_herb / S 
        new_predator_survival_rate = (R > 0 ? new_survived_pred / R : 0.0)
        new_H_biomass = sum(new_H_end)
        new_P_biomass = sum(new_P_end)
        new_biomass_at_the_end = new_H_biomass + new_P_biomass
        new_ratio = (new_H_biomass == 0.0 ? NaN : (new_P_biomass / new_H_biomass))

        new_single_run_results = DataFrame(
            cell_id                 = cell,
            i                       = local_i,
            j                       = local_j,
            survival_rate           = new_survival_rate,
            h_survival              = "$(new_survived_herb)/$S",
            p_survival              = "$(new_survived_pred)/$R",
            biomass_at_the_end      = new_biomass_at_the_end,
            H_vector                = [new_H_end],
            P_vector                = [new_P_end],
            H_eq                    = [H_eq],
            P_eq                    = [P_eq],
            survived_herbivores     = new_survived_herb,    
            survived_predators      = new_survived_pred,
            H_biomass               = new_H_biomass,
            P_biomass               = new_P_biomass,
            herb_pred_ratio         = new_ratio,
            total_survivors         = new_total_surv,
            total_species           = new_total_species,
            herbivore_survival_rate = new_herbivore_survival_rate,
            predator_survival_rate  = new_predator_survival_rate,
            mu                      = mu_val,
            nu                      = nu_val,
            epsilon_val             = eps_val,
            mean_m_alpha            = mean_m_alpha
        )
        
        append!(AAAA, new_single_run_results)
    end
    
    initial_biomass = sum(u0)
    final_biomass = sum(sol[:, end])
    println("Initial biomass was $initial_biomass and final biomass is $final_biomass")

    if perturbation_halfway
        if round(new_biomass_at_the_end, digits=0) == round(biomass_at_the_end, digits=0)
            println("The total biomass after the perturbation did not change.")
            if all(round.(vcat(new_H_end, new_P_end), digits=0) .== round.(vcat(H_end, P_end), digits=0))
                println("The equilibrium biomass after the perturbation did not change.")
            end
        else
            println("The total biomass after the perturbation changed and obviously so did the equilibrium biomasses.")
            println("Original biomass: $biomass_at_the_end, new biomass: $new_biomass_at_the_end")
        end
        if new_herbivore_survival_rate != herbivore_survival_rate || new_predator_survival_rate != predator_survival_rate
            println("The survival rates after the perturbation changed.")
        end
        diff = []
        for i in 2:size(new_sol, 2)
            # println("The biomass respect to equilibrium is ", abs(sum(new_sol[:, i]) - biomass_at_the_end))
            push!(diff, abs(sum(new_sol[:, i]) - biomass_at_the_end))
        end
        println("The maximum absolute difference is ", maximum(diff), " and it was found at time ", new_sol.t[argmax(diff)])    
    end
    if do_you_want_params && do_you_want_sol
        return AAAA, params, sol
    elseif do_you_want_params || do_you_want_sol
        return AAAA, do_you_want_params ? params : sol
    else
        return AAAA
    end

end

A_run = general_run(
    1, # cell
    # 0.2, 0.5, 0.1; # STABLE;
    # 0.1, 0.45, 0.3; # semi-STABLE;
    0.2, 0.05, 0.1; # UNSTABLE;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 1000.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = true,
    include_omnivores = true,
    plot = true,
    sp_removed_name = nothing,
    artificial_pi = true, pi_size = 10.0,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = true,
    log = false,
    initial_deviation = 0.0, # leave at 0.0 for no effect
    sp_affected = 38, # leave at 0 for no effect
    extinguishing_species = 0, # leave at 0 for no effect
    include_n_omnivores = 0, # leave at 0 for no effect
    nu_omni_proportion = 1.0, # leave at 1.0 for no effect
    nu_b_proportion = 1.0, # leave at 1.0 for no effect
    r_omni_proportion = 1.0, # leave at 1.0 for no effect
    force_nu_to = nothing, # leave at nothing for no effect
    use_cb = false,
    perturbation_halfway = true,
    species_to_perturb = 38, # leave at 0 for no effect
    removal_fraction = 1.0,
    randomise_attack_rates = true, # If you want you can keep it at true and just play with the sd
    attack_rate_sd = 0.0
)

# cb_no_trigger, cb_trigger = build_callbacks(37, 8, EXTINCTION_THRESHOLD, T_ext, 1)
# for mu in 0.0:0.1:0.9
#     for epsilon in 0.0:0.5:1.0
#         for m_alpha in 0.0:0.1:0.2
#             A_run = general_run(
#                 1, # cell
#                 mu, epsilon, m_alpha; # mu, epsilon, m_alpha
#                 delta_nu = 0.05,
#                 d_alpha = 1.0, d_i = 1.0,
#                 time_end = 1000.0,
#                 do_you_want_params = false,
#                 do_you_want_sol = false,
#                 include_predators = false,
#                 include_omnivores = true,
#                 plot = true,
#                 sp_removed_name = nothing,
#                 artificial_pi = true, pi_size = 10.0,
#                 H_init = nothing,
#                 P_init = nothing,
#                 ignore_inf_error = true,
#                 log = false,
#                 initial_deviation = 0.0, # leave at 0.0 for no effect
#                 extinguishing_species = 0, # leave at 0 for no effect
#                 include_n_omnivores = 0, # leave at 0 for no effect
#                 nu_omni_proportion = 1.0, # leave at 1.0 for no effect
#                 nu_b_proportion = 1.0, # leave at 1.0 for no effect
#                 r_omni_proportion = 0.01, # leave at 1.0 for no effect
#                 force_nu_to = nothing, # leave at nothing for no effect
#                 use_cb = false,
#                 perturbation_halfway = false,
#                 species_to_perturb = 0, # leave at 0 for no effect
#                 removal_fraction = 0.0
#             )
#         end
#     end
# end


