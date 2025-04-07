function g_abstract_run(
    S::Int, O::Int, R::Int,
    mu_val,
    eps_val,
    mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    time_end = 500.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    plot = false,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 0,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 0.01,
    force_nu_to = nothing,
    use_cb = false,
    perturbation_halfway = false,
    species_to_perturb = 0,
    removal_fraction = 0.0,
    conn = 0.2,  # target connectance for synthetic interaction matrices
    cell_abundance_h::Vector{Float64} = Float64[],  # if provided, else defaults will be ones
    cell_abundance_p::Vector{Float64} = Float64[]
)
    AAAA = DataFrame()

    # Build the synthetic community using the new abstract function.
    params_setup = g_abstract_parametrise_the_community(
        S, O, R;
        mu = mu_val,
        epsilon_val = eps_val,
        mean_m_alpha = mean_m_alpha,
        conn = conn,
        cell_abundance_h = cell_abundance_h,
        cell_abundance_p = cell_abundance_p,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end

    # Destructure the returned parameters.
    S_total    = params_setup.S       # Total number of herbivores+omnivores
    R          = params_setup.R
    H_eq       = params_setup.H_eq      # Baseline herbivore (and omnivore) abundances
    P_eq       = params_setup.P_eq      # Baseline predator abundances
    r_i        = params_setup.r_i       # Herbivore intrinsic growth rates
    K_i        = params_setup.K_i       # Herbivore carrying capacities
    nu_val     = params_setup.nu        # Effective predation rate
    A_matrix   = params_setup.A_matrix   # Trophic interaction matrix
    C_matrix   = params_setup.C_matrix   # Competition interaction matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    predator_list  = params_setup.predator_list

    # println("Herbivore list length: ", length(herbivore_list))
    
    # Initial conditions: if not provided, use the baseline equilibrium.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        if extinguishing_species == 0
            u0 = abs.(vcat(H_init, P_init) .+ randn(S_total+R) .* initial_deviation)
        else
            if extinguishing_species > S_total+R
                @error "Error: extinguishing_species > total species"
            end
            u0 = vcat(H_init, P_init) .* [i == extinguishing_species ? 0.0 : 1.0 for i in 1:(S_total+R)]
        end
    else
        u0 = H_init
    end

    # Build the parameter tuple for the ODE.
    if !isnothing(force_nu_to)
        nu_val = force_nu_to
    end
    params = (S+O, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)

    # Solve the ODE using omnivore_dynamics!
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

    # --- PERTURBATION HALFWAY (optional) ---
    if perturbation_halfway
        if species_to_perturb > 0 && species_to_perturb <= S_total+R
            u0_perturbed = sol[:, end]
            u0_perturbed[species_to_perturb] -= u0_perturbed[species_to_perturb] * removal_fraction
            prob2 = ODEProblem(general_dynamics!, u0_perturbed, (time_end, time_end+time_end), params)
            if use_cb
                new_sol = solve(prob2, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
            else
                new_sol = solve(prob2, Tsit5(); abstol = 1e-8, reltol = 1e-6)
            end
        elseif species_to_perturb > S_total+R
            error("Error: species_to_perturb > total species")
        else
            error("Error: perturbation_halfway activated but species_to_perturb is not specified")
        end
    end

    # --- Plotting using Makie ---
    if plot    
        fig = Figure(; size = (900, 600))
        ax = Axis(fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for abstract community",
            yscale = log ? log10 : identity)
        times = sol.t
        for i in 1:S_total
            # For simplicity, use blue for herbivores and green for omnivores.
            color = (herbivore_list[i] in herbivore_list && !occursin("Omnivore", herbivore_list[i])) ? :blue : :green
            lines!(ax, times, sol[i, :], label = "$(herbivore_list[i])", color = color)
        end
        for alpha in 1:R
            lines!(ax, times, sol[S_total+alpha, :], label = "$(predator_list[alpha])", linestyle = :dash, color = :red)
        end

        if perturbation_halfway
            ax2 = Axis(fig[1, 2],
                xlabel = "Time",
                ylabel = "Biomass",
                title = "After perturbation",
                yscale = log ? log10 : identity)
            times2 = new_sol.t
            for i in 1:S_total
                color = (herbivore_list[i] in herbivore_list && !occursin("Omnivore", herbivore_list[i])) ? :blue : :green
                lines!(ax2, times2, new_sol[i, :], label = "$(herbivore_list[i])", color = color)
            end
            for alpha in 1:R
                lines!(ax2, times2, new_sol[S_total+alpha, :], label = "$(predator_list[alpha])", linestyle = :dash, color = :red)
            end

            # # Plot total biomass before and after perturbation (excluding data between time_end and time_end+20)
            # ax3 = Axis(fig[2, 1],
            #     xlabel = "Time",
            #     ylabel = "Total Biomass",
            #     title = "Total Biomass Dynamics")
            # filtered_old = sol[:, sol.t .< time_end]
            # filtered_new = new_sol[:, new_sol.t .> time_end + 20]
            # times_old = sol.t[sol.t .< time_end]
            # times_new = new_sol.t[new_sol.t .> time_end + 20]
            # combined_times = vcat(times_old, times_new)
            # combined_sol = hcat(filtered_old, filtered_new)
            # total_biomass = vec(sum(combined_sol, dims=1))
            # lines!(ax3, combined_times, total_biomass, label = "Total Biomass", color = :black)
        end
        display(fig)
    end

    # --- Post-Processing: Extract final densities and compute survival rates ---
    H_end = sol[1:S_total, end]
    P_end = sol[S_total+1:S_total+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S_total + R
    survival_rate = total_surv / total_species
    herbivore_survival_rate = survived_herb / S_total
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    println("Estimated nu = $nu_val")

    single_run_results = DataFrame(
        survival_rate           = survival_rate,
        h_survival              = "$(survived_herb)/$S_total",
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
        new_H_end = new_sol[1:S_total, end]
        new_P_end = new_sol[S_total+1:S_total+R, end]
        new_H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_H_end)
        new_P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_P_end)
        new_survived_herb = count(x -> x > EXTINCTION_THRESHOLD, new_H_end)
        new_survived_pred = count(x -> x > EXTINCTION_THRESHOLD, new_P_end)
        new_total_surv = new_survived_herb + new_survived_pred
        new_total_species = S_total + R
        new_survival_rate = new_total_surv / new_total_species
        new_herbivore_survival_rate = new_survived_herb / S_total 
        new_predator_survival_rate = (R > 0 ? new_survived_pred / R : 0.0)
        new_H_biomass = sum(new_H_end)
        new_P_biomass = sum(new_P_end)
        new_biomass_at_the_end = new_H_biomass + new_P_biomass
        new_ratio = (new_H_biomass == 0.0 ? NaN : (new_P_biomass / new_H_biomass))

        new_single_run_results = DataFrame(
            survival_rate           = new_survival_rate,
            h_survival              = "$(new_survived_herb)/$S_total",
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
        diff = []
        for i in 2:size(new_sol, 2)
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

begin
    S, O, R = 35, 5, 8
    # cb_no_trigger, cb_trigger = build_callbacks(S+O, R, EXTINCTION_THRESHOLD, T_ext, 1)
    @time A_run = g_abstract_run(
        S, O, R,    
        0.5, 0.29, 0.1; # mu, epsilon, m_alpha
        delta_nu = 0.05,
        d_alpha = 1.0, d_i = 1.0,
        time_end = 500.0,
        do_you_want_params = false,
        do_you_want_sol = false,
        plot = true,
        H_init = nothing,
        P_init = nothing,
        ignore_inf_error = true,
        log = false,
        initial_deviation = 0.0,
        extinguishing_species = 0,
        nu_omni_proportion = 1.0,
        nu_b_proportion = 1.0,
        r_omni_proportion = 1.0,
        force_nu_to = nothing,
        use_cb = false,
        perturbation_halfway = false,
        species_to_perturb = 10,
        removal_fraction = 0.1,
        conn = 0.1,  # target connectance for synthetic interaction matrices
        cell_abundance_h = [i <= S ? 10.0 : 5.0 for i in 1:S+O],  # if provided, else defaults will be ones
        cell_abundance_p = [1.0 for i in 1:R]
    )
end