# INSIGHTS: RUN THE SCRIPTS IN THE FOLLOWING ORDER:
# 0) o_anal_num_equilibrium
# 1) PlottingEigenvalues
# 2) Elasticity
# 3) MultipleCells
# 4) SensitivityVsEcosystemFunctioning
execute_code = false
function analytical_equilibrium(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    include_omnivores = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 1.0,
    callbacks = false,
    plot = false
)

    # Convert parameters to Float64:
    mu_val         = to_float(mu_val)
    eps_val        = to_float(eps_val)
    mean_m_alpha   = to_float(mean_m_alpha)
    delta_nu       = to_float(delta_nu)
    d_alpha        = to_float(d_alpha)
    d_i            = to_float(d_i)
    pi_size        = to_float(pi_size)
    nu_omni_proportion = to_float(nu_omni_proportion)
    nu_b_proportion    = to_float(nu_b_proportion)
    r_omni_proportion  = to_float(r_omni_proportion)

    local_i, local_j = idx[cell][1], idx[cell][2]
    # Extract species names from the cell.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species: if predators are not included, restrict to herbivores & omnivores.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names || sp in omnivore_names]
    end
    # Process omnivores:
    if !include_omnivores
        # Exclude all omnivores.
        sp_nm = [sp for sp in sp_nm if !(sp in omnivore_names)]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Use the new community setup function.
    params_setup = o_attempt_setup_community(
        local_i, local_j,
        mu_val, eps_val, mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure parameters from the new setup.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq
    P_eq       = params_setup.P_eq
    r_i        = params_setup.r_i
    K_i        = params_setup.K_i
    nu_val     = params_setup.nu
    P_matrix   = params_setup.P_matrix
    O_matrix   = params_setup.O_matrix
    T_matrix   = params_setup.T_matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    B_matrix   = params_setup.B_matrix
    D_matrix   = params_setup.D_matrix
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    predator_list  = params_setup.predator_list
    
    # Set initial conditions.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    # println("typeof(u0): ", typeof(u0))
    
    # Build the parameter tuple for omnivore_dynamics!
    # Order: (S, R, K_i, r_i, mu, nu, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_omni)
    p = (
        S, R, K_i, r_i,
        mu_val, nu_val,
        P_matrix, O_matrix, T_matrix,
        epsilon, m_alpha,
        K_alpha, B_matrix, D_matrix,
        nu_val*nu_omni_proportion, nu_val*nu_b_proportion
    )
    
    # Define an inner wrapper so we can compute the Jacobian at u0.
    function f_wrapper(u)
        du = similar(u)
        omnivore_dynamics!(du, u, p, 0.0)
        return du
    end
    
    J = ForwardDiff.jacobian(f_wrapper, u0)
    # println("Jacobian at equilibrium:")
    # println(J)
    
    # Define and solve the ODE using omnivore_dynamics!
    prob = ODEProblem(omnivore_dynamics!, u0, (0.0, 1000.0), p)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if callbacks
            solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
        else
            solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    end

    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (600, 500))
        ax = Axis(
            fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for cell"
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
        display(fig)
    end

    # Package results into a NamedTuple.
    return (
        equilibrium = (
            H_star = H_star, P_star = P_star,
            u0 = u0,
            H_vector = sol[1:S, end],
            P_vector = sol[(S+1):(S+R), end],
            herbivore_list = herbivore_list,
            predator_list = predator_list
            ),
        parameters  = (
            S = S,
            R = R,
            K_i = K_i,
            r_i = r_i,
            mu = mu_val,
            nu = nu_val,
            P_matrix = P_matrix,
            O_matrix = O_matrix,
            T_matrix = T_matrix,
            epsilon = epsilon,
            m_alpha = m_alpha,
            K_alpha = K_alpha,
            B_matrix = B_matrix,
            D_matrix = D_matrix,
            nu_omni = nu_val*nu_omni_proportion,
            nu_b = nu_val*nu_b_proportion
        ),
        Jacobian = J
    )
end

# --- Compute equilibrium and assess local stability ---
if execute_code
    Threads.@threads for j in 0.0:0.01:1.0
        for i in 0.0:0.1:1.0
            for k in 0.0:0.1:1.0
                result = analytical_equilibrium( # 0.29 for cell 1 and 0.32 for cell 2
                    cell, 
                    1.0, 0.12, 1.0;
                    # i, j, k;
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = true, pi_size = 10.0,
                    H_init = nothing,
                    P_init = nothing,
                    nu_omni_proportion = 1.0,
                    nu_b_proportion = 1.0,
                    r_omni_proportion = 1.0,
                    callbacks = false,
                    plot = true
                );

                A_eq = result.equilibrium       # Contains H_star, P_star, and u0 (equilibrium state)
                A_p  = result.parameters        # Contains the ODE parameters
                A_j  = result.Jacobian          # The Jacobian at equilibrium

                # println("Jacobian at equilibrium:")
                # println(A_j)

                # Compute eigenvalues to assess local stability.
                eigvals = eigen(A_j).values
                # println("Eigenvalues:")
                # println(eigvals)
                stable = all(real.(eigvals) .< 0)
                if stable
                    println("The equilibrium is locally stable for mu = $(i), eps = $(j), m_alpha = $(k)")
                elseif !stable && A_eq.H_star == A_eq.H_vector && A_eq.P_star == A_eq.P_vector
                    println("The equilibrium is locally unstable")# but works for mu = $(i), eps = $(j), m_alpha = $(k)")
                end
            end
        end
    end
end
##################################################################################
##################################################################################
################### INDIVIDUAL PERTURBATION RESPONSE #############################
##################################################################################
##################################################################################
function plot_perturbation_response(
    A_eq, A_p;
    perturbation=0.01, tspan=(0.0,50.0),
    callbacks=true
    )
    S = length(A_eq.H_star)
    R = length(A_eq.P_star)
    # A_eq.u0 is the equilibrium state vector.
    u0 = A_eq.u0
    n = length(u0)
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "Time", ylabel = "Difference from equilibrium", 
              title = "Response to $(perturbation*100)% Perturbations")
    
    sp_nm = vcat(A_eq.herbivore_list, A_eq.predator_list)

    # Loop over each species index.
    for i in 1:S+R
        # Create a copy of the equilibrium state and perturb the i-th component.
        u0_perturbed = copy(u0)
        u0_perturbed[i] += perturbation * u0[i]
        
        # Define and solve the ODE with the perturbed initial condition.
        prob = ODEProblem(omnivore_dynamics!, u0_perturbed, tspan, A_p)
        if !callbacks
            sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        else
            sol = solve(prob, Tsit5(); callback = cb_no_trigger, abstol=1e-8, reltol=1e-6)
        end

        t = sol.t
        # Compute the difference over time for the i-th species.
        diff = [sol(u)[i] - u0[i] for u in t]
        
        # Plot the response for species i.
        if sp_nm[i] in herbivore_names
            lines!(ax, t, diff, label = "Herbivore $i", color = :blue)
        elseif sp_nm[i] in omnivore_names
            lines!(ax, t, diff, label = "Omnivore $i", color = :green)
        elseif sp_nm[i] in predator_names
            lines!(ax, t, diff, label = "Predator $i", color = :red)
        end
    end

    display(fig)
    # axislegend(ax)
    fig
end

# Example usage:
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function:
if execute_code
    fig = plot_perturbation_response(A_eq, A_p; perturbation=0.01, callbacks=false)
end
#######################################################################################
#######################################################################################
############################ TOTAL BIOMASS RESPONSE ###################################
#######################################################################################
#######################################################################################
function plot_total_biomass_response(
    A_eq, A_p; 
    perturbation=0.01, tspan=(0.0, 50.0), separate_guilds=false,
    callbacks=true
)
    # Extract the equilibrium state and compute the total biomass.
    u0 = A_eq.u0
    equilibrium_total = sum(u0)
    n = length(u0)
    println("Equilibrium total biomass: ", equilibrium_total)
    fig = Figure(resolution = (800, 600))
    if !separate_guilds 
        ax = Axis(fig[1,1], xlabel = "Time", ylabel = "Difference in Total Biomass", 
              title = "Total Biomass Response to $(perturbation*100)% Perturbations")
    else
        ax = Axis(fig[1,1], xlabel = "Time", ylabel = "Difference in Total Biomass",
                   title = "Total Biomass Response to $(perturbation*100)% Perturbations")
        ax1 = Axis(fig[1,2], xlabel = "Time", ylabel = "Difference in Total Biomass",
                   title = "Herbivore Biomass Response to $(perturbation*100)% Perturbations")
        ax2 = Axis(fig[2,1], xlabel = "Time", ylabel = "Difference in Total Biomass",
                   title = "Omnivore Biomass Response to $(perturbation*100)% Perturbations")
        ax3 = Axis(fig[2,2], xlabel = "Time", ylabel = "Difference in Total Biomass",
                   title = "Predator Biomass Response to $(perturbation*100)% Perturbations")
    end

    sp_nm = vcat(A_eq.herbivore_list, A_eq.predator_list)
    # Loop over each species index.
    for i in 1:n
        # Perturb the i-th species.
        u0_perturbed = copy(u0)
        u0_perturbed[i] += perturbation * u0[i]
        p_tuple = (A_p.S, A_p.R, A_p.K_i, A_p.r_i, A_p.mu, A_p.nu, A_p.P_matrix, A_p.O_matrix, A_p.T_matrix, A_p.epsilon, A_p.m_alpha, A_p.K_alpha, A_p.B_matrix, A_p.D_matrix, A_p.nu_omni, A_p.nu_b)
        # Solve the ODE with the perturbed initial condition.
        prob = ODEProblem(omnivore_dynamics!, u0_perturbed, tspan, p_tuple)
        if !callbacks
            sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        else
            sol = solve(prob, Tsit5(); callback = cb_no_trigger, abstol=1e-8, reltol=1e-6)
        end
        
        # Compute the total biomass difference over time.
        t = sol.t
        diff = [equilibrium_total + (sum(sol(u)) - equilibrium_total) for u in t]
        println("Total biomass difference for species $i: ", diff[end])
        if !separate_guilds
            if sp_nm[i] in herbivore_names
                lines!(ax, t, diff, label = "Herbivore $i", color = :blue)
            elseif sp_nm[i] in omnivore_names
                lines!(ax, t, diff, label = "Omnivore $i", color = :green)
            elseif sp_nm[i] in predator_names
                lines!(ax, t, diff, label = "Predator $i", color = :red)
            end
        else
            if sp_nm[i] in herbivore_names
                lines!(ax, t, diff, label = "Herbivore $i", color = :blue)
            elseif sp_nm[i] in omnivore_names
                lines!(ax, t, diff, label = "Omnivore $i", color = :green)
            elseif sp_nm[i] in predator_names
                lines!(ax, t, diff, label = "Predator $i", color = :red)
            end
            if sp_nm[i] in herbivore_names
                lines!(ax1, t, diff, label = "Herbivore $i", color = :blue)
            elseif sp_nm[i] in omnivore_names
                lines!(ax2, t, diff, label = "Omnivore $i", color = :green)
            elseif sp_nm[i] in predator_names
                lines!(ax3, t, diff, label = "Predator $i", color = :red)
            end
        end
    end
    display(fig)
    # axislegend(ax)
    return fig
end

# Example usage:
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function:
if execute_code
    fig = plot_total_biomass_response(
    A_eq, A_p;
    perturbation=1.0, tspan=(0.0, 1000.0), separate_guilds=false,
    callbacks=false
    )
end

#######################################################################################
#######################################################################################
############################ SENSITIVITY VS METRICS ###################################
#######################################################################################
#######################################################################################
# --- Function to compute sensitivity metrics and species traits ---
function compute_sensitivity_metrics(
    A_eq, A_p;
    perturbation=0.01, tspan=(0.0, 50.0),
    callbacks=true,
    tolerance_factor=1e-3
)
    # Extract the equilibrium state vector and total equilibrium biomass.
    u0 = A_eq.u0
    total_eq = sum(u0)
    n = length(u0)
    sensitivity = zeros(n)
    resilience = fill(NaN, n)  # resilience: recovery time for each species
    
    # For each species, perturb its equilibrium value by a small fraction,
    # simulate the system using omnivore_dynamics! (the full model),
    # and compute the maximum absolute deviation in total biomass relative to equilibrium.
    for i in 1:n
         u0_perturbed = copy(u0)
         u0_perturbed[i] -= perturbation * u0[i] 
         prob = ODEProblem(omnivore_dynamics!, u0_perturbed, tspan, A_p)
         if !callbacks
             sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
         else
            sol = solve(prob, Tsit5(); callback = cb_no_trigger, abstol=1e-8, reltol=1e-6)
         end
         total_over_time = vec(sum(sol[:, :], dims=1))  # gives a vector of total biomass at each time
    deviations = abs.(total_over_time .- total_eq)
    sensitivity[i] = maximum(deviations)

    # Resilience: first time at which total biomass is within tolerance of total_eq.
    tolerance = tolerance_factor * total_eq
    rec_time = NaN
    for (j, t_val) in enumerate(sol.t)
        if abs(total_over_time[j] - total_eq) < tolerance
            rec_time = t_val
            break
        end
    end
    resilience[i] = rec_time

    end
    
    # Compute species traits.
    # Trait 1: Equilibrium biomass for each species.
    biomass = copy(u0)
    
    # Trait 2: "Degree" (number of interactions).
    # Here we combine several matrices:
    # - For herbivores (indices 1:S): we add:
    #     • P_matrix[i, :] : number of predators feeding on species i,
    #     • O_matrix[i, :] : number of omnivory interactions where species i acts as consumer,
    #     • T_matrix[i, :] : number of herbivore-herbivore interactions (losses) acting on species i.
    # - For predators (indices S+1:n): we add:
    #     • P_matrix[:, α] : number of herbivores preyed upon by predator α,
    #     • B_matrix[α, :] : number of beneficial predator–predator interactions (α consuming others),
    #     • D_matrix[α, :] : number of predator–predator interactions where α is prey.
    S = A_p.S
    R = A_p.R
    P_matrix = A_p.P_matrix
    O_matrix = A_p.O_matrix
    T_matrix = A_p.T_matrix
    B_matrix = A_p.B_matrix
    D_matrix = A_p.D_matrix
    sp_nm = vcat(A_eq.herbivore_list, A_eq.predator_list)
    guilds = [sp_nm[i] in herbivore_names ? "Herbivore" :
              sp_nm[i] in omnivore_names ? "Omnivore" :
              sp_nm[i] in predator_names ? "Predator" : "Unknown" for i in 1:n]

    degree = zeros(n)
    for i in 1:S
        # If the matrices are binary, this gives the number of interactions.
        degree[i] = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
    end

    connectance = zeros(n)
    
    # For herbivores (and omnivores; indices 1:S)
    for i in 1:S
        # Realized interactions for species i:
        #   - Predation losses from predators: P_matrix[i, :]
        #   - Omnivory interactions where i is consumer: O_matrix[i, :]
        #   - Herbivore-herbivore interactions (e.g., being consumed): T_matrix[i, :]
        d = sum(P_matrix[i, :]) + sum(O_matrix[i, :]) + sum(T_matrix[i, :])
        # Maximum possible interactions for a herbivore:
        #   - It could, in principle, be attacked by all R predators.
        #   - It could interact with all other (S - 1) herbivores in both directions (consume and be consumed).
        max_possible = R + 2 * (S - 1)
        connectance[i] = d / max_possible
    end
    
    # For predators (indices S+1:S+R)
    for alpha in 1:R
        # Realized interactions for predator alpha:
        #   - It can feed on all S herbivores (from P_matrix[:, alpha]).
        #   - It has predator–predator interactions: both beneficial (B_matrix[alpha, :]) 
        #     and losses (D_matrix[alpha, :]).
        d = sum(P_matrix[:, alpha]) + sum(B_matrix[alpha, :]) + sum(D_matrix[alpha, :])
        # Maximum possible for a predator:
        #   - It could prey on all S herbivores.
        #   - It could interact with all other (R - 1) predators in two ways.
        max_possible = S + 2 * (R - 1)
        connectance[S+alpha] = d / max_possible
    end
     
    return (sensitivity = sensitivity,  resilience = resilience, biomass = biomass, degree = degree, connectance = connectance, guild = guilds)
end

# A_caca = compute_sensitivity_metrics(A_eq, A_p; perturbation = 1.0, tspan = (0.0, 100000.0), callbacks = false).resilience
# --- Function to plot sensitivity vs species traits using Makie ---
function plot_sensitivity_traits(metrics)
    fig = Figure(resolution = (1100, 700))
    ax1 = Axis(fig[1,1], xlabel = "Equilibrium Biomass", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Equilibrium Biomass")
    colors = [metrics.guild[i] == "Herbivore" ? :blue : metrics.guild[i] == "Predator" ? :red : :green for i in 1:length(metrics.guild)]
    MK.scatter!(ax1, metrics.biomass, metrics.sensitivity, markersize = 10, color = colors)
    
    ax2 = Axis(fig[1,2], xlabel = "Degree (Number of Interactions)", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Degree")
    MK.scatter!(ax2, metrics.degree, metrics.sensitivity, markersize = 10, color = colors)
    
    ax3 = Axis(fig[2,1], xlabel = "Connectance", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Connectance")
    MK.scatter!(ax3, metrics.connectance, metrics.sensitivity, markersize = 10, color = colors)
    display(fig)
    return fig
end

# --- Example Usage ---
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function
if execute_code
    metrics = compute_sensitivity_metrics(A_eq, A_p; perturbation = 1.0, tspan = (0.0, 50.0), callbacks = false)
    fig_traits = plot_sensitivity_traits(metrics)
    display(fig_traits)
end