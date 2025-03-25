using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots

function analytical_equilibrium(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing
)
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species if needed.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Get community parameters.
    params_setup = b_attempt_setup_community(
        local_i, local_j,
        mu_val,
        eps_val,
        mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure parameters.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq
    P_eq       = params_setup.P_eq
    r_i        = params_setup.r_i
    K_i        = params_setup.K_i
    nu_val     = params_setup.nu
    P_matrix   = params_setup.P_matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    
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
    println("typeof(u0): ", typeof(u0))
    
    # Build parameter tuple for the ODE.
    p = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, epsilon, m_alpha, K_alpha)
    
    # Define f_wrapper as an inner function so it can access p.
    function f_wrapper(u)
        du = similar(u)
        bipartite_dynamics!(du, u, p, 0.0)
        return du
    end
    
    # Compute the Jacobian at equilibrium using ForwardDiff.
    J = ForwardDiff.jacobian(f_wrapper, u0)
    println("Jacobian at equilibrium:")
    println(J)
    
    # Package results into a NamedTuple.
    return (
        equilibrium = (H_star = H_star, P_star = P_star, u0 = u0),
        parameters  = (
            S = S,
            R = R,
            r_i = r_i,
            K_i = K_i,
            mu = mu_val,
            nu = nu_val,
            P_matrix = P_matrix,
            epsilon = epsilon,
            m_alpha = m_alpha,
            K_alpha = K_alpha
        ),
        Jacobian = J
    )
end

# --- Compute equilibrium and assess local stability ---
result = analytical_equilibrium(
    1, 
    0.5, 1.0, 0.1;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing
)

A_eq = result.equilibrium       # Contains H_star, P_star, and u0 (equilibrium state)
A_p  = result.parameters        # Contains the ODE parameters
A_j  = result.Jacobian          # The Jacobian at equilibrium

println("Jacobian at equilibrium:")
println(A_j)

# Compute eigenvalues to assess local stability.
eigvals = eigen(A_j).values
println("Eigenvalues:")
println(eigvals)
stable = all(real.(eigvals) .< 0)
println("Is the equilibrium locally stable? ", stable)

function plot_perturbation_response(A_eq, A_p; perturbation=0.01, tspan=(0.0,50.0))
    S = length(A_eq.H_star)
    R = length(A_eq.P_star)
    # A_eq.u0 is the equilibrium state vector.
    u0 = A_eq.u0
    n = length(u0)
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "Time", ylabel = "Difference from equilibrium", 
              title = "Response to $(perturbation*100)% Perturbations")
    
    # Loop over each species index.
    for i in 1:S
        # Create a copy of the equilibrium state and perturb the i-th component.
        u0_perturbed = copy(u0)
        u0_perturbed[i] += perturbation * u0[i]
        
        # Define and solve the ODE with the perturbed initial condition.
        prob = ODEProblem(bipartite_dynamics!, u0_perturbed, tspan, A_p)
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        
        t = sol.t
        # Compute the difference over time for the i-th species.
        diff = [sol(u)[i] - u0[i] for u in t]
        
        # Plot the response for species i.
        lines!(ax, t, diff, label = "Herbivore $i", color = :blue)
    end
    for i in S+1:S+R
        # Create a copy of the equilibrium state and perturb the i-th component.
        u0_perturbed = copy(u0)
        u0_perturbed[i] += perturbation * u0[i]
        
        # Define and solve the ODE with the perturbed initial condition.
        prob = ODEProblem(bipartite_dynamics!, u0_perturbed, tspan, A_p)
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        
        t = sol.t
        # Compute the difference over time for the i-th species.
        diff = [sol(u)[i] - u0[i] for u in t]
        
        # Plot the response for species i.
        lines!(ax, t, diff, label = "Predator $(i)", color = :red)
    end

    display(fig)
    # axislegend(ax)
    fig
end

# Example usage:
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function:
fig = plot_perturbation_response(A_eq, A_p; perturbation=0.01)

function plot_total_biomass_response(A_eq, A_p; perturbation=0.01, tspan=(0.0, 50.0))
    # Extract the equilibrium state and compute the total biomass.
    u0 = A_eq.u0
    equilibrium_total = sum(u0)
    n = length(u0)
    println("Equilibrium total biomass: ", equilibrium_total)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1,1], xlabel = "Time", ylabel = "Difference in Total Biomass", 
              title = "Total Biomass Response to $(perturbation*100)% Perturbations")
    
    # Loop over each species index.
    for i in 1:n
        # Perturb the i-th species.
        u0_perturbed = copy(u0)
        u0_perturbed[i] += perturbation * u0[i]
        p_tuple = (A_p.S, A_p.R, A_p.K_i, A_p.r_i, A_p.mu, A_p.nu, A_p.P_matrix, A_p.epsilon, A_p.m_alpha, A_p.K_alpha)
        # Solve the ODE with the perturbed initial condition.
        prob = ODEProblem(bipartite_dynamics!, u0_perturbed, tspan, p_tuple)
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        
        # Compute the total biomass difference over time.
        t = sol.t
        diff = [sum(sol(u)) - equilibrium_total for u in t]
        println("Total biomass difference for species $i: ", diff[end])
        # Plot the response for species i.
        lines!(ax, t, diff, label = "Species $i")
    end
    display(fig)
    # axislegend(ax)
    return fig
end

# Example usage:
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function:
fig = plot_total_biomass_response(A_eq, A_p; perturbation=0.01, tspan=(0.0, 50.0))

######### PIPELINE ########
# --- Function to compute sensitivity metrics and species traits ---
function compute_sensitivity_metrics(A_eq, A_p; perturbation=0.01, tspan=(0.0, 50.0))
    # Extract equilibrium state and total equilibrium biomass.
    u0 = A_eq.u0
    total_eq = sum(u0)
    n = length(u0)
    sensitivity = zeros(n)
    
    # For each species, perturb its equilibrium value by a small fraction,
    # simulate the system, and compute the maximum absolute deviation
    # in total biomass relative to the equilibrium total.
    for i in 1:n
         u0_perturbed = copy(u0)
         u0_perturbed[i] += perturbation * u0[i]
         prob = ODEProblem(bipartite_dynamics!, u0_perturbed, tspan, A_p)
         sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
         # Compute the absolute deviation in total biomass at each time point.
         deviations = [abs(sum(sol(u)) - total_eq) for u in sol.t]
         sensitivity[i] = maximum(deviations)
    end
    
    # Compute species traits.
    # Trait 1: Equilibrium biomass for each species.
    biomass = copy(u0)
    
    # Trait 2: "Degree" (number of interactions).
    # For herbivores (indices 1:S): degree is the number of predators feeding on them.
    # For predators (indices S+1:n): degree is the number of prey they feed on.
    S = A_p.S
    R = A_p.R
    P_matrix = A_p.P_matrix
    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha])
    end
     
    return (sensitivity = sensitivity, biomass = biomass, degree = degree)
end

# --- Function to plot sensitivity vs species traits using Makie ---
function plot_sensitivity_traits(metrics)
    # Create a figure with two panels:
    #   - One scatter plot: sensitivity vs equilibrium biomass.
    #   - Another scatter plot: sensitivity vs degree.
    fig = Figure(resolution = (900, 400))
    ax1 = Axis(fig[1,1], xlabel = "Equilibrium Biomass", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Equilibrium Biomass")
    MK.scatter!(ax1, metrics.biomass, metrics.sensitivity, markersize = 10)
    
    ax2 = Axis(fig[1,2], xlabel = "Degree (Number of Interactions)", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Degree")
    MK.scatter!(ax2, metrics.degree, metrics.sensitivity, markersize = 10)
    
    display(fig)
end

# --- Example Usage ---
# Assuming you have already obtained the equilibrium and parameters from your analytical_equilibrium function:
# For example:
# A_eq, A_p, A_j = analytical_equilibrium( ... )

# (Here we assume A_eq and A_p are already defined in your workflow.)
metrics = compute_sensitivity_metrics(A_eq, A_p; perturbation = 0.01, tspan = (0.0, 50.0))
fig_traits = plot_sensitivity_traits(metrics)
display(fig_traits)

