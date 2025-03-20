function random_parametrise_the_community(
    species_names::Vector{String};
    mu::Float64 = 0.5,
    epsilon_val::Float64 = 1.0,
    mean_m_alpha::Float64 = 0.1,
    iberian_interact_NA::Matrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    cell_abundance_h::Vector{Float64} = Float64[],  # Observed (assumed equilibrium) herbivore abundances
    cell_abundance_p::Vector{Float64} = Float64[],  # Observed predator equilibrium abundances
    delta_nu::Float64 = 0.05,
    d_alpha::Float64 = 1.0,
    d_i::Float64 = 1.0
)
    # Identify herbivores and predators.
    herbivore_list = [sp for sp in species_names if startswith(sp, "H")]
    predator_list  = [sp for sp in species_names if startswith(sp, "P")]
    S = length(herbivore_list)
    R = length(predator_list)
    
    # Process herbivore equilibrium abundances.
    if S > 0
        if length(cell_abundance_h) == S
            H_eq = copy(cell_abundance_h)
        else
            error("Expected cell_abundance_h to have length S=$S")
        end
    else
        H_eq = Float64[]
    end
    
    # Process predator equilibrium abundances.
    if R > 0
        if length(cell_abundance_p) == R
            P_eq = copy(cell_abundance_p)
        else
            error("Expected cell_abundance_p to have length R=$R")
        end
    else
        P_eq = Float64[]
    end
    
    # For herbs, let the observed equilibrium be our H*.
    H_star = copy(H_eq)
    r_i = ones(S)
    # Initially we set a preliminary K_i as the observed abundances.
    K_i_initial = copy(H_eq)
    
    # For predators, assign mortality and conversion efficiency.
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Scale epsilon by (d_i/d_alpha); here d_i and d_alpha are 1 by default.
        epsilon = fill(epsilon_val, R) * (d_i / d_alpha)
        # Initially, set a preliminary K_alpha equal to m_alpha.
        K_alpha_initial = copy(m_alpha)
    else
        m_alpha = Float64[]
        epsilon = Float64[]
        K_alpha_initial = Float64[]
    end
    
    # Build the predation incidence matrix (S x R) based on the provided iberian_interact_NA.
    P_matrix = zeros(S, R)
    for i in 1:S
        global_herb_idx = species_dict[herbivore_list[i]]
        for j in 1:R
            global_pred_idx = species_dict[predator_list[j]]
            if iberian_interact_NA[global_pred_idx, global_herb_idx] == 1
                P_matrix[i, j] = 1.0
            end
        end
    end
    
    # --- Compute candidate ν values for each predator using the equilibrium condition ---
    # For predator α: assume P_eq[α] = ε[α]*ν*(sum_i P_matrix[i,α]*H_eq[i]) - m_alpha[α].
    # Rearranged:
    #   ν_α = (P_eq[α] + m_alpha[α]) / (ε[α] * (sum_i P_matrix[i,α]*H_eq[i]))
    nu_candidates = Float64[]
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_eq)
        if H_alpha_tot > 0
            push!(nu_candidates, (P_eq[alpha] + m_alpha[alpha]) / (epsilon[alpha] * H_alpha_tot))
        end
    end
    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + delta_nu)  # Add safety margin.
    
    # --- Recalculate derived carrying capacities ---
    # For herbivores: K_i = (1-μ)*H*_i + μ*(sum_j H*_j) + ν*(P_i^tot),
    # where P_i^tot = sum over predators that affect herbivore i.
    K_i = zeros(S)
    total_H = sum(H_star)
    for i in 1:S
        P_i_tot = sum(P_matrix[i, :] .* P_eq)
        K_i[i] = (1 - mu)*H_star[i] + mu * total_H + nu * P_i_tot
    end
    
    # For predators: K_α = ε*ν*(sum_i P_matrix[i,α]*H*_i) - P_eq[α].
    K_alpha = zeros(R)
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_star)
        K_alpha[alpha] = epsilon[alpha] * nu * H_alpha_tot - P_eq[alpha]
    end
    
    return (
        S = S, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        P_matrix = P_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_eq  # Final equilibrium values are taken from inputs.
    )
end

# -- Main Function: Random Bipartite Run --
function random_bipartite_run(
    S::Int, R::Int;
    mu::Float64 = 0.5,
    epsilon_val::Float64 = 1.0,
    mean_m_alpha::Float64 = 0.1,
    connectance::Float64 = 0.2,
    d_i::Union{Float64, Vector{Float64}} = 1.0,
    d_alpha::Union{Float64, Vector{Float64}} = 1.0,
    artificial_pi::Bool = false,
    time_end::Float64 = 500.0,
    delta_nu::Float64 = 0.05,
    abstol::Float64 = 1e-8,
    reltol::Float64 = 1e-6,
    plot::Bool = true,
    do_you_want_params::Bool = false,
    do_you_want_sol::Bool = false,
    min_h_abundance::Float64 = 5.0, max_h_abundance::Float64 = 10.0,
    min_p_abundance::Float64 = 0.5, max_p_abundance::Float64 = 1.0
)
    # Generate random herbivore names and predator names.
    herbivore_names = ["H$(i)" for i in 1:S]
    predator_names  = ["P$(i)" for i in 1:R]
    species_names = vcat(herbivore_names, predator_names)
    
    # Build species dictionary (herbivores first, then predators).
    species_dict = Dict{String, Int}()
    for (i, sp) in enumerate(species_names)
        species_dict[sp] = i
    end
    
    Random.seed!(123)
    # Generate random herbivore abundances (observed data).
    # For instance, uniform between 0.1 and 10.
    H_i0 = [rand(Uniform(min_h_abundance, max_h_abundance)) for _ in 1:S]
    P_i0 = [rand(Uniform(min_p_abundance, max_p_abundance)) for _ in 1:R]
    if artificial_pi
        H_i0 .= 1.0
        P_i0 .= 0.1
    end
    
    # Generate a random iberian_interact matrix.
    # Dimensions: (S+R) x (S+R). For predator–herbivore interactions,
    # for row indices S+1:S+R and column indices 1:S, set entry = 1 with probability connectance.
    N = S + R
    iberian_interact = zeros(Float64, N, N)
    for j in 1:R  # predators are in rows S+1 to S+R
        global_pred_idx = S + j
        for i in 1:S  # herbivores are in columns 1:S
            if rand() < connectance
                iberian_interact[global_pred_idx, i] = 1.0
            end
        end
    end
    # (Other entries can remain zero; they are not used in the parameterisation.)
    
    # Call the modified parameterisation function.
    params = random_parametrise_the_community(
        species_names;
        mu = mu,
        epsilon_val = epsilon_val,
        mean_m_alpha = mean_m_alpha,
        d_i = d_i,
        d_alpha = d_alpha,
        iberian_interact_NA = iberian_interact,
        species_dict = species_dict,
        cell_abundance_h = H_i0,
        cell_abundance_p = P_i0,
        delta_nu = delta_nu
    )
    println("Parameters: ", params)
    
    # Extract returned parameters.
    S_param = params.S
    R_param = params.R
    H_star = params.H_star   # equilibrium herbivore abundances (set equal to H_i0)
    P_star = params.P_star   # computed predator equilibria
    r_i = params.r_i
    K_i = params.K_i
    nu = params.nu
    # (Other parameters: mu, P_matrix, epsilon, m_alpha, K_alpha)
    P_matrix = params.P_matrix
    epsilon = params.epsilon
    m_alpha = params.m_alpha
    K_alpha = params.K_alpha
    
    # For simulation, we set initial conditions to the computed equilibrium values.
    u0 = vcat(H_star, P_star)
    
    # Build parameter tuple for the ODE function.
    ode_params = (S_param, R_param, K_i, r_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha)
    
    # Define and solve the ODE.
    cb_no_trigger, cb_trigger = build_callbacks(S_param, R_param, EXTINCTION_THRESHOLD, time_end, 1)
    prob = ODEProblem(bipartite_dynamics!, u0, (0.0, time_end), ode_params)
    sol = solve(prob, Tsit5(); callback=cb_no_trigger, abstol=abstol, reltol=reltol)
    # println("sol = ", sol)
    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (600, 500))
        ax = Axis(
            fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for random cell"
        )
        times_combined = sol.t
        # Plot herbivore dynamics (indices 1:S) in blue.
        for i in 1:S
            lines!(ax, times_combined, sol[i, :], label = "H$(i)", color = :blue)
        end
        # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
        for alpha in 1:R
            lines!(ax, times_combined, sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
        end
        display(fig)
    end

    H_end = sol[1:S, end]
    P_end = sol[S+1:S+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    survival_rate = total_surv / total_species
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    pred_herb_ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))

    df = DataFrame(
        S = [S_param],
        survived_herb = [survived_herb],
        R = [R_param],
        survived_pred = [survived_pred],
        survival_rate = [survival_rate],
        H_biomass = [H_biomass],
        P_biomass = [P_biomass],
        pred_herb_ratio = [pred_herb_ratio],
        mu = [mu],
        nu = [nu],
        epsilon = [epsilon_val],
        mean_m_alpha = [mean_m_alpha],
        connectance = [connectance],
        d_i = [d_i],
        d_alpha = [d_alpha],
        H_equilibrium = [H_star],
        P_equilibrium = [P_star],
        H_final = [H_end],
        P_final = [P_end],
    )

    if do_you_want_params && do_you_want_sol
        return df, params, sol
    elseif do_you_want_params || do_you_want_sol
        return df, do_you_want_params ? params : sol
    else
        return df
    end
    
end

# Example usage:
A = random_bipartite_run(
    10, 10;
    mu = 0.9,
    epsilon_val = 10.0,
    mean_m_alpha = 0.1,
    connectance = 1.0,
    d_i = 1.0,
    d_alpha = 1.0,   # weaker self-regulation for predators (larger K_alpha)
    artificial_pi = false,
    min_h_abundance = 5.0, max_h_abundance = 10.0,
    time_end = 500.0,
    delta_nu = 0.05,
    plot = true,
    do_you_want_params = false,
    do_you_want_sol = false
)
