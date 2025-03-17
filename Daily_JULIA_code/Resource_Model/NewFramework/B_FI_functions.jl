function b_new_parametrise_the_community(
    species_names::Vector{String};
    mu::Float64 = 0.5,
    epsilon_val::Float64 = 1.0,
    mean_m_alpha::Float64 = 0.1,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    cell_abundance::Vector{Float64} = Float64[],
    delta_nu::Float64 = 0.05,
    d_alpha::Float64 = 1.0,
    d_i::Float64 = 1.0
)
    # Identify herbivores and predators in the cell.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)
    
    # For herbivores, set baseline abundances H_i0 (from observed data).
    if S > 0
        if length(cell_abundance) == S
            H_i0 = copy(cell_abundance)
        else
            println("WARNING: cell_abundance length != S, defaulting to uniform abundances.")
            H_i0 = ones(S)
        end
    else
        H_i0 = Float64[]
    end
    
    # Set herbivore intrinsic growth rates. For simplicity, we assume r_i = 1 for all.
    r_i = ones(S)
    
    # Define carrying capacities for herbivores.
    # Here we set K_i = H_i0 (i.e. the observed data are our equilibrium herbivore abundances).
    K_i = copy(H_i0)
    
    # For predators, assign mortality and conversion efficiency.
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        epsilon = fill(epsilon_val, R) * (d_i / d_alpha)
        # We assume predator self-regulation factor d_α = 1 so that K_α = m_α.
        K_alpha = copy(m_alpha)
    else
        m_alpha = Float64[]
        epsilon = Float64[]
        K_alpha = Float64[]
    end
    
    # Build the predation incidence matrix (S x R).
    # P_matrix[i, α] = 1 if predator α preys on herbivore i, 0 otherwise.
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
    
    # --- Equilibrium Conditions for Predators ---
    # For each predator α, we want:
    #    P*_α = ε * ν * (sum over i of P_matrix[i,α]*H_i0) - m_α > 0.
    # This implies: ν > m_α / (ε * (sum over i of P_matrix[i,α]*H_i0)).
    nu_candidates = Float64[]
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_i0)
        if H_alpha_tot > 0
            push!(nu_candidates, m_alpha[alpha] / (epsilon[alpha] * H_alpha_tot))
        end
    end
    if !isempty(nu_candidates)
        nu_min = maximum(nu_candidates)
    else
        nu_min = 0.0
    end
    
    # Compute equilibrium predator abundances using the computed nu_min.
    # For each predator α: P*_α = ε * nu_min * (sum over i of P_matrix[i,α]*H_i0) - m_α.
    P_star = zeros(R)
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_i0)
        P_star[alpha] = epsilon[alpha] * nu_min * H_alpha_tot - m_alpha[alpha]
    end
    if any(P_star .<= 0.0)
        for alpha in 1:R
            H_alpha_tot = sum(P_matrix[:, alpha] .* H_i0)
            P_star[alpha] = epsilon[alpha] * nu_min*delta_nu * H_alpha_tot - m_alpha[alpha]
        end
        if any(P_star .<= 0.0)
            @info "Some predator equilibrium abundances are still negative or zero"
            for alpha in 1:R 
                H_alpha_tot = sum(P_matrix[:, alpha] .* H_i0)
                P_star[alpha] = epsilon[alpha] * nu_min * H_alpha_tot - m_alpha[alpha]
            end
        else
            @info "We had to use a higher nu_min (1+delta, where delta = $delta_nu) to avoid negative predator abundances"
        end
        
    end
    
    # Assume herbivore equilibrium equals the observed data.
    H_star = copy(H_i0)
    
    return (
        S = S, R = R,
        H_i0 = H_i0, r_i = r_i, K_i = K_i,
        mu = mu, nu = nu_min,
        P_matrix = P_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_star
    )
end

function b_new_setup_community_from_cell(
    i::Int, j::Int;
    mu::Float64 = 0.5,
    # nu::Float64 = 0.01,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    species_names::Vector{String} = String[],
    artificial_pi::Bool = false,
    delta_nu::Float64 = 0.05,
    d_alpha::Float64 = 1.0,
    d_i::Float64 = 1.0
)
    # 1) Retrieve the cell and extract species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    
    if isempty(species_names)
        species_names = o_extract_species_names_from_a_cell(cell)
    end
    S, R = o_identify_n_of_herbs_and_preds(species_names)
    
    # 2) Build herbivore abundance vector.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    
    cell_abundance_herbs = Float64[]
    for sp in herbivore_list
        local_idx = species_dict_herbivores_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_herbs, val)
        end
    end
    if artificial_pi
        cell_abundance_herbs = ones(length(cell_abundance_herbs))*100.0
    end
    
    # 3) Call the new parameterisation function for the new framework.
    params = b_new_parametrise_the_community(
         species_names;
         mu = mu,
         epsilon_val = epsilon_val,
         mean_m_alpha = mean_m_alpha,
         iberian_interact_NA = iberian_interact_NA,
         species_dict = species_dict,
         cell_abundance = cell_abundance_herbs,
         delta_nu = delta_nu,
         d_alpha = d_alpha,
         d_i = d_i
    )

    S, R, H_i0, r_i, K_i, mu, nu,
    P_matrix, epsilon, m_alpha, K_alpha,
    herbivore_list, predator_list, species_names,
    H_star, P_star = params
    
    return (
        S = S, R = R,
        H_i0 = H_i0, r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        P_matrix = P_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_star
    )
end

# Example usage:
A = b_new_setup_community_from_cell(
    18, 1;
    mu = 0.5,
    # nu = 0.01,
    mean_m_alpha = 0.1,
    epsilon_val = 1.0,
    iberian_interact_NA = iberian_interact_NA,
    species_dict = species_dict,
    # species_names = species_names,
    artificial_pi = false
)

S, R, H_i0, r_i, K_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha, herbivore_list, predator_list, species_names, H_star, P_star = A

