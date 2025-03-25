function b_parametrise_the_community(
    species_names::Vector{String};
    mu::Float64 = 0.5,
    epsilon_val::Float64 = 1.0,
    mean_m_alpha::Float64 = 0.1,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    cell_abundance_h::Vector{Float64} = Float64[],  # Observed (assumed equilibrium) herbivore abundances
    cell_abundance_p::Vector{Float64} = Float64[],  # Observed predator equilibrium abundances
    delta_nu::Float64 = 0.05,
    d_alpha::Float64 = 1.0,
    d_i::Float64 = 1.0
)
    # Identify herbivores and predators in the cell.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
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
    new_di = r_i ./ K_i
    new_da = m_alpha ./ K_alpha
    # println("new_di = ", new_di)
    # println("new_da = ", new_da)
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

function b_setup_community_from_cell(
    i::Int, j::Int;
    mu::Float64 = 0.5,
    # nu::Float64 = 0.01,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    species_names::Vector{String} = String[],
    artificial_pi::Bool = false,
    pi_size::Float64 = 1.0,
    delta_nu::Float64 = 0.05,
    d_alpha::Float64 = 1.0,
    d_i::Float64 = 1.0
)
    # 1) Retrieve the cell and extract species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    
    if isempty(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    
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
        cell_abundance_herbs = fill(pi_size, length(cell_abundance_herbs))
    end
    
    # 3) Build predator abundance vector.
    predator_list  = [sp for sp in species_names if sp in predator_names]
    
    cell_abundance_preds = Float64[]
    for sp in predator_list
        local_idx = species_dict_predators_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_preds, val)
        end
    end
    if artificial_pi
        cell_abundance_preds = fill(pi_size*0.1, length(cell_abundance_preds))
    end
    
    # 3) Call the new parameterisation function for the new framework.
    params = b_parametrise_the_community(
         species_names;
         mu = mu,
         epsilon_val = epsilon_val,
         mean_m_alpha = mean_m_alpha,
         iberian_interact_NA = iberian_interact_NA,
         species_dict = species_dict,
         cell_abundance_h = cell_abundance_herbs,
         cell_abundance_p = cell_abundance_preds,
         delta_nu = delta_nu,
         d_alpha = d_alpha,
         d_i = d_i
    )    

    S, R, H_eq, P_eq,
    r_i, K_i, mu, nu,
    P_matrix, epsilon, m_alpha, K_alpha,
    herbivore_list, predator_list, species_names,
    H_star, P_star = params
    
    return (
        S = S, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        P_matrix = P_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_star
    )
end

# # Example usage:
# A = b_setup_community_from_cell(
#     18, 1;
#     mu = 0.5,
#     # nu = 0.01,
#     mean_m_alpha = 0.1,
#     epsilon_val = 1.0,
#     iberian_interact_NA = iberian_interact_NA,
#     species_dict = species_dict,
#     # species_names = species_names,
#     artificial_pi = true,
#     delta_nu = 0.05
# )

# S, R, H_i0, r_i, K_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha, herbivore_list, predator_list, species_names, H_star, P_star = A

# b_attempt_setup_community(
#     18, 1, 
#     0.5, 1.0, 0.1;
#     species_names = nothing,
#     artificial_pi = true, pi_size = 1.0,
#     delta_nu = 0.05,
#     d_alpha = 1.0, d_i = 1.0
# )