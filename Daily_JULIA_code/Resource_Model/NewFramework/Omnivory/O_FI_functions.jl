################################################################################
# Basic Functions for Extracting Species Information
################################################################################
function o_extract_species_names_from_a_cell(cell::MyBirmmals)
    names = birmmals_biomass_fixed[:, :species]
    species_names = String[]
    for i in 1:length(cell.a)
        if !iszero(cell.a[i])
            push!(species_names, names[i])
        end
    end
    return species_names
end

function o_identify_n_of_herbs_and_preds(species_names::Vector{String})
    S = 0
    R = 0
    for name in species_names
        if name in herbivore_names
            S += 1
        elseif name in predator_names
            R += 1
        end
    end
    return S, R
end

function o_check_predator_has_prey(species_names::Vector{String})
    cell_predators = [sp for sp in species_names if sp in predator_names]
    predator_has_prey = Dict{String,Bool}()
    for pred in cell_predators
        global_pred_idx = species_dict[pred]
        found_prey = false
        for candidate in species_names
            if candidate == pred
                continue
            end
            global_candidate_idx = species_dict[candidate]
            if iberian_interact_NA[global_pred_idx, global_candidate_idx] == 1
                found_prey = true
                break
            end
        end
        predator_has_prey[pred] = found_prey
    end
    allHavePrey = all(values(predator_has_prey))
    if allHavePrey
        return (true, 0, String[])
    else
        predators_without_prey = [pred for (pred, hasprey) in predator_has_prey if !hasprey]
        return (false, length(predators_without_prey), predators_without_prey)
    end
end

function o_parametrise_the_community(
    species_names::Vector{String};
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    cell_abundance_h::Vector{Float64} = Float64[],  # Observed (assumed equilibrium) herbivore/omnivore abundances
    cell_abundance_p::Vector{Float64} = Float64[],  # Observed carnivore equilibrium abundances
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 0.01
)
    
    # Convert parameters to plain Float64:
    mu             = to_float(mu)
    epsilon_val    = to_float(epsilon_val)
    mean_m_alpha   = to_float(mean_m_alpha)
    delta_nu       = to_float(delta_nu)
    d_alpha        = to_float(d_alpha)
    d_i            = to_float(d_i)
    r_omni_proportion = to_float(r_omni_proportion)
    
    # Identify herbivores and omnivores (treated with the herbivore equations) 
    # and carnivores (treated with the predator equations).
    herbivore_list = [sp for sp in species_names if sp in herbivore_names || sp in omnivore_names]
    predator_list  = [sp for sp in species_names if sp in carnivore_names]
    S = length(herbivore_list)
    R = length(predator_list)
    
    # Process herbivore/omnivore equilibrium abundances.
    if S > 0
        if length(cell_abundance_h) == S
            H_eq = copy(cell_abundance_h)
        else
            error("Expected cell_abundance_h to have length S=$S")
        end
    else
        H_eq = Float64[]
    end
    
    # Process carnivore equilibrium abundances.
    if R > 0
        if length(cell_abundance_p) == R
            P_eq = copy(cell_abundance_p)
        else
            error("Expected cell_abundance_p to have length R=$R")
        end
    else
        P_eq = Float64[]
    end
    
    # For herbivores/omnivores, the observed equilibrium is our H*.
    H_star = copy(H_eq)
    r_i = [herbivore_list[i] in omnivore_names ? r_omni_proportion : 1.0 for i in 1:S]
    K_i_initial = copy(H_eq)
    
    # For carnivores, assign mortality and conversion efficiency.
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Scale epsilon by (d_i/d_alpha)
        epsilon = fill(epsilon_val, R) * (d_i / d_alpha)
        K_alpha_initial = copy(m_alpha)
    else
        m_alpha = Float64[]
        epsilon = Float64[]
        K_alpha_initial = Float64[]
    end
    
    # --- Build interaction matrices using iberian_interact_NA ---
    # 1. P_matrix (S x R): predation interactions where carnivores consume herbivores.
    P_matrix = zeros(S, R)
    for i in 1:S
        global_h_idx = species_dict[herbivore_list[i]]
        for j in 1:R
            global_p_idx = species_dict[predator_list[j]]
            if iberian_interact_NA[global_p_idx, global_h_idx] == 1
                P_matrix[i, j] = 1.0
            end
        end
    end
    
    # 2. O_matrix (S x S): for omnivory; only rows for species in omnivore_names have nonzero entries.
    #    O_matrix[i,j] = 1 if omnivore species i consumes herbivore species j.
    O_matrix = zeros(S, S)
    for i in 1:S
        if herbivore_list[i] in omnivore_names
            global_consumer_idx = species_dict[herbivore_list[i]]
            for j in 1:S
                global_prey_idx = species_dict[herbivore_list[j]]
                if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                    O_matrix[i, j] = 1.0
                end
            end
        end
    end
    
    # 3. T_matrix (S x S): for herbivore-herbivore interactions representing consumption.
    #    T_matrix[i,j] = 1 if herbivore (or omnivore) species j consumes species i.
    T_matrix = zeros(S, S)
    for i in 1:S
        global_prey_idx = species_dict[herbivore_list[i]]
        for j in 1:S
            global_consumer_idx = species_dict[herbivore_list[j]]
            if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                T_matrix[i, j] = 1.0
            end
        end
    end
    
    # 4. B_matrix (R x R): predator-predator beneficial interactions.
    #    B_matrix[α,β] = 1 if carnivore α consumes carnivore β.
    B_matrix = zeros(R, R)
    for alpha in 1:R
        global_consumer_idx = species_dict[predator_list[alpha]]
        for beta in 1:R
            global_prey_idx = species_dict[predator_list[beta]]
            if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                B_matrix[alpha, beta] = 1.0
            end
        end
    end
    
    # 5. D_matrix (R x R): predator-predator detrimental interactions.
    #    D_matrix[α,β] = 1 if carnivore β consumes carnivore α.
    D_matrix = zeros(R, R)
    for alpha in 1:R
        global_prey_idx = species_dict[predator_list[alpha]]
        for beta in 1:R
            global_consumer_idx = species_dict[predator_list[beta]]
            if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                D_matrix[alpha, beta] = 1.0
            end
        end
    end
    
    # --- Compute candidate ν values for each carnivore using the equilibrium condition ---
    # For each predator α, the equilibrium for the predator equation reads:
    #   εν*(H_alpha_tot + B_alpha_tot) - ν*(D_alpha_tot) = P_eq[α] + m_alpha[α],
    # where
    #   H_alpha_tot = ∑₍i₎ P_matrix[i,α]*H_star[i],
    #   B_alpha_tot = ∑₍β₎ B_matrix[α,β]*P_eq[β],
    #   D_alpha_tot = ∑₍β₎ D_matrix[α,β]*P_eq[β].
    nu_candidates = Float64[]
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_star)
        B_alpha_tot = sum(B_matrix[alpha, :] .* P_eq)
        D_alpha_tot = sum(D_matrix[alpha, :] .* P_eq)
        denominator = epsilon[alpha] * (H_alpha_tot + B_alpha_tot) - D_alpha_tot
        if denominator > 0
            push!(nu_candidates, (P_eq[alpha] + m_alpha[alpha]) / denominator)
        end
    end
    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + delta_nu)  # Add safety margin.
    
    # --- Recalculate derived carrying capacities ---
    # For herbivores: 
    #   K_i = (1-μ)*H*_i + μ*(∑₍j₎ H*_j) + ν*(P_i^tot) - εν*(O_i^tot) + ν*(T_i^tot),
    # where for each herbivore i:
    #   P_i^tot = ∑₍predators₎ [if iberian_interact_NA indicates predation on i] * P_eq,
    #   O_i^tot = ∑₍j₎ [if species i (an omnivore) consumes j] * H_star[j] (0 if i is not omnivore),
    #   T_i^tot = ∑₍j₎ [if species j consumes i] * H_star[j].
    K_i = zeros(S)
    total_H = sum(H_star)
    for i in 1:S
        # Compute predation pressure from carnivores.
        P_i_tot = sum(P_matrix[i, :] .* P_eq)
        # For omnivory: if species i is an omnivore, sum the biomass of its prey.
        O_i_tot = 0.0
        if herbivore_list[i] in omnivore_names
            O_i_tot = sum(O_matrix[i, :] .* H_star)
        end
        # Trophic loss from being consumed by other herbivores.
        T_i_tot = sum(T_matrix[i, :] .* H_star)
        K_i[i] = (1 - mu) * H_star[i] + mu * total_H + nu * P_i_tot - epsilon_val * nu * O_i_tot + nu * T_i_tot
    end
    
    # For carnivores, the carrying capacity is given by:
    #   K_α = εν*(H_alpha_tot) - P_eq[α] + εν*(B_alpha_tot) - ν*(D_alpha_tot).
    K_alpha = zeros(R)
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_star)
        B_alpha_tot = sum(B_matrix[alpha, :] .* P_eq)
        D_alpha_tot = sum(D_matrix[alpha, :] .* P_eq)
        K_alpha[alpha] = epsilon[alpha] * nu * H_alpha_tot - P_eq[alpha] + epsilon[alpha] * nu * B_alpha_tot - nu * D_alpha_tot
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
        O_matrix = O_matrix,
        T_matrix = T_matrix,
        B_matrix = B_matrix,
        D_matrix = D_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_eq  # Final equilibrium values are taken from inputs.
    )
end

function o_setup_community_from_cell(
    i::Int, j::Int;
    mu = 0.5,
    # nu = 0.01,
    mean_m_alpha = 0.1,
    epsilon_val = 1.0,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    species_names::Vector{String} = String[],
    artificial_pi::Bool = false,
    pi_size = 1.0,
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 0.01
)
    # 1) Retrieve the cell and extract species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    
    if isempty(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    
    # 2) Build herbivore (and omnivore) abundance vector.
    # Note: species in herbivore_names OR omnivore_names are treated via the herbivore equations.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names || sp in omnivore_names]
    
    cell_abundance_herbs = Float64[]
    for sp in herbivore_list
        local_idx = species_dict_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_herbs, val)
        end
    end
    if artificial_pi
        cell_abundance_herbs = [sp in herbivore_names ? pi_size : sp in omnivore_names ? pi_size * 0.5 : pi_size * 0.1 for sp in herbivore_list]
    end
    
    # 3) Build carnivore (predator) abundance vector.
    predator_list  = [sp for sp in species_names if sp in carnivore_names]
    
    cell_abundance_preds = Float64[]
    for sp in predator_list
        local_idx = species_dict_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_preds, val)
        end
    end
    if artificial_pi
        cell_abundance_preds = fill(pi_size * 0.1, length(cell_abundance_preds))
    end
    
    # 4) Call the new parameterisation function for the new framework.
    params = o_parametrise_the_community(
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
         d_i = d_i,
         r_omni_proportion = r_omni_proportion
    )

    S, R, H_eq, P_eq,
    r_i, K_i, mu, nu,
    P_matrix, O_matrix, T_matrix, B_matrix, D_matrix,
    epsilon, m_alpha, K_alpha,
    herbivore_list, predator_list, species_names,
    H_star, P_star = params
    
    return (
        S = S, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        P_matrix = P_matrix, O_matrix = O_matrix, T_matrix = T_matrix, B_matrix = B_matrix, D_matrix = D_matrix,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names,
        H_star = H_star, P_star = P_star
    )
end

# A_s = o_setup_community_from_cell(
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
