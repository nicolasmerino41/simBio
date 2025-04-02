################################################################################
# Basic Functions for Extracting Species Information
################################################################################
function a_parametrise_the_community(
    S::Int, O::Int, R::Int;
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    conn = 0.2,  # target connectance (probability of an interaction)
    cell_abundance_h::Vector{Float64} = Float64[],  # Observed (or assumed) herbivore/omnivore abundances
    cell_abundance_p::Vector{Float64} = Float64[],  # Observed carnivore abundances
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 1.0
)
    # Convert parameters to Float64:
    mu             = Float64(mu)
    epsilon_val    = Float64(epsilon_val)
    mean_m_alpha   = Float64(mean_m_alpha)
    delta_nu       = Float64(delta_nu)
    d_alpha        = Float64(d_alpha)
    d_i            = Float64(d_i)
    r_omni_proportion = Float64(r_omni_proportion)
    conn           = Float64(conn)

    # Total number of "herbivores" is the sum of true herbivores and omnivores.
    S_total = S + O

    # Define species lists.
    herbivore_list = vcat(["Herbivore $sp" for sp in 1:S],
                          ["Omnivore $sp" for sp in 1:O])
    predator_list  = ["Predator $sp" for sp in 1:R]
    sp_nm = vcat(herbivore_list, predator_list)  # synthetic species names

    # Process equilibrium abundances.
    if S_total > 0
        if length(cell_abundance_h) == S_total
            H_eq = copy(cell_abundance_h)
        else
            # If not provided, default to ones.
            H_eq = ones(S_total)
        end
    else
        H_eq = Float64[]
    end

    if R > 0
        if length(cell_abundance_p) == R
            P_eq = copy(cell_abundance_p)
        else
            # If not provided, default to ones.
            P_eq = ones(R)
        end
    else
        P_eq = Float64[]
    end

    # For herbivores/omnivores, H* is H_eq.
    H_star = copy(H_eq)
    # Assume intrinsic growth rate = 1 for each herbivore.
    r_i = ones(S_total)
    # Initial carrying capacities for herbivores: use observed abundances.
    K_i_initial = copy(H_eq)

    # For predators, assign mortality and conversion efficiency.
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Scale epsilon by (d_i/d_alpha).
        epsilon = fill(epsilon_val, R) * (d_i / d_alpha)
        K_alpha_initial = copy(m_alpha)
    else
        m_alpha = Float64[]
        epsilon = Float64[]
        K_alpha_initial = Float64[]
    end

    # --- Build synthetic interaction matrices based on connectance ---
    # 1. P_matrix (S_total x R): predation interactions (carnivores preying on herbivores).
    P_matrix = zeros(Float64, S_total, R)
    for i in 1:S_total
        for j in 1:R
            P_matrix[i, j] = rand() < conn ? 1.0 : 0.0
        end
    end

    # 2. O_matrix (S_total x S_total): omnivory interactions.
    # Only rows corresponding to omnivores (i > S) have nonzero entries.
    O_matrix = zeros(Float64, S_total, S_total)
    for i in 1:S_total
        if i > S   # i.e. species from the Omnivore block
            for j in 1:S_total
                if i != j  # no self-interaction
                    O_matrix[i, j] = rand() < conn ? 1.0 : 0.0
                end
            end
        end
    end

    # 3. T_matrix (S_total x S_total): herbivore-herbivore (or omnivore-herbivore) interactions.
    # For simplicity, generate random interactions (excluding self-interactions).
    T_matrix = zeros(Float64, S_total, S_total)
    for i in 1:S_total
        for j in 1:S_total
            if i != j
                T_matrix[i, j] = rand() < conn ? 1.0 : 0.0
            end
        end
    end

    # 4. B_matrix (R x R): beneficial predator-predator interactions.
    B_matrix = zeros(Float64, R, R)
    for alpha in 1:R
        for beta in 1:R
            if alpha != beta
                B_matrix[alpha, beta] = rand() < conn ? 1.0 : 0.0
            end
        end
    end

    # 5. D_matrix (R x R): detrimental predator-predator interactions.
    D_matrix = zeros(Float64, R, R)
    for alpha in 1:R
        for beta in 1:R
            if alpha != beta
                D_matrix[alpha, beta] = rand() < conn ? 1.0 : 0.0
            end
        end
    end

    # --- Compute candidate ν values for each predator using the equilibrium condition ---
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
    nu *= (1.0 + delta_nu)  # apply safety margin

    # --- Recalculate derived carrying capacities ---
    # For herbivores:
    K_i = zeros(Float64, S_total)
    total_H = sum(H_star)
    for i in 1:S_total
        # Predation pressure from carnivores.
        P_i_tot = sum(P_matrix[i, :] .* P_eq)
        # Omnivory: only if species i is an omnivore.
        O_i_tot = (i > S) ? sum(O_matrix[i, :] .* H_star) : 0.0
        # Loss from being consumed by other herbivores.
        T_i_tot = sum(T_matrix[i, :] .* H_star)
        K_i[i] = (1 - mu) * H_star[i] + mu * total_H + nu * P_i_tot - epsilon_val * nu * O_i_tot + nu * T_i_tot
    end

    # For predators:
    K_alpha = zeros(Float64, R)
    for alpha in 1:R
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_star)
        B_alpha_tot = sum(B_matrix[alpha, :] .* P_eq)
        D_alpha_tot = sum(D_matrix[alpha, :] .* P_eq)
        K_alpha[alpha] = epsilon[alpha] * nu * H_alpha_tot - P_eq[alpha] + epsilon[alpha] * nu * B_alpha_tot - nu * D_alpha_tot
    end

    # (Optional) Compute derived ratios if needed:
    new_di = r_i ./ K_i
    new_da = m_alpha ./ K_alpha

    return (
        S = S_total, R = R,
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
        species_names = sp_nm,
        H_star = H_star, P_star = P_eq  # Final equilibrium values are taken from inputs.
    )
end


function a_setup_community_from_cell(
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
