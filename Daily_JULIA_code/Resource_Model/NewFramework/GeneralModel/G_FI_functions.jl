# G_FI_functions.jl
# This file contains functions to set up the community model for the new framework.
using Distributions  # for the Truncated Normal distribution

function g_parametrise_the_community(
    species_names::Vector{String};
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    cell_abundance_h::Vector{Float64} = Float64[],  # Herbivore/omnivore equilibrium abundances
    cell_abundance_p::Vector{Float64} = Float64[],  # Predator equilibrium abundances
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 1.0,
    randomise_attack_rates::Bool = false,
    attack_rate_sd::Float64 = 0.2
)
    # Convert parameters to Float64.
    mu          = Float64(mu)
    epsilon_val = Float64(epsilon_val)
    mean_m_alpha = Float64(mean_m_alpha)
    delta_nu    = Float64(delta_nu)
    d_alpha     = Float64(d_alpha)
    d_i         = Float64(d_i)
    r_omni_proportion = Float64(r_omni_proportion)

    # Identify species groups.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names || sp in omnivore_names]
    predator_list  = [sp for sp in species_names if sp in carnivore_names]
    S_total = length(herbivore_list)
    R = length(predator_list)
    sp_nm = vcat(herbivore_list, predator_list)

    # Process equilibrium abundances.
    if S_total > 0
        if length(cell_abundance_h) == S_total
            H_eq = copy(cell_abundance_h)
        else
            error("Expected cell_abundance_h to have length S_total=$S_total")
        end
    else
        H_eq = Float64[]
    end
    if R > 0
        if length(cell_abundance_p) == R
            P_eq = copy(cell_abundance_p)
        else
            error("Expected cell_abundance_p to have length R=$R")
        end
    else
        P_eq = Float64[]
    end
    
    H_star = copy(H_eq)
    # For herbivores, assign intrinsic growth rates (with omnivores getting a lower rate).
    r_i = [sp in omnivore_names ? r_omni_proportion : 1.0 for sp in herbivore_list]

    # For predators.
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Scale conversion efficiency by (d_i/d_alpha).
        epsilon = fill(epsilon_val, R) * (d_i / d_alpha)
    else
        m_alpha = Float64[]
        epsilon = Float64[]
    end

    # --- Build the merged interaction matrix A ---
    # A will be of size (S_total+R) x (S_total+R).
    # For each interaction (as indicated by the binary matrix iberian_interact_NA), assign:
    #   - If a predator consumes a herbivore: the herbivore (target) loses biomass (–1)
    #     and the predator (consumer) gains biomass (+ε).
    #   - If an omnivore consumes a herbivore: similar rules apply.
    #   - For predator–predator interactions (or herbivore consuming predator if indicated), use the same rule.
    n_total = S_total + R
    A = zeros(Float64, n_total, n_total)
    for i in 1:n_total
        for j in 1:n_total
            global_target = species_dict[sp_nm[i]]
            global_consumer = species_dict[sp_nm[j]]
            if iberian_interact_NA[global_consumer, global_target] == 1
                # Determine types.
                target_is_herb = (i <= S_total)
                consumer_is_pred = (j > S_total)
                consumer_is_omnivore = (j <= S_total) && (sp_nm[j] in omnivore_names)
                consumer_is_herb = (j <= S_total) && !(sp_nm[j] in omnivore_names)
                if target_is_herb && consumer_is_pred && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = 1.0         # Herbivore loses biomass.
                    A[j,i] = epsilon_val   # Predator gains biomass.
                elseif target_is_herb && consumer_is_omnivore && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0         # Herbivore loses biomass.
                    A[j,i] = -epsilon_val   # Omnivore gains biomass.
                elseif (!target_is_herb) && consumer_is_pred && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0
                    A[j,i] = epsilon_val
                elseif (!target_is_herb) && (consumer_is_herb || consumer_is_omnivore) && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0
                    A[j,i] = -epsilon_val
                end
            end
        end
    end

    # Scale each row of A: herbivores (rows 1:S_total) divide by d_i; predators (rows S_total+1:n_total) by d_alpha.
    for i in 1:S_total
        A[i, :] .= A[i, :] ./ d_i
    end
    for i in S_total+1:n_total
        A[i, :] .= A[i, :] ./ d_alpha
    end

    # --- Build the competition matrix C (for herbivores only) ---
    C = fill(mu, S_total, S_total)
    for i in 1:S_total
        C[i, i] = 0.0
    end
    C .= C ./ d_i

    # --- Compute the effective predation scaling factor ν ---
    # For each predator (indexed i = S_total+α), calculate:
    #   ν_α = (P_eq[α] + m_α) / (∑_{j=1}^{S_total+R} A[i,j] * X_j),
    # where X = [H_star; P_eq].
    nu_candidates = Float64[]
    X = vcat(H_star, P_eq)
    for alpha in 1:R
        idx_pred = S_total + alpha
        gain = sum(A[idx_pred, :] .* X)
        if gain > 0
            push!(nu_candidates, (P_eq[alpha] + m_alpha[alpha]) / gain)
        end
    end
    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + delta_nu)  # Apply safety margin.

    # --- Embed the calculated ν into A ---
    A .= nu .* A

    # --- Recalculate carrying capacities ---
    # For herbivores:
    K_i = zeros(Float64, S_total)
    total_H = sum(H_star)
    for i in 1:S_total
        K_i[i] = (1 - mu) * H_star[i] + sum(A[i, :] .* X) + sum(C[i, :] .* H_star)
    end
    # For predators:
    K_alpha = zeros(Float64, R)
    for alpha in 1:R
        idx_pred = S_total + alpha
        K_alpha[alpha] = sum(A[idx_pred, :] .* X) - P_eq[alpha]
    end

    new_di = r_i ./ K_i
    new_da = m_alpha ./ K_alpha

    # --- Optionally randomise the attack rates ---
    if randomise_attack_rates
        # Use a truncated Normal to ensure positive multipliers.
        multiplier_dist = Truncated(Normal(1.0, attack_rate_sd), 0, Inf)
        for i in axes(A,1), j in axes(A,2)
            A[i, j] *= rand(multiplier_dist)
        end
    end

    return (
        S = S_total, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        A = A,
        C = C,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = sp_nm,
        H_star = H_star, P_star = P_eq
    )
end

function g_setup_community_from_cell(
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
    r_omni_proportion = 1.0,
    randomise_attack_rates = false,
    attack_rate_sd = 0.2
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
    params = g_parametrise_the_community(
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
        r_omni_proportion = r_omni_proportion,
        randomise_attack_rates = randomise_attack_rates,
        attack_rate_sd = attack_rate_sd
    )

    # 5) Extract the parameters from the returned structure.
    S, R, H_eq, P_eq, r_i, K_i, mu, nu, A, C,
    epsilon, m_alpha, K_alpha, herbivore_list, predator_list, sp_nm,
    H_star, P_star = params
    
    return (
        S = S, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_i, K_i = K_i,
        mu = mu, nu = nu,
        A_matrix = A, C_matrix = C,
        epsilon = epsilon, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = sp_nm,
        H_star = H_star, P_star = P_star
    )
end

function g_attempt_setup_community(
    i, j, 
    mu_val, eps_val, mean_m_alpha;
    species_names = nothing,
    artificial_pi = false, pi_size = 1.0,
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    r_omni_proportion = 1.0,
    randomise_attack_rates = false,
    attack_rate_sd = 0.2
)
    # Convert parameters:
    mu_val       = to_float(mu_val)
    eps_val      = to_float(eps_val)
    mean_m_alpha = to_float(mean_m_alpha)
    pi_size      = to_float(pi_size)
    delta_nu   = to_float(delta_nu)
    d_alpha    = to_float(d_alpha)
    d_i        = to_float(d_i)
    r_omni_proportion = to_float(r_omni_proportion)
    
    try
        params = g_setup_community_from_cell(
            i, j;
            mu = mu_val,
            mean_m_alpha = mean_m_alpha,
            epsilon_val = eps_val,
            iberian_interact_NA = iberian_interact_NA,
            species_dict = species_dict,
            species_names = isnothing(species_names) ? String[] : species_names,
            artificial_pi = artificial_pi, pi_size = pi_size,
            delta_nu = delta_nu,
            d_alpha = d_alpha,
            d_i = d_i,
            r_omni_proportion = r_omni_proportion,
            randomise_attack_rates = randomise_attack_rates,
            attack_rate_sd = attack_rate_sd
        )
        return (
            S = params.S, R = params.R,
            H_eq = params.H_eq, P_eq = params.P_eq,
            r_i = params.r_i, K_i = params.K_i,
            mu = params.mu, nu = params.nu,
            A_matrix = params.A_matrix, C_matrix = params.C_matrix,
            epsilon = params.epsilon, m_alpha = params.m_alpha, K_alpha = params.K_alpha,
            herbivore_list = params.herbivore_list, predator_list = params.predator_list,
            species_names = params.species_names,
            H_star = params.H_star, P_star = params.P_star
        )
    catch e
        println("hey error in a_attempt_setup_community")
        println(e)
        return nothing
    end
end
