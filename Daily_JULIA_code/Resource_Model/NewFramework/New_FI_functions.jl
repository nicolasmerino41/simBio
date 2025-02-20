# FI_functions.jl – Revised Version Incorporating New Parametrization

################################################################################
# Basic Functions for Extracting Species Information
################################################################################
function extract_species_names_from_a_cell(cell::MyBirmmals)
    names = birmmals_biomass_fixed[:, :species]
    species_names = String[]
    for i in 1:length(cell.a)
        if !iszero(cell.a[i])
            push!(species_names, names[i])
        end
    end
    return species_names
end

function identify_n_of_herbs_and_preds(species_names::Vector{String})
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

function check_predator_has_prey(species_names::Vector{String})
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

################################################################################
# 1) PARAMETRISE THE COMMUNITY (Revised for New Model)
################################################################################
function new_parametrise_the_community(
    species_names::Vector{String};
    # Basic arguments:
    NPP::Float64 = 1000.0,
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    symmetrical_competition::Bool = true,
    
    # Predator parameters:
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    
    # Global NamedMatrix & dictionary for who-eats-whom
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    
    # Standard deviations for random draws (for m_i if not empirical)
    m_standard_deviation::Float64 = 0.0,
    h_standard_deviation::Float64 = 0.0,
    
    # Cell-level herbivore abundances (empirical values, one per herbivore)
    cell_abundance::Vector{Float64} = Float64[],
    
    # Instead of real_H0, we now use the empirical abundances directly:
    real_H0::Bool = true,
    
    # Optional: provide an empirical H0 vector (if available) otherwise, cell_abundance is used.
    H0_vector::Vector{Float64} = nothing,
    
    # New parameter: allometric exponent
    alpha::Float64 = 0.25
)
    #### (A) Identify herbivores vs. predators
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)
    
    #### (B) Define H_i^0 and m_i for herbivores
    # In the new framework, we set H_i^0 = H_i^obs, which are the empirical abundances.
    if S > 0
        if length(cell_abundance) == S
            H_i0 = copy(cell_abundance)
        else
            println("WARNING: cell_abundance length != S, defaulting to uniform abundances.")
            H_i0 = ones(S)
        end
        # For m_i, we can still use random draws (or empirical values if available).
        m_i = [abs(rand(Normal(M_mean, m_standard_deviation))) for _ in 1:S]
    else
        H_i0 = Float64[]
        m_i  = Float64[]
    end
    
    #### (C) Compute baseline competition among herbivores: generate mu_matrix
    if S > 0
        V, mu_matrix = generate_competition_matrix(S, mu, symmetrical_competition; check_condition=true)
    else
        V = zeros(0,0)
        mu_matrix = zeros(0,0)
    end
    
    #### (D) Predator-related matrices (a_matrix and A)
    a_matrix = zeros(S, R)
    A = Matrix{Float64}(I, R, R)
    for row_pred in 1:R
        A[row_pred, row_pred] = -1.0
    end
    for herb_i in 1:S
        global_herb_idx = species_dict[herbivore_list[herb_i]]
        for pred_j in 1:R
            global_pred_idx = species_dict[predator_list[pred_j]]
            if iberian_interact_NA[global_pred_idx, global_herb_idx] == 1
                a_matrix[herb_i, pred_j] = mu_predation
            end
        end
    end
    for row_pred in 1:R
        global_pred_row = species_dict[predator_list[row_pred]]
        for col_pred in 1:R
            if row_pred != col_pred
                global_pred_col = species_dict[predator_list[col_pred]]
                if iberian_interact_NA[global_pred_row, global_pred_col] == 1
                    A[row_pred, col_pred] = mu_predation
                end
            end
        end
    end
    m_alpha = fill(mean_m_alpha, R)
    epsilon = fill(epsilon_val, R)
    
    #### (E) Compute predator-mediated coefficients C and G, and the modified interaction matrix
    if R > 0
        A_inv = inv(A)
    else
        A_inv = zeros(0,0)
    end
    C = zeros(S,S)
    G = zeros(S)
    for r in 1:S
        for c in 1:S
            val = 0.0
            for alpha in 1:R
                for beta in 1:R
                    val += epsilon[alpha] * a_matrix[r, alpha] * A_inv[alpha, beta] * a_matrix[c, beta]
                end
            end
            C[r,c] = val
        end
    end
    for r in 1:S
        val = 0.0
        for alpha in 1:R
            for beta in 1:R
                val += a_matrix[r, alpha] * A_inv[alpha, beta] * m_alpha[beta]
            end
        end
        G[r] = val
    end
    M_modified = copy(mu_matrix)
    for r in 1:S, c in 1:S
        M_modified[r,c] += C[r,c] * (H_i0[r] / m_i[r])
    end
    
    #### (F) --- In the new framework, we use the empirical abundances as H_i0, so no need to compute hat_p.
    #### Instead, we now compute the intrinsic growth rate directly.
    
    # Retrieve body masses for herbivores from birmmals_biomass_fixed.
    bodyMasses = [birmmals_biomass_fixed[birmmals_biomass_fixed.species .== sp, :bodyMass][1] for sp in herbivore_list]
    
    # Compute x using the equilibrium condition:
    # x = NPP / sum_{i=1}^{S} (M_i^{-alpha} * H_i0[i])
    numerator_x = NPP
    denominator_x = sum([bodyMasses[i]^-alpha * H_i0[i] for i in 1:S])
    x = numerator_x / denominator_x
    
    # Compute raw growth rates:
    raw_g = [bodyMasses[i]^-alpha * x for i in 1:S]
    
    # Compute the niche parameter for each herbivore: beta_i = (g_i_raw / m_i) - 1.
    beta = [raw_g[i] / m_i[i] - 1 for i in 1:S]
    
    # Adjust the growth rate according to the new model:
    # g_i = raw_g[i] * (beta[i] / (1 + beta[i]))
    g_i = [raw_g[i] * (beta[i] / (1 + beta[i])) for i in 1:S]
    
    #### (G) --- The new model now does not compute p_vec or hat_p.
    #### Instead, we use the computed g_i in the herbivore dynamics.
    #### (The subsequent dynamical simulation will use these g_i values.)
    
    # (H) Return the computed parameters
    return (
        S, R,
        H_i0, m_i,
        g_i,          # final adjusted growth rates
        x,            # scaling parameter x
        beta,         # niche parameters
        G,            # predator feedback term
        M_modified,   # modified interaction matrix (includes predator effects)
        a_matrix, A,
        epsilon, m_alpha,
        raw_g         # raw growth rates before adjustment
    )
end


################################################################################
# 2) SETUP COMMUNITY FROM CELL
################################################################################
function new_setup_community_from_cell(
    i::Int, j::Int;
    NPP::Float64 = 1000.0,
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    m_standard_deviation::Float64 = 0.0,
    h_standard_deviation::Float64 = 0.0,
    artificial_pi = false,
    real_H0::Bool = false,
    H0_vector::Vector{Float64} = nothing,
    species_names::Vector{String} = nothing
)
    # 1) Grab the cell and determine species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    
    # 2) Build cell_abundance for herbivores (order must match herbivore_list)
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    
    cell_abundance_herbs = Float64[]
    herbivore_list_cell = String[]
    for sp in herbivore_list
        local_idx = species_dict_herbivores_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_herbs, val)
            push!(herbivore_list_cell, sp)
        end
    end
    
    if artificial_pi
        cell_abundance_herbs = ones(length(cell_abundance_herbs))
    end
    
    # 3) Parametrise the community using the updated model.
    # The new parametrisation function returns:
    #   S, R, H_i0, m_i, g_i, x_final, beta, G, M_modified, a_matrix, A, epsilon, m_alpha, raw_g
    S, R, H_i0, m_i, g_i, x_final, beta, G, M_modified, a_matrix, A, epsilon, m_alpha, raw_g =
    new_parametrise_the_community(
            species_names;
            NPP = NPP,
            M_mean = M_mean,
            mu = mu,
            symmetrical_competition = symmetrical_competition,
            mean_m_alpha = mean_m_alpha,
            epsilon_val = epsilon_val,
            mu_predation = mu_predation,
            iberian_interact_NA = iberian_interact_NA,
            species_dict = species_dict,
            m_standard_deviation = m_standard_deviation,
            h_standard_deviation = h_standard_deviation,
            cell_abundance = cell_abundance_herbs,
            real_H0 = real_H0,
            H0_vector = H0_vector
        )
    
    return (
        S, R,
        species_names, herbivore_list, predator_list,
        H_i0, m_i,
        g_i, x_final, beta,   # adjusted growth rates, scaling parameter, niche parameters
        G,
        M_modified, a_matrix, A,
        epsilon, m_alpha,
        raw_g
    )
end


################################################################################
# A FUNCTION TO EXTRACT INFO OF A SPECIES
################################################################################
function extract_info_of_a_species(species_name::String)
    if species_name in birmmals_biomass_fixed.species
        pi_value = birmmals_biomass_fixed[birmmals_biomass_fixed.species .== species_name, :biomass][1]
        rounded_pi_value = round(pi_value[1], digits=2)
        sorted_biomass = sort(birmmals_biomass_fixed[:, :biomass], rev=true)
        rank = findfirst(x -> x == pi_value, sorted_biomass)
    else
        error("Species $species_name not found in birmmals_biomass_fixed")
    end

    position_in_birmmals = findfirst(x -> x == species_name, birmmals_names)
    position_in_spain_names = findfirst(x -> x == species_name, spain_names)

    if species_name in herbivore_names_as_birmmals
        what = "herbivore"
    elseif species_name in predator_names_as_birmmals
        what = "predator"
    else
        error("Species $species_name is neither herbivore nor predator")
    end

    if species_name in herbivore_names_as_birmmals
        H0_value_vector = []
        for i in idx
            push!(H0_value_vector, H0_DA[i].a[species_dict_herbivores_in_birmmals[species_name]])
        end
        average_H0_value = mean(H0_value_vector)
        H0_value_sd = std(H0_value_vector)
    end

    println("$species_name is a $what with p_i $rounded_pi_value, which falls $rank/$(length(birmmals_biomass_fixed[:, :biomass]))")
    println("Its position in birmmals_names is $position_in_birmmals and spain_names $position_in_spain_names")

    if species_name in herbivore_names_as_birmmals
        println("Its average H0 value is $average_H0_value")
        println("Its standard deviation is $H0_value_sd")
    end
end