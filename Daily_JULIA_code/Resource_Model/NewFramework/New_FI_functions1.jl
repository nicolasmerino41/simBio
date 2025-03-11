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
    # (We now drop externally supplied predator abundances.)
    
    # New keyword arguments for Holling type II functional response:
    hollingII::Bool = false,
    h::Float64 = 0.1,
    
    # New parameter: allometric exponent
    alpha::Float64 = 0.25
)
    #### (A) Identify herbivores vs. predators
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)
    
    #### (B) Define H_i^0 and m_i for herbivores
    if S > 0
        if length(cell_abundance) == S
            H_i0 = copy(cell_abundance)
        else
            println("WARNING: cell_abundance length != S, defaulting to uniform abundances.")
            H_i0 = ones(S)
        end
        m_i = [abs(rand(Normal(M_mean, m_standard_deviation))) for _ in 1:S]
    else
        H_i0 = Float64[]
        m_i  = Float64[]
    end
    
    #### (C) Baseline competition: mu_matrix
    if S > 0
        V, mu_matrix = generate_competition_matrix(S, mu, symmetrical_competition; check_condition=true)
    else
        V = zeros(0,0)
        mu_matrix = zeros(0,0)
    end
    
    #### (D) Construct predator-related matrices (a_matrix and A)
    a_matrix = zeros(S, R)
    A = Matrix{Float64}(I, R, R)
    for row_pred in 1:R
        A[row_pred, row_pred] = 1.0
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
    
    #### (E) Compute predator-mediated coefficients C and G (using the raw a_matrix)
    if R > 0
        A_inv = inv(A)
    else
        A_inv = zeros(0,0)
    end
    C_raw = zeros(S, S)
    G_raw = zeros(S)
    for r in 1:S
        for c in 1:S
            val = 0.0
            for alpha_idx in 1:R
                for beta_idx in 1:R
                    val += epsilon[alpha_idx] * a_matrix[r, alpha_idx] *
                           A_inv[alpha_idx, beta_idx] * a_matrix[c, beta_idx]
                end
            end
            C_raw[r, c] = val
        end
    end
    for r in 1:S
        val = 0.0
        for alpha_idx in 1:R
            for beta_idx in 1:R
                val += a_matrix[r, alpha_idx] * A_inv[alpha_idx, beta_idx] * m_alpha[beta_idx]
            end
        end
        G_raw[r] = val
    end
    
    #### (F) Incorporate predator-mediated effects into the competition matrix:
    M_modified = copy(mu_matrix)
    for r in 1:S, c in 1:S
        M_modified[r, c] += C_raw[r, c] * (H_i0[r] / m_i[r])
    end
    
    #### (G) Compute intrinsic growth rates via allometric scaling
    bodyMasses = [birmmals_biomass_fixed[birmmals_biomass_fixed.species .== sp, :bodyMass][1] for sp in herbivore_list]
    norm_sq = sum((bodyMasses[i]^(-alpha) * H_i0[i])^2 for i in 1:S)
    x = (norm_sq > 0) ? NPP / norm_sq : (println("WARNING: norm_sq is zero."), 0.0)[2]
    raw_g = [ (NPP * bodyMasses[i]^(-2*alpha) * H_i0[i]) / norm_sq for i in 1:S ]
    beta = [ raw_g[i] / m_i[i] - 1 for i in 1:S ]
    g_i = [ raw_g[i] * (beta[i] / (1 + beta[i])) for i in 1:S ]
    
    #### (H) Compute A_star: nondimensional attack rates for herbivores
    d_i = [ m_i[i] / H_i0[i] for i in 1:S ]
    A_star = zeros(S, R)
    for i in 1:S
        for alph in 1:R
            A_star[i, alph] = a_matrix[i, alph] / d_i[i]
        end
    end

    #### (I) Compute predator contribution using the model derivation
    C_new = zeros(S, S)
    G_new = zeros(S)
    if R > 0
        for i in 1:S
            for j in 1:S
                val = 0.0
                for alph in 1:R
                    for bet in 1:R
                        val += A_star[i, alph] * A_inv[alph, bet] * a_matrix[j, bet]
                    end
                end
                C_new[i, j] = val
            end
            valG = 0.0
            for alph in 1:R
                for bet in 1:R
                    valG += A_star[i, alph] * A_inv[alph, bet] * m_alpha[bet]
                end
            end
            G_new[i] = valG
        end
    end

    #### (J) Define effective observed herbivore abundance incorporating predation
    H_i0_eff = zeros(S)
    for i in 1:S
        comp_sum = 0.0
        for j in 1:S
            if i != j
                comp_sum += mu_matrix[i, j] * H_i0[j]
            end
        end
        # Compute the linear predator contribution
        linear_pred_contrib = (R > 0) ? (sum(C_new[i, j] * H_i0[j] for j in 1:S) - G_new[i]) : 0.0
        # If Holling type II is requested, saturate the predator contribution:
        if hollingII
            pred_contrib = linear_pred_contrib / (1 + h * linear_pred_contrib)
        else
            pred_contrib = linear_pred_contrib
        end
        H_i0_eff[i] = H_i0[i] + comp_sum + pred_contrib
    end
    
    #### (K) Return the computed parameters (using H_i0_eff instead of H_i0)
    return (
        S, R,
        H_i0_eff, m_i,
        g_i,          # adjusted effective growth rates
        beta,         # niche parameters
        G_raw,        # predator feedback term (raw version)
        M_modified,   # modified competition matrix (including predator effects)
        A_star,       # nondimensional attack rates for herbivores
        a_matrix, A,
        epsilon, m_alpha,
        x,            # scaling parameter x
        raw_g,         # raw growth rates before adjustment
        h
    )
end

################################################################################
# 2) SETUP COMMUNITY FROM CELL (Revised for New Model)
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
    species_names::Vector{String} = nothing,
    alpha = 0.25,
    hollingII::Bool = false,
    h::Float64 = 0.1
)
    # 1) Retrieve the cell and extract species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    
    # 2) Build herbivore abundance vector from the cell.
    # (The order must match the herbivore list used in the parametrisation.)
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    
    cell_abundance_herbs = Float64[]
    for sp in herbivore_list
        local_idx = species_dict_herbivores_in_birmmals[sp]
        val = cell.a[local_idx]
        if val > 0.0
            push!(cell_abundance_herbs, val)
        end
    end
    if artificial_pi
        cell_abundance_herbs = ones(length(cell_abundance_herbs))
    end
    
    # 3) Call the new parametrisation function.
    # Note: the updated new_parametrise_the_community now computes the effective herbivore abundance,
    # incorporating predator contributions from the model derivation.
    S, R, H_i0, m_i, g_i, beta, G, M_modified, A_star, a_matrix, A, epsilon, m_alpha, x_final, raw_g, h =
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
            alpha = alpha,
            hollingII = hollingII,
            h = h
        )
    
    return (
        S, R,
        species_names, herbivore_list, predator_list,
        H_i0, m_i,
        g_i, x_final, beta,   # adjusted effective growth rates, scaling parameter, niche parameters
        G,
        M_modified, A_star, a_matrix, A,
        epsilon, m_alpha,
        raw_g, h
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