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

function o_new_parametrise_the_community(
    species_names::Vector{String};
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    cell_abundance::Vector{Float64} = Float64[],
    o_pred_avg::Float64 = 1.0,
    A_loss_avg::Float64 = 1.0
)
    # Identify herbivores and predators.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)
    
    if S > 0
        if length(cell_abundance) == S
            H_i0 = copy(cell_abundance)
        else
            println("WARNING: cell_abundance length != S, defaulting to uniform abundances.")
            H_i0 = ones(S)
        end
        m_i = [abs(rand(Normal(M_mean, 0.0))) for _ in 1:S]
    else
        H_i0 = Float64[]
        m_i = Float64[]
    end
    
    if S > 0
        V, mu_matrix = generate_competition_matrix(S, mu, symmetrical_competition; check_condition=true)
    else
        V = zeros(0,0)
        mu_matrix = zeros(0,0)
    end
    
    # Construct predator-related matrices.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
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
    
    # Compute d_i.
    d_i = [ m_i[i] / H_i0[i] for i in 1:S ]
    
    # Compute effective herbivore-herbivore interaction matrices.
    O_loss = zeros(S, S)
    O_gain = zeros(S, S)
    for i in 1:S
        for j in 1:S
            if i != j
                O_loss[i,j] = o_pred_avg * (mu_matrix[i,j] / d_i[i])
                O_gain[i,j] = o_pred_avg * epsilon_val * (mu_matrix[j,i] / d_i[j])
            end
        end
    end
    
    # Set herbivore intrinsic growth rates to 1.
    g_i = ones(S)
    
    # Compute predator-induced loss matrix.
    A_star = zeros(S, R)
    for i in 1:S
        for alph in 1:R
            A_star[i, alph] = A_loss_avg * (mu_predation / d_i[i])
        end
    end
    
    # Compute effective observed herbivore abundance.
    H_i0_eff = zeros(S)
    for i in 1:S
        comp_sum = 0.0
        for j in 1:S
            if i != j
                comp_sum += mu_matrix[i,j] * H_i0[j]
            end
        end
        net_O = 0.0
        for j in 1:S
            if i != j
                net_O += (O_gain[i,j] - O_loss[i,j]) * H_i0[j]
            end
        end
        H_i0_eff[i] = H_i0[i] + comp_sum + net_O
    end
    
    beta = zeros(S)
    
    return (
        S = S, R = R,
        H_i0_eff = H_i0_eff, m_i = m_i,
        g_i = g_i, beta = beta,
        O_loss = O_loss, O_gain = O_gain,
        A_star = A_star, a_matrix = a_matrix, A = A,
        epsilon = epsilon, m_alpha = m_alpha,
        x = 0.0, raw_g = zeros(S), h = 0.0,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = species_names
    )
end

function o_new_setup_community_from_cell(
    i::Int, j::Int;
    M_mean::Float64 = 0.1,            # average herbivore mortality scale
    mu::Float64 = 0.5,                # competition parameter among herbivores
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,      # predator mortality scale
    epsilon_val::Float64 = 1.0,       # conversion efficiency for predation gains
    mu_predation::Float64 = 0.01,     # average strength for predator-herbivore interactions (for both A and A*)
    iberian_interact_NA::NamedMatrix{Float64} = iberian_interact_NA,
    species_dict::Dict{String,Int} = species_dict,
    # cell_abundance::Vector{Float64} = Float64[],
    species_names::Vector{String} = nothing,
    o_pred_avg::Float64 = 1.0         # average strength for herbivore-herbivore interactions
)
    # 1) Retrieve the cell and extract species present.
    cell = DA_birmmals_with_pi_corrected[i, j]
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    
    # 2) Build herbivore abundance vector.
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
    # 'artificial_pi' is assumed to be defined globally.
    if artificial_pi
        cell_abundance_herbs = ones(length(cell_abundance_herbs))
    end
    println("We at least got here")
    # 3) Call the new parametrisation function for the new framework.
    return o_new_parametrise_the_community(
         species_names;
         M_mean = M_mean,
         mu = mu,
         symmetrical_competition = symmetrical_competition,
         mean_m_alpha = mean_m_alpha,
         epsilon_val = epsilon_val,
         mu_predation = mu_predation,
         iberian_interact_NA = iberian_interact_NA,
         species_dict = species_dict,
         cell_abundance = cell_abundance_herbs,
         o_pred_avg = o_pred_avg,
         A_loss_avg = 1.0,     # Use 1.0 as default for predator-induced loss strength
    )
end

o_new_setup_community_from_cell(
    
)