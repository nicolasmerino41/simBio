##############################################################################################
##############################################################################################
############################ SAME CODE BUT INSIDE A FUNCTION #################################

# Assuming the following globals:
# - birmmals_biomass_fixed::DataFrame with columns:
#   species (String), mean_density (Float64), bodyMass (Float64), biomass (Float64)
# - herbivore_names::Vector{String} containing names of herbivores
# - predator_names::Vector{String} containing names of predators
# - DA_birmmals::Matrix{MyBirmmals} containing the DA cells
# - MyBirmmals is a struct with a field 'a' that indicates presence/absence (or abundance)
#   of each species (length matches rows in birmmals_biomass_fixed)

# Example:
# struct MyBirmmals
#     a::SVector{207,Float64}
# end

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
    # Determine the predators present in this cell
    cell_predators = [sp for sp in species_names if sp in predator_names]
    
    # Initialize the output dictionary
    predator_has_prey = Dict{String,Bool}()
    
    # Iterate over each predator in the cell
    for pred in cell_predators
        # Get the global index of the predator from species_dict
        global_pred_idx = species_dict[pred]
        found_prey = false

        # Check among all species in the cell whether there is at least one prey
        # (Here we consider both herbivores and predators as potential prey.)
        for candidate in species_names
            # Skip checking self-predation if desired
            if candidate == pred
                continue
            end

            # Get the global index for the candidate prey
            global_candidate_idx = species_dict[candidate]
            # Check if the interaction matrix indicates that `pred` eats `candidate`
            if iberian_interact_NA[global_pred_idx, global_candidate_idx] == 1
                found_prey = true
                break  # Stop checking once we find a prey
            end
        end

        predator_has_prey[pred] = found_prey
    end

    # Check if all predators have at least one prey
    allHavePrey = all(values(predator_has_prey))
    
    if allHavePrey
        # All predators have prey: Return true along with 0 and an empty vector
        return (true, 0, String[])
    else
        # Identify predators that do not have any prey
        predators_without_prey = [pred for (pred, hasprey) in predator_has_prey if !hasprey]
        num_without_prey = length(predators_without_prey)
        return (false, num_without_prey, predators_without_prey)
    end
end

###############################################################################
# 1) PARAMETRISE THE COMMUNITY
###############################################################################
function parametrise_the_community(
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

    # Standard deviations for random draws:
    m_standard_deviation::Float64 = 0.0,
    h_standard_deviation::Float64 = 0.0,

    # Cell-level herbivore abundance (length S) to build hat_H
    cell_abundance::Vector{Float64} = Float64[],

    real_H0::Bool = false,
    H0_vector::Vector{Float64} = nothing
)

    ############################################################################
    # (A) Identify herbivores vs. predators from species_names
    ############################################################################
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)

    ############################################################################
    # (B) Build H_i^0 and m_i from random draws, if herbivores exist
    ############################################################################
    if S > 0
        if !real_H0
            # fallback approach: random draws
            H0_mean = NPP / S
            H_i0 = [abs(rand(Normal(H0_mean, h_standard_deviation))) for _ in 1:S]
            m_i  = [abs(rand(Normal(M_mean,   m_standard_deviation))) for _ in 1:S]
        else
            # real_H0 == true => must get from H0_vector
            if isnothing(H0_vector)
                error("H0_vector must be provided if real_H0 is true")
            end
            # We'll build a temporary array of length S to store the herbivores' H0
            H_i0_temp = Float64[]  # will hold the subset from H0_vector
    
            # For each herbivore 'sp' in herbivore_list, find its index in species_dict_herbivores_in_birmmals
            for sp in herbivore_list
                local_idx = species_dict_herbivores_in_birmmals[sp]
                push!(H_i0_temp, H0_vector[local_idx])
            end
    
            # Then scale by baseline factor (NPP / S)
            H_i0_baseline = NPP / S
            # e.g. multiply each value in H_i0_temp by H_i0_baseline
            H_i0 = H_i0_baseline .* H_i0_temp
    
            # for mortality, just do random draws or do something else
            m_i  = [abs(rand(Normal(M_mean, m_standard_deviation))) for _ in 1:S]
        end
    else
        # No herbivores => empty arrays
        H_i0 = Float64[]
        m_i  = Float64[]
    end

    ############################################################################
    # (C) Baseline competition among herbivores => mu_matrix
    ############################################################################
    if S > 0
        V, mu_matrix = generate_competition_matrix(S, mu, symmetrical_competition; check_condition=true)
    else
        V          = zeros(0,0)
        mu_matrix  = zeros(0,0)
    end

    ############################################################################
    # PredatorVSpredator and herbivoreVspred feeding matrices =(a_matrix & A)
    ############################################################################
    a_matrix = zeros(S, R)
    A        = Matrix{Float64}(I, R, R)

    # negative diag for predator self-regulation
    for row_pred in 1:R
        A[row_pred, row_pred] = -1.0
    end

    # fill a_matrix if pred_j eats herb_i
    for herb_i in 1:S
        global_herb_idx = species_dict[herbivore_list[herb_i]]
        for pred_j in 1:R
            global_pred_idx = species_dict[predator_list[pred_j]]
            if iberian_interact_NA[global_pred_idx, global_herb_idx] == 1
                a_matrix[herb_i, pred_j] = mu_predation
            end
        end
    end

    # fill A if pred_j eats pred_k
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

    # predator mortality, assimilation
    m_alpha = fill(mean_m_alpha, R)
    epsilon = fill(epsilon_val, R)

    ############################################################################
    # (E) Compute C_{ij}, G_i => M_modified
    ############################################################################
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
                    val += epsilon[alpha]*a_matrix[r, alpha]*A_inv[alpha, beta]*a_matrix[c, beta]
                end
            end
            C[r,c] = val
        end
    end

    for r in 1:S
        val = 0.0
        for alpha in 1:R
            for beta in 1:R
                val += a_matrix[r, alpha]*A_inv[alpha, beta]*m_alpha[beta]
            end
        end
        G[r] = val
    end

    M_modified = copy(mu_matrix)
    for r in 1:S, c in 1:S
        M_modified[r,c] += C[r,c]*(H_i0[r]/m_i[r])
    end

    ############################################################################
    # (F) Incorporate the "cell_abundance" => define hat_H
    ############################################################################
    # According to doc:
    #   hat_p[i] = [hat_H[i] + ∑ (M_modified[i,j]*hat_H[j])] / h_i[i]
    # then p_i = hat_p[i] / ∑ hat_p.

    # We'll define barH from H_i0
    barH = (S>0) ? mean(H_i0) : 1.0

    # define h_i from doc
    local_h = zeros(S)
    for i in 1:S
        local_h[i] = H_i0[i]/(S*barH)
    end

    # define a localHatH vector from cell_abundance if length(cell_abundance) == S
    # otherwise fallback
    localHatH = zeros(S)
    if (length(cell_abundance) == S)
        for i in 1:S
            localHatH[i] = cell_abundance[i]
        end
        # println("good")
    else
        # fallback => uniform
        for i in 1:S
            localHatH[i] = 1.0
        end
        println("ERROR cell abundance != S")
    end

    # compute hat_p
    hat_p = zeros(S)
    for i in 1:S
        numerator = localHatH[i]
        for j in 1:S
            numerator += M_modified[i,j]*localHatH[j]
        end
        hat_p[i] = numerator / local_h[i]
    end

    # normalise
    sum_hat_p = sum(hat_p)
    p_vec = zeros(S)
    if sum_hat_p>0
        p_vec .= hat_p ./ sum_hat_p
    else
        p_vec .= 1.0/S
    end

    ############################################################################
    # (G) Solve the doc's NPP eq => x => define g_i
    ############################################################################
    IplusM = I + M_modified
    IM_inv = inv(IplusM)

    A_vec = zeros(S)
    B_vec = zeros(S)
    for r in 1:S
        A_val = 0.0
        B_val = 0.0
        for c in 1:S
            A_val += IM_inv[r,c]*p_vec[c]
            B_val += IM_inv[r,c]*(G[c]/m_i[c] - 1.0)
        end
        A_vec[r] = A_val
        B_vec[r] = B_val
    end

    A_coef = 0.0
    B_coef = 0.0
    C_coef = 0.0
    for r in 1:S
        A_coef += p_vec[r]*m_i[r]*H_i0[r]*A_vec[r]
        B_coef += p_vec[r]*m_i[r]*H_i0[r]*B_vec[r] + G[r]*H_i0[r]*A_vec[r]
        C_coef += G[r]*H_i0[r]*B_vec[r]
    end
    C_coef -= NPP

    disc = B_coef^2 - 4A_coef*C_coef
    if disc < 0
        error("No real solution for x found!")
    end
    x_candidates = [
        (-B_coef + sqrt(disc)) / (2A_coef),
        (-B_coef - sqrt(disc)) / (2A_coef)
    ]
    pos_x = filter(x->x>0, x_candidates)
    if isempty(pos_x)
        error("No positive x found!")
    end
    x_final = maximum(pos_x)

    # define g_i
    g_i = [x_final * p_vec[i] * m_i[i] for i in 1:S]

    return (
        S, R,
        H_i0, m_i,
        p_vec,
        localHatH,       # final resource allocation proportions
        x_final, g_i,
        G,           # from predator feedback
        M_modified,
        a_matrix, A,
        epsilon, m_alpha
    )
end


################################################################################
# 2) SETUP COMMUNITY FROM CELL
################################################################################

function setup_community_from_cell(
    i::Int, j::Int;
    NPP::Float64 = 1000.0,
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    iberian_interact_NA::NamedMatrix{Float64}=iberian_interact_NA,
    species_dict::Dict{String,Int}=species_dict,
    m_standard_deviation::Float64 = 0.0,
    h_standard_deviation::Float64 = 0.0,
    artificial_pi = false,
    real_H0::Bool = false,
    H0_vector::Vector{Float64} = nothing,
    species_names::Vector{String} = nothing
)

    # 1) Grab the cell => find which species are present
    cell = DA_birmmals_with_pi_corrected[i, j]
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)  # list of species with .a[k] != 0
    end
    S, R = identify_n_of_herbs_and_preds(species_names)

    # 2) Build cell_abundance specifically for the herbivores
    #    so that param. function can interpret it as hat_H
    #    We'll do it in the same order as herbivore_list => only length S
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]

    cell_abundance_herbs = Float64[]
    herbivore_list_cell = String[]  # if you also want the names in parallel

    for sp in herbivore_list
        local_idx = species_dict_herbivores_in_birmmals[sp]  # position in birmmals_names
        val = cell.a[local_idx]                              # abundance from cell

        if val > 0.0
            push!(cell_abundance_herbs, val)
            push!(herbivore_list_cell, sp)
        end
    end
    
    if artificial_pi
        cell_abundance_herbs = ones(length(cell_abundance_herbs))
    end

    # 3) Parametrise the community => pass "cell_abundance_herbs"
    S, R,
    H_i0, m_i,
    p_vec, 
    localHatH,
    x_final, g_i,
    G, M_modified,
    a_matrix, A,
    epsilon, m_alpha =
        parametrise_the_community(
            species_names;
            NPP=NPP,
            M_mean=M_mean,
            mu=mu,
            symmetrical_competition=symmetrical_competition,
            mean_m_alpha=mean_m_alpha,
            epsilon_val=epsilon_val,
            mu_predation=mu_predation,
            iberian_interact_NA=iberian_interact_NA,
            species_dict=species_dict,
            m_standard_deviation=m_standard_deviation,
            h_standard_deviation=h_standard_deviation,
            cell_abundance=cell_abundance_herbs,
            real_H0=real_H0,
            H0_vector=H0_vector
        )

    return (
        S, R,
        species_names, herbivore_list, predator_list,
        H_i0, m_i,
        p_vec, x_final, g_i,  # final resource allocation & growth
        localHatH,
        G,
        M_modified, a_matrix, A,
        epsilon, m_alpha
    )
end

####### A FUNCTION TO EXTRACT INFO OF A SPECIES (I.E. pi_value, whether it's herbivore or predator, etc)
function extract_info_of_a_species(species_name::String)
    # 1) Extract the pi value from the cell
    if species_name in birmmals_biomass_fixed.species
        pi_value = birmmals_biomass_fixed[birmmals_biomass_fixed.species .== species_name, :biomass][1]
        rounded_pi_value = round(pi_value[1], digits=2)
        # Correctly compute the rank of the species' biomass
        sorted_biomass = sort(birmmals_biomass_fixed[:, :biomass], rev=true)  # Sort biomass values in descending order
        rank = findfirst(x -> x == pi_value, sorted_biomass)
    else
        error("Species $species_name not found in birmmals_biomass_fixed")
    end

    # 2) Get the species' positions in `birmmals_names` and `spain_names`
    position_in_birmmals = findfirst(x -> x == species_name, birmmals_names)
    position_in_spain_names = findfirst(x -> x == species_name, spain_names)

    # 3) Check if the species is a herbivore or predator
    if species_name in herbivore_names_as_birmmals
        what = "herbivore"
    elseif species_name in predator_names_as_birmmals
        what = "predator"
    else
        error("Species $species_name is neither herbivore nor predator")
    end
    
    # 4) Check its average H0 value from H0_DA
    if species_name in herbivore_names_as_birmmals
        H0_value_vector = []
        for i in idx
            push!(H0_value_vector, H0_DA[i].a[species_dict_herbivores_in_birmmals[species_name]])
        end
        average_H0_value = mean(H0_value_vector)
        H0_value_sd = std(H0_value_vector)
    end 

    # 5) Print the output
    println("$species_name is a $what with p_i $rounded_pi_value, which falls $rank/$(length(birmmals_biomass_fixed[:, :biomass]))")
    println("Its position in birmmals_names is $position_in_birmmals and spain_names $position_in_spain_names")

    if species_name in herbivore_names_as_birmmals
        println("Its average H0 value is $average_H0_value")
        println("Its standard deviation is $H0_value_sd")
    end
end

# extract_info_of_a_species("Apodemus sylvaticus")
