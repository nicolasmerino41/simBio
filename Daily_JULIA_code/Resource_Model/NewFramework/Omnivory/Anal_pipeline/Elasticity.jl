############################################################################
############################################################################
########################## ELASTICITIES ####################################
############################################################################
to_float(x) = ForwardDiff.value(x)

function precompute_static_data(cell; species_names=nothing)
    # Assume cell is available and species_names can be computed.
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    S, R = identify_n_of_herbs_and_preds(species_names)
    # println("cell is $cell and species_names are $species_names")
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

    # If species names are not provided, extract them from the cell.
    if isnothing(species_names)
        species_names = extract_species_names_from_a_cell(cell)
    end
    
    # Define the species lists.
    herbivore_list = [sp for sp in species_names if sp in herbivore_names || sp in omnivore_names]
    predator_list  = [sp for sp in species_names if sp in carnivore_names]
    
    # Number of herbivores/omnivores and predators.
    S = length(herbivore_list)
    R = length(predator_list)
    
    # --- Precompute Interaction Matrices ---
    
    # 1. P_matrix (S x R): predation interactions where a predator consumes a herbivore.
    P_matrix = zeros(Float64, S, R)
    for i in 1:S
        # Look up the global index for herbivore species.
        global_h_idx = species_dict[herbivore_list[i]]
        for j in 1:R
            # Look up the global index for predator species.
            global_p_idx = species_dict[predator_list[j]]
            if iberian_interact_NA[global_p_idx, global_h_idx] == 1
                P_matrix[i, j] = 1.0
            end
        end
    end
    
    # 2. O_matrix (S x S): omnivory interactions.
    #    Only rows corresponding to omnivore species have nonzero entries.
    O_matrix = zeros(Float64, S, S)
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
    
    # 3. T_matrix (S x S): herbivore-herbivore (or omnivore-herbivore) consumption interactions.
    T_matrix = zeros(Float64, S, S)
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
    B_matrix = zeros(Float64, R, R)
    for α in 1:R
        global_consumer_idx = species_dict[predator_list[α]]
        for β in 1:R
            global_prey_idx = species_dict[predator_list[β]]
            if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                B_matrix[α, β] = 1.0
            end
        end
    end
    
    # 5. D_matrix (R x R): predator-predator detrimental interactions.
    D_matrix = zeros(Float64, R, R)
    for α in 1:R
        global_prey_idx = species_dict[predator_list[α]]
        for β in 1:R
            global_consumer_idx = species_dict[predator_list[β]]
            if iberian_interact_NA[global_consumer_idx, global_prey_idx] == 1
                D_matrix[α, β] = 1.0
            end
        end
    end

    # Compute indices, species lists, and interaction matrices,
    # all using non-differentiable operations.
    static_data = (
        species_names = species_names,
        herbivore_list = [sp for sp in species_names if sp in herbivore_names || sp in omnivore_names],
        predator_list = [sp for sp in species_names if sp in carnivore_names],
        iberian_interact_NA = iberian_interact_NA,
        H_star = cell_abundance_herbs,
        P_eq = cell_abundance_preds,
        P_matrix = P_matrix,
        O_matrix = O_matrix,
        T_matrix = T_matrix,
        B_matrix = B_matrix,
        D_matrix = D_matrix,
        S = S,
        R = R
    )
    return static_data
end

function differentiable_equilibrium(mu, eps, mean_m_alpha, static_data)
    # Extract precomputed data for herbivores.
    H_star = static_data.H_star         # Vector for herbivores/omnivores
    P_eq   = static_data.P_eq           # Vector for predators
    P_matrix = static_data.P_matrix
    B_matrix = static_data.B_matrix
    D_matrix = static_data.D_matrix

    # Compute candidate nu values from herbivore-predator interactions:
    nu_candidates = map(1:size(P_matrix,2)) do alpha
        H_alpha_tot = sum(P_matrix[:, alpha] .* H_star)
        B_alpha_tot = sum(B_matrix[alpha, :] .* P_eq)
        D_alpha_tot = sum(D_matrix[alpha, :] .* P_eq)
        denominator = eps * (H_alpha_tot + B_alpha_tot) - D_alpha_tot
        denominator > 0 ? (P_eq[alpha] + mean_m_alpha) / denominator : 0.0
    end

    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + 0.05)  # applying a safety factor

    # Compute herbivore equilibrium (K_i) as before.
    T_type = promote_type(eltype(H_star), typeof(mu))
    total_H = sum(H_star)
    S = length(H_star)
    K_i = similar(H_star, T_type)
    for i in 1:S
        P_i_tot = sum(P_matrix[i, :] .* P_eq)
        O_i_tot = sum(static_data.O_matrix[i, :] .* H_star)
        T_i_tot = sum(static_data.T_matrix[i, :] .* H_star)
        K_i[i] = (1 - mu)*H_star[i] + mu*total_H + nu*P_i_tot - eps*nu*O_i_tot + nu*T_i_tot
    end

    # Now compute predator equilibrium (K_alpha)
    R = static_data.R
    K_alpha = zeros(T_type, R)
    for α in 1:R
        H_alpha_tot = sum(P_matrix[:, α] .* H_star)
        B_alpha_tot = sum(B_matrix[α, :] .* P_eq)
        D_alpha_tot = sum(D_matrix[α, :] .* P_eq)
        # Here, we assume the predator equilibrium is given by:
        K_alpha[α] = eps * nu * H_alpha_tot - P_eq[α] + eps * nu * B_alpha_tot - nu * D_alpha_tot
    end

    # Concatenate herbivore and predator equilibria.
    u_star = vcat(K_i, K_alpha)
    return u_star
end

function com_log_eq(p, cell)
    # p contains [mu, eps, mean_m_alpha] as Dual numbers.
    # First, precompute the static data outside of the differentiation:
    local_i, local_j = idx[cell][1], idx[cell][2]  # or whatever cell identifier you use
    static_data = precompute_static_data(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Now call the differentiable core function.
    u_star = differentiable_equilibrium(p[1], p[2], p[3], static_data)
    return NaNMath.log.(u_star)
end

p0 = [0.5, 1.0, 0.1]
J_elasticity = ForwardDiff.jacobian(com_log_eq, p0)

############################################################################
############################################################################
########################## PLOTTING ELASTICITIES ###########################
############################################################################
# --- Prepare species type information ---
# Assume sp_nm corresponds to the rows in J_elasticity.
# For each species, determine its type.
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[idx[1][1], idx[1][2]])
species_types = map(sp -> sp in omnivore_names ? "omnivore" :
                           sp in herbivore_names ? "herbivore" :
                           sp in carnivore_names ? "predator" : "unknown", sp_nm)

# Define colors for each type.
type_colors = Dict("herbivore" => :blue, "omnivore" => :green, "predator" => :red, "unknown" => :gray)
param_names = ["μ", "ε", "m_α"]
# --- Create bar plots for each parameter ---
for param in 1:3
    # Extract sensitivities for the current parameter.
    sens = J_elasticity[:, param]
    
    # Sort species by sensitivity in descending order.
    sorted_idx = sortperm(sens, rev = true)
    sorted_sens = sens[sorted_idx]
    sorted_species = sp_nm[sorted_idx]
    sorted_types = species_types[sorted_idx]
    sorted_colors = [type_colors[t] for t in sorted_types]
    
    # Create a figure and axis.
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1],
        title = "Sensitivity to Parameter $(param_names[param])",
        ylabel = "Sensitivity",
        xlabel = "Species (sorted)",
        xticklabelrotation = π/4,
        xticklabelalign = (:right, :center),
    )
    
    # Plot the bars.
    barplot!(ax, 1:length(sorted_sens), sorted_sens, color = sorted_colors)
    
    # Optionally, add species names as x-tick labels.
    # For many species these may overlap; adjust as needed.
    ax.xticks = (1:length(sorted_species), sorted_species)
    
    display(fig)
end

