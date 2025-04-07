function g_abstract_parametrise_the_community(
    S::Int, O::Int, R::Int;
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    conn = 0.2,                # target connectance (probability an interaction exists)
    cell_abundance_h::Vector{Float64} = Float64[],  # default: not provided → ones
    cell_abundance_p::Vector{Float64} = Float64[],  # default: not provided → ones
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 1.0,
    randomise_attack_rates::Bool = false,
    attack_rate_sd::Float64 = 0.2
)
    # Convert parameters to Float64.
    mu             = Float64(mu)
    epsilon_val  = Float64(epsilon_val)
    mean_m_alpha = Float64(mean_m_alpha)
    delta_nu     = Float64(delta_nu)
    d_alpha      = Float64(d_alpha)
    d_i          = Float64(d_i)
    r_omni_proportion = Float64(r_omni_proportion)
    conn         = Float64(conn)

    # Total number in the herbivore compartment.
    S_total = S + O

    # Define synthetic species names.
    herbivore_list = vcat(["Herbivore $i" for i in 1:S],
                          ["Omnivore $i" for i in 1:O])
    predator_list  = ["Predator $i" for i in 1:R]
    sp_nm = vcat(herbivore_list, predator_list)

    # Process equilibrium abundances:
    H_eq = (S_total > 0) ? (length(cell_abundance_h) == S_total ? copy(cell_abundance_h) : ones(S_total)) : Float64[]
    P_eq = (R > 0) ? (length(cell_abundance_p) == R ? copy(cell_abundance_p) : ones(R)) : Float64[]
    H_star = copy(H_eq)

    # Intrinsic growth rates for herbivores:
    r_vec = vcat(ones(S), ones(O) * r_omni_proportion)  # pure herbivores grow at rate 1, omnivores at r_omni_proportion

    # For predators, assign mortality and conversion efficiency:
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Conversion efficiency scaled by (d_i/d_alpha)
        epsilon_vec = fill(epsilon_val, R) * (d_i / d_alpha)
    else
        m_alpha = Float64[]
        epsilon_vec = Float64[]
    end

    # --- Build the merged interaction matrix A ---
    # A is of size n_total x n_total, where n_total = S_total + R.
    n_total = S_total + R
    A = zeros(Float64, n_total, n_total)
    # We loop over all pairs (i,j). For each potential interaction, we decide (with probability conn)
    # whether an interaction exists. If yes, then we assign:
    #   - The effect on the target (row i) is +1.0 (representing a loss in biomass when attacked).
    #   - The effect on the consumer (row j) is +epsilon_val (a benefit from consumption).
    #
    # Here we distinguish by the types of the species:
    #   - If the target is herbivore (i ≤ S_total) and the consumer is predator (j > S_total), or
    #   - If the target is herbivore and the consumer is omnivore (j ≤ S_total and j > S),
    #   - If both are predators, or if a herbivore (even a pure herbivore) consumes a predator,
    # we use the same rule.
    for i in 1:n_total
        for j in 1:n_total
            if rand() < conn  # an interaction exists
                # Determine species types.
                target_is_herb = (i ≤ S_total)
                consumer_is_pred = (j > S_total)
                consumer_is_omnivore = (j ≤ S_total) && (j > S)  # indices S+1..S_total are omnivores
                consumer_is_herb = (j ≤ S_total) && (j ≤ S)
                # We assign values only if the entry is not already set (to avoid double counting).
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

    # Scale each row of A:
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
    # Let X = [H_star; P_eq]
    X = vcat(H_star, P_eq)
    nu_candidates = Float64[]
    # For each predator (indices S_total+alpha), compute:
    for alpha in 1:R
        idx_pred = S_total + alpha
        gain = sum(A[idx_pred, :] .* X)
        if gain > 0
            push!(nu_candidates, (P_eq[alpha] + m_alpha[alpha]) / gain)
        end
    end
    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + delta_nu)  # apply safety margin

    # --- Embed ν into A (i.e. multiply all entries by ν) ---
    A .= nu .* A

    # --- Optionally randomise attack rates ---
    if randomise_attack_rates
        # Use a truncated Normal distribution with mean 1 and sd attack_rate_sd.
        multiplier_dist = Truncated(Normal(1.0, attack_rate_sd), 0, Inf)
        for i in axes(A,1), j in axes(A,2)
            A[i,j] *= rand(multiplier_dist)
        end
    end

    # --- Recalculate carrying capacities ---
    # For herbivores:
    K_i = zeros(Float64, S_total)
    total_H = sum(H_star)
    for i in 1:S_total
        # Compute trophic effect from all interactions.
        trophic_effect = sum(A[i, :] .* X)
        comp_effect = sum(C[i, :] .* H_star)
        K_i[i] = (1 - mu) * H_star[i] + trophic_effect + comp_effect
    end
    # For predators:
    K_alpha = zeros(Float64, R)
    for alpha in 1:R
        idx_pred = S_total + alpha
        trophic_effect_pred = sum(A[idx_pred, :] .* X)
        K_alpha[alpha] = trophic_effect_pred - P_eq[alpha]
    end

    # (Optional) Compute derived ratios (not used here but could be returned)
    new_di = r_vec ./ K_i
    new_da = m_alpha ./ K_alpha

    return (
        S = S_total, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_vec, K_i = K_i,
        mu = mu, nu = nu,
        A_matrix = A,
        C_matrix = C,
        epsilon = epsilon_vec, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = sp_nm,
        H_star = H_star, P_star = P_eq
    )
end
