function new_dynamics!(du, u, p, t)
    # Unpack the first 12 parameters:
    # S       : number of herbivore species
    # R       : number of predator species
    # H_i0    : effective observed herbivore abundances (vector of length S)
    # m_i     : herbivore mortality rates (vector, length S)
    # g_i     : effective herbivore growth rates (vector, length S)
    # beta    : niche parameters for herbivores (vector, length S)
    # M_mod   : modified competition matrix among herbivores (S×S)
    # A_star  : nondimensional attack rates (S×R)
    # A_pred  : predator attack rate matrix (R×S)
    # P0      : predator overpopulation thresholds (vector of length R)
    # B       : predator interaction matrix (R×R)
    # m_alpha : predator mortality rates (vector of length R)
    S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha = p[1:12]
    
    # Check if we have an extra parameter h (handling time) for Holling type II:
    holling = false
    if length(p) >= 13
        h = p[13]
        holling = true
    end

    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics:
    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            comp_term = 0.0
            for j in 1:S
                comp_term += M_mod[i, j] * H[j]
            end
            pred_term = 0.0
            for alph in 1:R
                pred_term += A_star[i, alph] * P[alph]
            end
            effective_density = H[i] + comp_term + pred_term
            duH[i] = H[i] * g_i[i] * (beta[i] / (1 + beta[i])) *
                     (1 - effective_density / H_i0[i])
        elseif H[i] < EXTINCTION_THRESHOLD
            duH[i] = H[i] < 0.0 ? abs(H[i]) : -H[i]
        elseif iszero(H[i])
            duH[i] = 0.0
        end
    end

    # Predator dynamics:
    for alph in 1:R
        if P[alph] > EXTINCTION_THRESHOLD  
            attack_sum = 0.0
            for i in 1:S
                attack_sum += A_pred[alph, i] * H[i]
            end
            # Apply Holling type II if requested, otherwise use linear feeding:
            sat_attack = holling ? (attack_sum / (1 + h * attack_sum)) : attack_sum
            
            interact_sum = 0.0
            for bet in 1:R
                interact_sum += B[alph, bet] * P[bet]
            end
            duP[alph] = P[alph] * m_alpha[alph] * ((sat_attack - P0[alph] - interact_sum) / P0[alph])
        
        elseif P[alph] < EXTINCTION_THRESHOLD
            duP[alph] = P[alph] < 0.0 ? abs(P[alph]) : -P[alph]
        elseif iszero(P[alph])
            duP[alph] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end
