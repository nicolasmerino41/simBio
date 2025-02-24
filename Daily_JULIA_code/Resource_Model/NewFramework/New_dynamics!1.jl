function new_dynamics!(du, u, p, t)
    # Unpack parameters for the new model:
    # S       : number of herbivore species
    # R       : number of predator species
    # H_i0    : effective observed herbivore abundances (H_i^(obs)(μ,C,G); vector of length S)
    # m_i     : natural mortality rates for herbivores (vector, length S)
    # g_i     : effective (observed) herbivore growth rates (vector, length S)
    # beta    : niche parameters for herbivores (vector, length S), where beta = g_i/m_i - 1
    # M_mod   : competition matrix among herbivores (S×S), representing μ_{ij} (with zeros on the diagonal)
    # A_star  : effective predation matrix for herbivores (S×R), representing nondimensionalized attack rates A*_{iα}
    # A_pred  : predator attack rate matrix (R×S), where A_pred[α,i] is the attack rate of predator α on herbivore i
    # P0      : predator overpopulation thresholds (vector of length R, i.e. P_{α0} = m_α/d_α)
    # B       : predator interaction matrix (R×R), representing nondimensional interactions B_{αβ}
    # m_alpha : predator mortality rates (vector of length R)
    S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, B, m_alpha = p

    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics:
    # New equation:
    #   dH_i/dt = g_i * (beta_i/(1+beta_i)) *
    #             [ 1 - ( H_i + Σ_j μ_{ij}H_j + Σ_α A*_{iα} P_α ) / H_i0 ] * H_i
    for i in 1:S
        if H[i] > 0.0
            # Compute the total density affecting species i:
            competition_term = 0.0
            for j in 1:S
                competition_term += M_mod[i, j] * H[j]
            end
            predation_term = 0.0
            for α in 1:R
                predation_term += A_star[i, α] * P[α]
            end
            effective_density = H[i] + competition_term + predation_term
            duH[i] = H[i] * g_i[i] * (beta[i] / (1 + beta[i])) *
                     (1 - effective_density / H_i0[i])
        else
            duH[i] = 0.0
        end
    end

    # Predator dynamics:
    # New equation:
    # dP_α/dt = m_alpha * [ (Σ_i A_pred[α,i]*H_i - P0[α] - Σ_β B[α,β]*P_β) / P0[α] ] * P_α
    for α in 1:R
        if P[α] > 0.0
            attack_sum = 0.0
            for i in 1:S
                attack_sum += A_pred[α, i] * H[i]
            end
            interaction_sum = 0.0
            for β in 1:R
                interaction_sum += B[α, β] * P[β]
            end
            duP[α] = P[α] * m_alpha[α] * ((attack_sum - P0[α] - interaction_sum) / P0[α])
        else
            duP[α] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end
