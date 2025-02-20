function new_dynamics!(du, u, p, t)
    # Here, we assume the parameters are packed as:
    # S, R, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon, m_alpha
    S, R, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon, m_alpha = p

    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics (new formulation)
    # For each herbivore species, the per capita growth is:
    #   (g_i + G_i) * (beta_i/(1+beta_i)) * (1 - (H_i + sum_j(M_modified[i,j]*H_j)) / H_i^(obs)(mu))
    for i in 1:S
        if H[i] > 0.0
            duH[i] = H[i] * (g_i[i] + G[i]) * (beta[i] / (1 + beta[i])) *
                     (1 - (H[i] + sum(M_modified[i, x] * H[x] for x in 1:S)) / H_i0[i])
        else
            duH[i] = 0.0
        end
    end

    # Predator dynamics (unchanged)
    for α in 1:R
        if P[α] > 0.0
            predation_sum = 0.0
            for j in 1:S
                predation_sum += a_matrix[j, α] * H[j]
            end
            predator_interactions = 0.0
            for β in 1:R
                predator_interactions += A[α, β] * P[β]
            end
            duP[α] = P[α] * (epsilon[α] * predation_sum - m_alpha[α] + predator_interactions)
        else
            duP[α] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end
