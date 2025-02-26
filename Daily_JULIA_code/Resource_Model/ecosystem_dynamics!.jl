# ODE definition:
function ecosystem_dynamics!(du, u, p, t)
    S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = p

    ext_threshold = 1.0e-6
    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics
    for i in 1:S
        if H[i] > 0.0
            m_ii = m_i[i]
            numerator = (g_i[i] + G[i])/m_ii - 1.0
            interaction_sum = H[i]
            for x in 1:S
                interaction_sum += M_modified[i,x]*H[x]
            end
            duH[i] = H[i]*m_ii*(numerator - interaction_sum/H_i0[i])
        else
            duH[i] = 0.0
        end
    end

    # Predator dynamics
    for α in 1:R
        if P[α] > 0.0
            predation_sum = 0.0
            for j in 1:S
                predation_sum += a_matrix[j, α]*H[j]
            end
            predator_interactions = 0.0
            for β in 1:R
                predator_interactions += A[α, β]*P[β]
            end
            duP[α] = P[α]*(epsilon[α]*predation_sum - m_alpha[α] + predator_interactions)
        else
            duP[α] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end