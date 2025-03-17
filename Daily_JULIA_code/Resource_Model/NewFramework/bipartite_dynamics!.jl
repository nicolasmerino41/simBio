function bipartite_dynamics!(du, u, p, t)
    # p is a tuple containing:
    # S         : number of herbivores
    # R         : number of predators
    # K_i       : vector of herbivore carrying capacities (K_i)
    # r_i       : herbivore intrinsic growth rates
    # mu        : competition parameter for herbivores
    # nu        : effective predation rate (loss rate for herbivores)
    # P_matrix  : S x R predation incidence matrix (1 if predator alpha preys on herbivore i)
    # epsilon   : vector of predator conversion efficiencies
    # m_alpha   : vector of predator mortality rates
    # K_alpha   : vector of predator carrying capacities (K_alpha)
    S, R, K_i, r_i, mu, nu, P_matrix, epsilon, m_alpha, K_alpha = p
    
    # Partition state vector: first S entries are herbivore densities, next R are predator densities.
    H = @view u[1:S]
    P = @view u[S+1:S+R]
    
    duH = zeros(S)
    duP = zeros(R)
    
    # Total herbivore biomass.
    H_tot = sum(H)
    
    # Herbivore dynamics:
    # dH_i/dt = H_i * r_i * [ 1 - ((1-mu)*H_i + mu*H_tot + nu * P_i_tot) / K_i ]
    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            # P_i_tot: total predation pressure on herbivore i.
            P_i_tot = 0.0
            for alpha in 1:R
                P_i_tot += P_matrix[i, alpha] * P[alpha]
            end
            effective_density = (1 - mu) * H[i] + mu * H_tot + nu * P_i_tot
            duH[i] = H[i] * r_i[i] * (1 - effective_density / K_i[i])
        # elseif H[i] < EXTINCTION_THRESHOLD
        #     duH[i] = H[i] < 0.0 ? abs(H[i]) : -H[i]
        # elseif iszero(H[i])
        #     duH[i] = 0.0
        end
    end
    
    # Predator dynamics:
    # dP_alpha/dt = P_alpha * m_alpha * [ (epsilon*nu*H_alpha_tot - P_alpha) / K_alpha - 1 ]
    for alpha in 1:R
        if P[alpha] > EXTINCTION_THRESHOLD
            # H_alpha_tot: total available prey for predator alpha.
            H_alpha_tot = 0.0
            for i in 1:S
                H_alpha_tot += P_matrix[i, alpha] * H[i]
            end
            duP[alpha] = P[alpha] * m_alpha[alpha] * ( (epsilon[alpha] * nu * H_alpha_tot - P[alpha]) / K_alpha[alpha] - 1 )
        # elseif P[alpha] < EXTINCTION_THRESHOLD
        #     duP[alpha] = P[alpha] < 0.0 ? abs(P[alpha]) : -P[alpha]
        # elseif iszero(P[alpha])
        #     duP[alpha] = 0.0
        end
    end
    
    # Write derivatives back into du.
    for i in 1:S
        du[i] = duH[i]
    end
    for alpha in 1:R
        du[S+alpha] = duP[alpha]
    end
end
