function omnivore_dynamics!(du, u, p, t)
    # p is a tuple containing:
    # S, R, K_i, r_i, mu, nu, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix
    S, R, K_i, r_i, mu, nu, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_omni = p

    # Partition state vector: first S entries are herbivore (and omnivore) densities, next R are predator densities.
    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(eltype(u), S)
    duP = zeros(eltype(u), R)

    # Total herbivore biomass.
    H_tot = sum(H)

    # Herbivore/Omnivore dynamics:
    # dH_i/dt = H_i * r_i * { 1 - [ (1-μ)H_i + μ H_tot + ν P_i^tot - εν O_i^tot + ν T_i^tot ] / K_i }
    # where:
    #   P_i^tot = ∑₍α₎ P_matrix[i,α]*P[α]
    #   O_i^tot = ∑₍j₎ O_matrix[i,j]*H[j]   (biomass of prey available to omnivore i; zero if species i is not omnivore)
    #   T_i^tot = ∑₍j₎ T_matrix[i,j]*H[j]   (biomass of herbivores consuming species i)
    # Note: The conversion efficiency for omnivory (ε) is assumed uniform; here we use ε = epsilon[1] if available.
    eps_herb = length(epsilon) > 0 ? epsilon[1] : 1.0

    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            # Compute total predation pressure on herbivore i.
            P_i_tot = 0.0
            for alpha in 1:R
                P_i_tot += P_matrix[i, alpha] * P[alpha]
            end

            # Compute total biomass of prey for omnivory (only nonzero for omnivores).
            O_i_tot = 0.0
            for j in 1:S
                O_i_tot += O_matrix[i, j] * H[j]
            end

            # Compute total biomass of herbivores that consume species i.
            T_i_tot = 0.0
            for j in 1:S
                T_i_tot += T_matrix[i, j] * H[j]
            end

            effective_density = (1 - mu) * H[i] + mu * H_tot + nu * P_i_tot - eps_herb * nu_omni * O_i_tot + nu_omni * T_i_tot
            duH[i] = H[i] * r_i[i] * (1 - effective_density / K_i[i])
        end
    end

    # Predator dynamics:
    # dP_α/dt = P_α * m_α * { (εν H_α^tot - P_α + εν B_α^tot - ν D_α^tot) / K_α - 1 }
    # where:
    #   H_α^tot = ∑₍i₎ P_matrix[i,α]*H[i]
    #   B_α^tot = ∑₍β₎ B_matrix[α,β]*P[β]
    #   D_α^tot = ∑₍β₎ D_matrix[α,β]*P[β]
    for alpha in 1:R
        if P[alpha] > EXTINCTION_THRESHOLD
            H_alpha_tot = 0.0
            for i in 1:S
                H_alpha_tot += P_matrix[i, alpha] * H[i]
            end

            B_alpha_tot = 0.0
            for beta in 1:R
                B_alpha_tot += B_matrix[alpha, beta] * P[beta]
            end

            D_alpha_tot = 0.0
            for beta in 1:R
                D_alpha_tot += D_matrix[alpha, beta] * P[beta]
            end

            duP[alpha] = P[alpha] * m_alpha[alpha] * ((epsilon[alpha] * nu * H_alpha_tot - P[alpha] + epsilon[alpha] * nu * B_alpha_tot - nu * D_alpha_tot) / K_alpha[alpha] - 1)
        end
    end

    # Write the computed derivatives back into du.
    for i in 1:S
        du[i] = duH[i]
    end
    for alpha in 1:R
        du[S+alpha] = duP[alpha]
    end
end
