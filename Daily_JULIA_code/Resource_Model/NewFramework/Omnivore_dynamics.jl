function omnivore_dynamics!(du, u, p, t)
    # p is a tuple with:
    # S         : number of herbivores
    # R         : number of predators
    # H_i0_eff  : effective (baseline) herbivore abundances
    # g_i       : herbivore intrinsic growth rates
    # mu        : competition parameter for herbivores
    # O_loss    : herbivore-herbivore loss matrix (O_{ij}^*)
    # O_gain    : herbivore-herbivore gain matrix (O_{ji})
    # A_star    : effective predation loss matrix (A_{αi}^*)
    # a_matrix  : raw predator-herbivore interaction matrix (a_{αi})
    # epsilon   : conversion efficiencies for predators
    # B         : predator-predator interaction matrix
    # m_alpha   : predator mortality (and sets P_{α}^{(0)} = m_alpha)
    S, R, H_i0_eff, g_i, mu, O_loss, O_gain, A_star, a_matrix, epsilon, B, m_alpha = p

    # Partition state vector into herbivore (H) and predator (P) densities.
    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics
    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            H_tot = sum(H)
            # Baseline: mean-field mixture of self and total herbivore density.
            baseline = (1 - mu) * H[i] + mu * H_tot

            # Sum over herbivore–herbivore interactions.
            loss_herb = 0.0
            gain_herb = 0.0
            for j in 1:S
                if i != j
                    loss_herb += O_loss[i, j] * H[j]
                    gain_herb += O_gain[i, j] * H[j]
                end
            end

            # Sum predation losses from all predators.
            loss_pred = 0.0
            for α in 1:R
                loss_pred += A_star[i, α] * P[α]
            end

            # Combine the effects into an effective density.
            effective_density = baseline - loss_herb + gain_herb - loss_pred

            # Herbivore dynamics follow a logistic-like form.
            duH[i] = H[i] * g_i[i] * (1 - effective_density / H_i0_eff[i])
        elseif H[i] < EXTINCTION_THRESHOLD
            duH[i] = H[i] < 0.0 ? abs(H[i]) : -H[i]
        elseif iszero(H[i])
            duH[i] = 0.0
        end
    end

    # Predator dynamics
    for α in 1:R
        if P[α] > EXTINCTION_THRESHOLD
            # Compute total effective attack from herbivores.
            attack_sum = 0.0
            for i in 1:S
                # Using the conversion efficiency to scale the raw attack rate.
                attack_sum += epsilon[α] * a_matrix[i, α] * H[i]
            end

            # Sum of predator–predator interactions.
            interact_sum = 0.0
            for β in 1:R
                interact_sum += B[α, β] * P[β]
            end

            # Assume reference density P₀ = m_alpha[α] (i.e. d_α = 1).
            P0 = m_alpha[α]
            duP[α] = P[α] * m_alpha[α] * ( (attack_sum - P[α] + interact_sum) / P0 - 1 )
            # Equivalently, this is:
            # duP[α] = P[α] * (attack_sum - P[α] + interact_sum - m_alpha[α])
        elseif P[α] < EXTINCTION_THRESHOLD
            duP[α] = P[α] < 0.0 ? abs(P[α]) : -P[α]
        elseif iszero(P[α])
            duP[α] = 0.0
        end
    end

    # Write back the derivatives.
    for i in 1:S
        du[i] = duH[i]
    end
    for α in 1:R
        du[S+α] = duP[α]
    end
end
