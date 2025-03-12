function omnivore_dynamics!(du, u, p, t)
    # p = (S, R, H_i0_eff, m_i, g_i, beta, O_loss, O_gain, A_loss, B, m_alpha)
    S, R, H_i0_eff, m_i, g_i, beta, O_loss, O_gain, A_loss, B, m_alpha = p
    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)
    
    # Herbivore dynamics:
    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            H_tot = sum(H)
            baseline = (1 - mu)*H[i] + mu*H_tot
            loss_herb = 0.0
            gain_herb = 0.0
            for j in 1:S
                if i != j
                    loss_herb += O_loss[i,j] * H[j]
                    gain_herb += O_gain[i,j] * H[j]
                end
            end
            loss_pred = 0.0
            for alpha in 1:R
                loss_pred += A_loss[i, alpha] * P[alpha]
            end
            effective_density = baseline - loss_herb + gain_herb - loss_pred
            duH[i] = H[i] * g_i[i] * (1 - effective_density / H_i0_eff[i])
        else
            duH[i] = 0.0
        end
    end
    
    # Predator dynamics (linear):
    for alph in 1:R
        if P[alph] > EXTINCTION_THRESHOLD
            attack_sum = 0.0
            for i in 1:S
                attack_sum += A[alph, alph] * H[i]  # Here A[alph,alph] = 1
            end
            interact_sum = 0.0
            for bet in 1:R
                interact_sum += B[alph, bet] * P[bet]
            end
            duP[alph] = P[alph] * m_alpha[alph] * ((attack_sum - interact_sum) / m_alpha[alph])
        else
            duP[alph] = 0.0
        end
    end
    
    du[1:S] = duH
    if R > 0
        du[S+1:S+R] = duP
    end
end