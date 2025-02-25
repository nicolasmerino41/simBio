function new_attempt_setup_community(
    i, j, 
    mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha; 
    localNPP,
    # localH0_vector, 
    species_names = nothing, artificial_pi = false,
    alpha = 0.25
)
    try
        S2, R2, species_names, herb_list, pred_list,
        H_i0, m_i, g_i, x_final, beta,
        G, M_modified, A_star, a_matrix, A,
        epsilon, m_alpha, raw_g =
            new_setup_community_from_cell(
                i, j;
                NPP = localNPP,
                M_mean = 0.1,
                mu = mu_val,
                symmetrical_competition = sym_competition,
                mean_m_alpha = mean_m_alpha,
                epsilon_val = eps_val,
                mu_predation = mu_pred_val,
                iberian_interact_NA = iberian_interact_NA,
                species_dict = species_dict,
                m_standard_deviation = 0.0,
                h_standard_deviation = 0.0,
                artificial_pi = artificial_pi,
                # real_H0 = true,
                # H0_vector = localH0_vector,
                species_names = species_names,
                alpha = alpha
            )
        return (
            S = S2, 
            R = R2,
            species_names = species_names,
            herbivore_list = herb_list,
            predator_list = pred_list,
            H_i0 = H_i0,
            m_i = m_i,
            g_i = g_i,
            x_final = x_final,
            beta = beta,
            G = G,
            M_modified = M_modified,
            A_star = A_star,
            a_matrix = a_matrix,
            A = A,
            epsilon = epsilon,
            m_alpha = m_alpha,
            raw_g = raw_g
        )
    catch e
        # Optionally log the error here.
        # @warn "Failed in new_attempt_setup_community: $e\n$(catch_backtrace())"
        return nothing
    end
end