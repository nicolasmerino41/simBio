function o_new_attempt_setup_community(
    i, j, 
    mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha; 
    localNPP,   # not used in new framework
    species_names = nothing, artificial_pi = false,
    o_pred_avg::Float64 = 1.0,
)
    try
        S2, R2, species_names, herb_list, pred_list,
        H_i0, m_i, g_i, x_final, beta,
        O_loss, O_gain, A_star, a_matrix, A,
        epsilon, m_alpha, raw_g, h = o_new_setup_community_from_cell(
                i, j;
                M_mean = 0.1,
                mu = mu_val,
                symmetrical_competition = sym_competition,
                mean_m_alpha = mean_m_alpha,
                epsilon_val = eps_val,
                mu_predation = mu_pred_val,
                iberian_interact_NA = iberian_interact_NA,
                species_dict = species_dict,
                cell_abundance = [],  # assume cell abundances come from cell data
                o_pred_avg = o_pred_avg
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
            O_loss = O_loss,
            O_gain = O_gain,
            A_star = A_star,
            a_matrix = a_matrix,
            A = A,
            epsilon = epsilon,
            m_alpha = m_alpha,
            raw_g = raw_g,
            h = h
        )
    catch e
        println("hey error in o_new_attempt_setup_community")
        return nothing
    end
end