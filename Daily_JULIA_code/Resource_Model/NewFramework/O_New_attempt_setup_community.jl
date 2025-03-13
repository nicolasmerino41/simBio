function o_new_attempt_setup_community(
    i, j, 
    mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha; 
    species_names = nothing, artificial_pi = false,
    o_pred_avg::Float64 = 0.0,
)
    try
        S, R, H_i0_eff, m_i, g_i, O_loss, O_gain, A_star,
        a_matrix, A, epsilon, m_alpha, herbivore_list, predator_list, species_names =
         o_new_setup_community_from_cell(
                i, j;
                M_mean = 0.1,
                mu = mu_val,
                symmetrical_competition = sym_competition,
                mean_m_alpha = mean_m_alpha,
                epsilon_val = eps_val,
                mu_predation = mu_pred_val,
                iberian_interact_NA = iberian_interact_NA,
                species_dict = species_dict,
                species_names = species_names,
                artificial_pi = artificial_pi,
                o_pred_avg = o_pred_avg
            )
        return (
            S, R, H_i0_eff, m_i, g_i, O_loss, O_gain, A_star,
            a_matrix, A, epsilon, m_alpha, herbivore_list, predator_list, species_names
        )
    catch e
        println("hey error in o_new_attempt_setup_community")
        # println(e)
        return nothing
    end
end