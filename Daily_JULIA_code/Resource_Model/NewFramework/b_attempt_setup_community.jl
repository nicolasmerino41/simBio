function b_attempt_setup_community(
    i, j, 
    mu_val, eps_val, mean_m_alpha;
    species_names = nothing, artificial_pi = false
)
    try
        params = b_new_setup_community_from_cell(
            i, j;
            mu = mu_val,
            # nu = nu_val,
            mean_m_alpha = mean_m_alpha,
            epsilon_val = eps_val,
            iberian_interact_NA = iberian_interact_NA,
            species_dict = species_dict,
            species_names = isnothing(species_names) ? String[] : species_names,
            artificial_pi = artificial_pi
        )
        return (
            S = params.S, R = params.R,
            H_i0 = params.H_i0, r_i = params.r_i, K_i = params.K_i,
            mu = params.mu, nu = params.nu, 
            P_matrix = params.P_matrix, epsilon = params.epsilon, m_alpha = params.m_alpha,
            K_alpha = params.K_alpha,
            herbivore_list = params.herbivore_list, predator_list = params.predator_list,
            species_names = params.species_names,
            H_star = params.H_star, P_star = params.P_star
        ) 
    catch e
        println("hey error in b_attempt_setup_community")
        return nothing
    end
end
