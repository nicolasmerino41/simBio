function o_attempt_setup_community(
    i, j, 
    mu_val, eps_val, mean_m_alpha;
    species_names = nothing,
    artificial_pi = false, pi_size = 1.0,
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    r_omni_proportion = 1.0
)
    # Convert parameters:
    mu_val       = to_float(mu_val)
    eps_val      = to_float(eps_val)
    mean_m_alpha = to_float(mean_m_alpha)
    pi_size      = to_float(pi_size)
    delta_nu   = to_float(delta_nu)
    d_alpha    = to_float(d_alpha)
    d_i        = to_float(d_i)
    r_omni_proportion = to_float(r_omni_proportion)
    
    try
        params = o_setup_community_from_cell(
            i, j;
            mu = mu_val,
            mean_m_alpha = mean_m_alpha,
            epsilon_val = eps_val,
            iberian_interact_NA = iberian_interact_NA,
            species_dict = species_dict,
            species_names = isnothing(species_names) ? String[] : species_names,
            artificial_pi = artificial_pi, pi_size = pi_size,
            delta_nu = delta_nu,
            d_alpha = d_alpha,
            d_i = d_i,
            r_omni_proportion = r_omni_proportion
        )
        return (
            S = params.S, R = params.R,
            H_eq = params.H_eq, P_eq = params.P_eq,
            r_i = params.r_i, K_i = params.K_i,
            mu = params.mu, nu = params.nu, 
            P_matrix = params.P_matrix, O_matrix = params.O_matrix, T_matrix = params.T_matrix,
            B_matrix = params.B_matrix, D_matrix = params.D_matrix,
            epsilon = params.epsilon, m_alpha = params.m_alpha,
            K_alpha = params.K_alpha,
            herbivore_list = params.herbivore_list, predator_list = params.predator_list,
            species_names = params.species_names,
            H_star = params.H_star, P_star = params.P_star
        )
    catch e
        println("hey error in o_attempt_setup_community")
        println(e)
        return nothing
    end
end