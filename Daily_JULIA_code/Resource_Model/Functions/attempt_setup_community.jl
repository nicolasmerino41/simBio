function attempt_setup_community(
    i, j, 
    mu_val, mu_pred_val, eps_val, sym_competition; 
    localNPP = nothing, localH0_vector = nothing, 
    species_names=nothing, artificial_pi=false,
    herbivore_m = 0.1,
    predator_m = 0.1
)
    # println("i = $i, j = $j", "mu_val = $mu_val, mu_pred_val = $mu_pred_val, eps_val = $eps_val, localNPP = $localNPP, localH0_vector = ", length(localH0_vector), "art_pi = $artificial_pi")
    try
        S2, R2,
        species_names, herb_list, pred_list,
        H_i0, m_i,
        p_vec, x_final, g_i,
        localHatH, G,
        M_modified, a_matrix, A,
        epsilon_vector, m_alpha = setup_community_from_cell(
            i, j;
            NPP = localNPP,
            M_mean = herbivore_m,
            mu = mu_val,
            symmetrical_competition = sym_competition,
            mean_m_alpha = predator_m,
            epsilon_val = eps_val,
            mu_predation = mu_pred_val,
            iberian_interact_NA = iberian_interact_NA,
            species_dict = species_dict,
            m_standard_deviation = 0.0,
            h_standard_deviation = 0.0,
            artificial_pi = artificial_pi,
            real_H0 = true,
            H0_vector = localH0_vector,
            species_names = species_names
        )
        # If success, return NamedTuple
        return (
            S=S2, 
            R=R2,
            species_names=species_names,
            herbivore_list=herb_list,
            predator_list=pred_list,
            H_i0=H_i0,
            m_i=m_i,
            p_vec=p_vec,
            x_final=x_final, 
            g_i=g_i,
            localHatH=localHatH,
            G=G,
            M_modified=M_modified,
            a_matrix=a_matrix,
            A=A,
            epsilon_vector=epsilon_vector,
            m_alpha=m_alpha
        )
    catch e
        # If anything fails, we log it or ignore it, then return nothing
        # @warn "Failed inside attempt_setup_community with error: $e\n$(catch_backtrace())"
        return nothing
    end
end