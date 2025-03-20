# using NLsolve

# Define the equilibrium function.
# "params" is a tuple: (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
function equilibrium_system!(F, x, params)
    # Unpack the first 12 parameters:
    # S       : number of herbivore species
    # R       : number of predator species
    # H_i0    : effective observed herbivore abundances (vector of length S)
    # m_i     : herbivore mortality rates (vector, length S)
    # g_i     : effective herbivore growth rates (vector, length S)
    # beta    : niche parameters for herbivores (vector, length S)
    # M_mod   : modified competition matrix among herbivores (S×S)
    # A_star  : nondimensional attack rates for herbivores (S×R)
    # A_pred  : predator attack rate matrix (R×S)
    # P0      : predator overpopulation thresholds (vector of length R)
    # A       : predator interaction matrix (R×R)
    # m_alpha : predator mortality rates (vector of length R)
    S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha = params[1:12]
    
    # Check for extra parameter h (handling time) for Holling type II.
    holling = false
    h = 0.0
    if length(params) >= 13
        h = params[13]
        holling = true
    end

    # x[1:S] are equilibrium herbivore abundances; x[S+1:S+R] are predator abundances.
    H = x[1:S]
    P = x[S+1:S+R]
    
    # Herbivore equilibrium equations:
    for i in 1:S
        if H[i] > EXTINCTION_THRESHOLD
            comp_term = 0.0
            for j in 1:S
                comp_term += M_mod[i, j] * H[j]
            end
            # Compute the linear predator contribution.
            linear_pred_term = 0.0
            for α in 1:R
                linear_pred_term += A_star[i, α] * P[α]
            end
            # If Holling type II is active, saturate the predator contribution.
            pred_term = holling ? (linear_pred_term / (1 + h * linear_pred_term)) : linear_pred_term
            F[i] = H[i] + comp_term + pred_term - H_i0[i]
        else
            F[i] = 0.0
        end
    end

    # Predator equilibrium equations:
    for α in 1:R
        if P[α] > EXTINCTION_THRESHOLD
            attack_sum = 0.0
            for i in 1:S
                attack_sum += A_pred[α, i] * H[i]
            end
            interact_sum = 0.0
            for β in 1:R
                interact_sum += A[α, β] * P[β]
            end 
            F[S + α] = attack_sum - m_alpha[α] - interact_sum
        else
            F[S + α] = 0.0
        end
    end
end

# Function to run the equilibrium solver using parameters from herbivore_run.
# Here we assume that "herbivore_run" returns both the results DataFrame and the parameters tuple.
function herbivore_equilibrium_run(
    cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
    include_predators=true,
    plot=false, sp_removed_name=nothing,
    artificial_pi=false, NPP=nothing,
    alpha = 0.25, ignore_inf_error=false,
    hollingII=false, h=0.1
)
    # First, run your herbivore_run function to get the parameters.
    results, params_from_run = herbivore_run(
        cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
        include_predators=include_predators,
        plot=plot, sp_removed_name=sp_removed_name, NPP=NPP, artificial_pi=artificial_pi,
        alpha=alpha, do_you_want_params=true,
        ignore_inf_error=ignore_inf_error,
        hollingII=hollingII, h=h
    )
    
    if isnothing(params_from_run)
        @error "Herbivore run failed. Cannot compute equilibrium."
        return nothing, nothing, nothing
    end
    
    # Here, "params_from_run" is assumed to be the tuple:
    # (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
    # In our previous code, we had:
    # params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
    # We'll use that for our equilibrium system.
    
    # For the equilibrium solver we need an initial guess.
    # A natural initial guess is to start at the effective herbivore abundance for H,
    # and for predators, a small positive value.
    S = params_from_run[1]
    R = params_from_run[2]
    H0 = params_from_run[3]
    # Use H0 as the initial guess for herbivores, and a small constant (e.g., 0.1) for predators.
    x0 = vcat(H0, fill(0.1, R))
    
    # Now call the root finder.
    sol = nlsolve((F, x) -> equilibrium_system!(F, x, params_from_run), x0)
    
    if sol.f_converged && any(sol.zero[S+1:end] .< 0.0)
        H_eq = sol.zero[1:S]
        P_eq = sol.zero[S+1:end]
        @info "Converged with negative predator abundances"
        # println("Equilibrium herbivore abundances: ", H_eq)
        println("Equilibrium predator abundances: ", P_eq)
        println("Equilibrium herbivore abundances: ", H_eq)
        return H_eq, P_eq, params_from_run
    elseif sol.f_converged && !any(sol.zero[S+1:end] .< 0.0)
        @info "Converged without negative predator abundances"
        H_eq = sol.zero[1:S]
        P_eq = sol.zero[S+1:end]
        println("Equilibrium predator abundances: ", P_eq)
        println("Equilibrium herbivore abundances: ", H_eq)

        return H_eq, P_eq, params_from_run
    else
        @error "Equilibrium solver did not converge."
        return nothing, nothing, params_from_run
    end
end

# begin
    
#     mu = 0.0
#     mu_pred = 0.01
#     epsilon = 1.0
#     m_alpha = 0.01
#     alpha = 0.25
#     h = 0.1
#     deviation_from_equilibrium = 0.0

#     H_eq, P_eq, params = herbivore_equilibrium_run(
#         1, mu, mu_pred, epsilon, true, m_alpha;
#         include_predators=true,
#         plot=false, sp_removed_name=nothing,
#         NPP=nothing, artificial_pi=false,
#         alpha=alpha,
#         ignore_inf_error=false,
#         hollingII=true, h=h
#     )
    
#     if !isnothing(H_eq) || !isnothing(P_eq)
#         if all(H_eq .> 0.0) && all(P_eq .> 0.0)
#             @info "Equilibrium exists with non-negative herbivores and predators."
#         elseif all(H_eq .> 0.0) && !any(P_eq .> 0.0)
#             @info "Equilibrium exists with non-negative herbivores and negative predators."
#         elseif !any(H_eq .> 0.0) && all(P_eq .> 0.0)
#             @info "Equilibrium exists with negative herbivores and non-negative predators."
#         elseif !any(H_eq .> 0.0) && !any(P_eq .> 0.0)
#             @info "Equilibrium exists with negative herbivores and predators."
#         else
#             @info "Equilibrium exists with mixed signs (some positive, some negative) herbivores and predators."
#         end   
#         H_eq[H_eq .< 0.0] .= 0.0
#         P_eq[P_eq .< 0.0] .= 0.0
#         AAAA = run_right_away(
#             params;
#             time_end=500.0,
#             plot=true,
#             H_init=H_eq*(1+deviation_from_equilibrium), P_init=P_eq*(1+deviation_from_equilibrium)
#             )
#         println("Herbivore abundances: ", AAAA.H_vector[1])
#         println("Predator abundances: ", AAAA.P_vector[1])
#     end
# end