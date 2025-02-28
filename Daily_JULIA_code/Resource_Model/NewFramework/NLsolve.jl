using NLsolve

# Define the equilibrium function.
# "params" is a tuple: (S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
function equilibrium_system!(F, x, params)
    # Unpack parameters
    S, R, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha = params
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
            pred_term = 0.0
            for α in 1:R
                pred_term += A_star[i, α] * P[α]
            end
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
    alpha = 0.25, disable_inf_error=false
)
    # First, run your herbivore_run function to get the parameters.
    results, params_from_run = herbivore_run(
        cell, mu_val, mu_pred_val, eps_val, sym_competition, mean_m_alpha;
        include_predators=include_predators,
        plot=plot, sp_removed_name=sp_removed_name, NPP=NPP, artificial_pi=artificial_pi,
        alpha=alpha, do_you_want_params=true,
        disable_inf_error=disable_inf_error
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
        return H_eq, P_eq, params_from_run
    elseif sol.f_converged && !any(sol.zero[S+1:end] .< 0.0)
        @info "Converged without negative predator abundances"
        H_eq = sol.zero[1:S]
        P_eq = sol.zero[S+1:end]
        println("Equilibrium predator abundances: ", P_eq)

        return H_eq, P_eq, params_from_run
    else
        @error "Equilibrium solver did not converge."
        return nothing, nothing, params_from_run
    end
end

begin
    
    mu = 0.5
    mu_pred = 0.0005
    epsilon = 0.05
    m_alpha = 0.5
    alpha = 0.25
    deviation_from_equilibrium = 0.0

    H_eq, P_eq, params = herbivore_equilibrium_run(
        1, mu, mu_pred, epsilon, true, m_alpha;
        include_predators=true,
        plot=false, sp_removed_name=nothing, 
        NPP=nothing, artificial_pi=false,
        alpha = alpha,
        disable_inf_error=false
    )
    if !isnothing(H_eq) || !isnothing(P_eq)
        AAAA = run_right_away(
            params;
            time_end=500.0,
            plot=true,
            H_init=H_eq*(1+deviation_from_equilibrium), P_init=P_eq*(1+deviation_from_equilibrium)
        )
        println("Herbivore abundances: ", AAAA.H_vector[1])
        println("Predator abundances: ", AAAA.P_vector[1])
    end
end