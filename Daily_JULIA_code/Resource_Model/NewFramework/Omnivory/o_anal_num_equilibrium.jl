using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots

function analytical_equilibrium(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    nu_omni_proportion = 1.0 
)
    local_i, local_j = idx[cell][1], idx[cell][2]
    # Extract species names from the cell.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species: if predators are not included, restrict to herbivores & omnivores.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names || sp in omnivore_names]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Use the new community setup function.
    params_setup = o_attempt_setup_community(
        local_i, local_j,
        mu_val, eps_val, mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure parameters from the new setup.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq
    P_eq       = params_setup.P_eq
    r_i        = params_setup.r_i
    K_i        = params_setup.K_i
    nu_val     = params_setup.nu
    P_matrix   = params_setup.P_matrix
    O_matrix   = params_setup.O_matrix
    T_matrix   = params_setup.T_matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    B_matrix   = params_setup.B_matrix
    D_matrix   = params_setup.D_matrix
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    
    # Set initial conditions.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    println("typeof(u0): ", typeof(u0))
    
    # Build the parameter tuple for omnivore_dynamics!
    # Order: (S, R, K_i, r_i, mu, nu, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_omni)
    p = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, O_matrix, T_matrix, epsilon, m_alpha, K_alpha, B_matrix, D_matrix, nu_val*nu_omni_proportion)
    
    # Define an inner wrapper so we can compute the Jacobian at u0.
    function f_wrapper(u)
        du = similar(u)
        omnivore_dynamics!(du, u, p, 0.0)
        return du
    end
    
    J = ForwardDiff.jacobian(f_wrapper, u0)
    println("Jacobian at equilibrium:")
    println(J)
    
    # Package results into a NamedTuple.
    return (
        equilibrium = (H_star = H_star, P_star = P_star, u0 = u0),
        parameters  = (
            S = S,
            R = R,
            r_i = r_i,
            K_i = K_i,
            mu = mu_val,
            nu = nu_val,
            P_matrix = P_matrix,
            O_matrix = O_matrix,
            T_matrix = T_matrix,
            epsilon = epsilon,
            m_alpha = m_alpha,
            K_alpha = K_alpha,
            B_matrix = B_matrix,
            D_matrix = D_matrix,
            nu_omni = nu_val*nu_omni_proportion
        ),
        Jacobian = J
    )
end

# --- Compute equilibrium and assess local stability ---
result = analytical_equilibrium(
    1, 
    0.5, 1.0, 0.1;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    nu_omni_proportion = 1.0
)

A_eq = result.equilibrium       # Contains H_star, P_star, and u0 (equilibrium state)
A_p  = result.parameters        # Contains the ODE parameters
A_j  = result.Jacobian          # The Jacobian at equilibrium

println("Jacobian at equilibrium:")
println(A_j)

# Compute eigenvalues to assess local stability.
eigvals = eigen(A_j).values
println("Eigenvalues:")
println(eigvals)
stable = all(real.(eigvals) .< 0)
println("Is the equilibrium locally stable? ", stable)

