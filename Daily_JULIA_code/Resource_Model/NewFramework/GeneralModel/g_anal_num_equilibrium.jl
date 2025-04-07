execute_code = false
function g_analytical_equilibrium(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    include_omnivores = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 1.0,
    callbacks = false,
    plot = false
)

    # Convert parameters to Float64:
    mu_val         = to_float(mu_val)
    eps_val        = to_float(eps_val)
    mean_m_alpha   = to_float(mean_m_alpha)
    delta_nu       = to_float(delta_nu)
    d_alpha        = to_float(d_alpha)
    d_i            = to_float(d_i)
    pi_size        = to_float(pi_size)
    nu_omni_proportion = to_float(nu_omni_proportion)
    nu_b_proportion    = to_float(nu_b_proportion)
    r_omni_proportion  = to_float(r_omni_proportion)

    local_i, local_j = idx[cell][1], idx[cell][2]
    # Extract species names from the cell.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species: if predators are not included, restrict to herbivores & omnivores.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names || sp in omnivore_names]
    end
    # Process omnivores:
    if !include_omnivores
        # Exclude all omnivores.
        sp_nm = [sp for sp in sp_nm if !(sp in omnivore_names)]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Use the new community setup function.
    params_setup = g_attempt_setup_community(
        local_i, local_j,
        mu_val, eps_val, mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end
    
    # Destructure the returned parameters.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq        # Baseline herbivore (and omnivore) abundances
    P_eq       = params_setup.P_eq        # Baseline predator abundances
    r_i        = params_setup.r_i          # Herbivore intrinsic growth rates
    K_i        = params_setup.K_i          # Herbivore carrying capacities
    mu_val     = params_setup.mu           # Herbivore competition
    nu_val     = params_setup.nu           # Effective predation rate
    A_matrix   = params_setup.A_matrix     # R x S omnivore-herbivore interaction matrix
    C_matrix   = params_setup.C_matrix     # R x S omnivore-omnivore interaction matrix
    epsilon    = params_setup.epsilon      # Predator conversion efficiencies
    m_alpha    = params_setup.m_alpha      # Predator mortality rates
    K_alpha    = params_setup.K_alpha      # Predator carrying capacities
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    predator_list   = params_setup.predator_list
    
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
    # println("typeof(u0): ", typeof(u0))
    
    # Build the parameter tuple for general_dynamics!
    # Order: (S, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)
    
    p = (S, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)
    
    # Define an inner wrapper so we can compute the Jacobian at u0.
    function f_wrapper(u)
        du = similar(u)
        general_dynamics!(du, u, p, 0.0)
        return du
    end
    
    J = ForwardDiff.jacobian(f_wrapper, u0)
    # println("Jacobian at equilibrium:")
    # println(J)
    
    # Define and solve the ODE using general_dynamics!
    prob = ODEProblem(general_dynamics!, u0, (0.0, 1000.0), p)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if callbacks
            solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
        else
            solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    end

    # --- Plotting the Dynamics using Makie ---
    if plot    
        fig = Figure(; size = (1000, 600))
        ax = Axis(
            fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for cell"
        )
        times_combined = sol.t
        # Plot herbivore/omnivore dynamics: species that are in herbivore_names are blue; those in omnivore_names are green.
        for i in 1:S
            if herbivore_list[i] in omnivore_names
                lines!(ax, times_combined, sol[i, :], label = "O$(i)", color = :green)
                # println("O$(i)")
            else
                lines!(ax, times_combined, sol[i, :], label = "H$(i)", color = :blue)
                # println("H$(i)")
            end
        end
        # Plot predator dynamics (indices S+1:S+R) in red dashed lines.
        for alpha in 1:R
            lines!(ax, times_combined, sol[S+alpha, :], label = "P$(alpha)", linestyle = :dash, color = :red)
        end
        display(fig)
    end

    # Package results into a NamedTuple.
    return (
        equilibrium = (
            H_star = H_star, P_star = P_star,
            u0 = u0,
            H_vector = sol[1:S, end],
            P_vector = sol[(S+1):(S+R), end],
            herbivore_list = herbivore_list,
            predator_list = predator_list
            ),
        parameters  = (
            S = S,
            R = R,
            K_i = K_i,
            r_i = r_i,
            mu = mu_val,
            nu = nu_val,
            A_matrix = A_matrix,
            C_matrix = C_matrix,
            m_alpha = m_alpha,
            K_alpha = K_alpha
        ),
        Jacobian = J
    )
end

# --- Compute equilibrium and assess local stability ---
if execute_code
#     for j in 0.0:0.1:1.0
#         for i in 0.0:0.1:1.0
#             for k in 0.1:0.1:0.3
                begin
                i,j,k = 0.2, 0.5, 0.1 
                result = g_analytical_equilibrium(
                    1, 
                    i, j, k; # mu, eps, m_alpha
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = true, pi_size = 10.0,
                    H_init = nothing,
                    P_init = nothing,
                    r_omni_proportion = 1.0,
                    callbacks = false,
                    plot = true
                );

                A_eq = result.equilibrium       # Contains H_star, P_star, and u0 (equilibrium state)
                A_p  = result.parameters        # Contains the ODE parameters
                A_j  = result.Jacobian          # The Jacobian at equilibrium

                # println("Jacobian at equilibrium:")
                # println(A_j)
                if any(isnan, A_j) || any(isinf, A_j)
                    println("Jacobian contains NaN or Inf values for mu = $(i), eps = $(j), m_alpha = $(k)")
                else
                    # Compute eigenvalues to assess local stability.
                    eigvals = eigen(A_j).values
                    # println("Eigenvalues:")
                    # println(eigvals)
                    stable = all(real.(eigvals) .< 0)
                    if stable
                        println("STABLEEEE")# for mu = $(i), eps = $(j), m_alpha = $(k)")
                        sr = (length(A_eq.H_vector[A_eq.H_vector .> EXTINCTION_THRESHOLD])+ length(A_eq.P_vector[A_eq.P_vector .> EXTINCTION_THRESHOLD]))/(A_p.S+A_p.R)
                        println("SR = $(sr)")
                    elseif !stable && isapprox(sum(A_eq.H_star), sum(A_eq.H_vector), atol = 1e-3) && isapprox(sum(A_eq.P_star), sum(A_eq.P_vector), atol = 1e-3)
                        println("Unstable but at equilibrium")#  for mu = $(i), eps = $(j), m_alpha = $(k)")
                        println("H_star_biomass = ", sum(A_eq.H_star), " and final biomass = ", sum(A_eq.H_vector))
                    elseif !stable
                        println("Unstable")#for mu = $(i), eps = $(j), m_alpha = $(k)")
                        println("H_star_biomass = ", sum(A_eq.H_star), " and final biomass = ", sum(A_eq.H_vector))
                    end
                end
                end
#             end
#         end
#     end
end
