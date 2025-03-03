using NLsolve, LinearAlgebra, ForwardDiff, DataFrames

#--------------------------------------------------------------------
# Define the equilibrium system for herbivores only.
# We assume the equilibrium condition is:
#
#    F_i(H) = H_i + sum_j M_mod[i,j]*H_j - H0[i] = 0,
#
# where H0 is the effective observed herbivore abundance (without predator effects)
# and M_mod is the (modified) competition matrix.
#--------------------------------------------------------------------
function herbivore_equilibrium_system!(F, H, params)
    M_mod, H0 = params  # M_mod is an S×S matrix, H0 is a vector of length S
    S = length(H0)
    for i in 1:S
        comp_sum = 0.0
        for j in 1:S
            comp_sum += M_mod[i, j] * H[j]
        end
        F[i] = H[i] + comp_sum - H0[i]
    end
end

#--------------------------------------------------------------------
# Compute the Jacobian of the equilibrium system numerically using ForwardDiff.
#--------------------------------------------------------------------
function compute_jacobian(f!, H, params)
    return ForwardDiff.jacobian(x -> begin
        F = similar(x)
        f!(F, x, params)
        F
    end, H)
end

#--------------------------------------------------------------------
# Define a function to run the herbivore-only equilibrium analysis for a given cell 
# and parameter configuration. Here we only vary mu while keeping NPP fixed.
#
# The function calls new_setup_community_from_cell (which returns:
# (S, R, H_i0_eff, m_i, g_i, beta, G_raw, M_modified, A_star, a_matrix, A, epsilon, m_alpha, x, raw_g))
# and then uses NLsolve to solve:
#
#    H + M_mod * H - H0 = 0.
#
# It then computes the Jacobian and checks its eigenvalues.
#--------------------------------------------------------------------
function explore_herbivore_equilibrium(cell::Int; 
    NPP::Float64 = 1000.0,
    mu_range = range(0.0, stop=1.0, length=20),  # vary mu over 20 points
    M_mean::Float64 = 0.1,
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,    # not used if predators are excluded
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.0,    # herbivore-only, so no predation effect
    artificial_pi = false,
    alpha::Float64 = 0.25,
    include_predators::Bool = true
)
    
    results_table = DataFrame(mu = Float64[], f_converged = Bool[], x_converged = Bool[], stable = Bool[], max_eig = Float64[])
    
    # @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]

    # Extract species names from the cell data.
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    
    if !include_predators 
        filter!(name -> !(name in predator_names), sp_nm)
    end

    for mu_val in mu_range
        # For a herbivore-only model, force predators to zero:
        # (This is done by filtering species_names in herbivore_run)
        results = new_attempt_setup_community(
            idx[cell][1], idx[cell][2],
            mu_val, mu_predation, epsilon_val, symmetrical_competition, mean_m_alpha;
            localNPP = NPP,
            species_names = sp_nm,  # let new_setup_community_from_cell extract from the cell
            artificial_pi = artificial_pi,
            alpha = alpha
        )
        # If results is nothing, skip this parameter set.
        if results === nothing
            continue
        end
        
        S = results.S
        if !include_predators
            R = results.R
            # For herbivore-only, we expect R == 0; 
            # however, if R > 0 we can simply ignore predator terms.
            if R > 0
                @warn "Herbivore-only run: ignoring R = $R predator species."
            end
        end
        
        # Extract number of species and effective herbivore abundance.
        S_val = results.S  # Number of herbivores
        # In a herbivore-only run, we expect R == 0.
        H0_eff = results.H_i0   # Effective observed herbivore abundances (vector of length S)
        M_mod = results.M_modified

        # Destructure the returned NamedTuple.
        S2          = results.S
        R2          = results.R
        H_i0        = results.H_i0
        m_i         = results.m_i
        g_i         = results.g_i
        beta        = results.beta
        M_mod       = results.M_modified    # Modified competition matrix (S×S)
        A_star      = results.A_star        # Nondimensional predation rates (S×R)
        a_matrix    = results.a_matrix      # Herbivore–predator interaction matrix (S×R)
        A           = results.A             # Predator interaction matrix (R×R)
        m_alpha     = results.m_alpha

        # Build predator attack matrix (R×S) by transposing a_matrix.
        A_pred = transpose(a_matrix)

        # Set predator overpopulation thresholds.
        # (Here we assume P0 = m_alpha, which implies a self-regulation coefficient of 1.)
        P0 = m_alpha
        params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
        # Use H0_eff as the initial guess for the equilibrium.
        x0 = float.(vcat(H0_eff, H0_eff[1:R2].*0.1))

        # Solve the equilibrium system: H + M_mod*H - H0_eff = 0.
        sol = nlsolve((F, x) -> equilibrium_system!(F, x, params), x0)
        # Check convergence:
        if !sol.f_converged
            @warn "nlsolve did not converge for mu = $mu_val"
            continue
        end
        
        H_eq = sol.zero
        
        # Compute the Jacobian at equilibrium.
        J = compute_jacobian(equilibrium_system!, H_eq, params)
        eigvals_J = eigvals(J)
        max_real = maximum(real.(eigvals_J))
        stable = max_real < 0.0
        push!(results_table, (mu_val, sol.f_converged, sol.x_converged, stable, max_real))
    end
    
    return results_table
end

#--------------------------------------------------------------------
# Now, run the parameter exploration for a given cell (say cell 1)
#--------------------------------------------------------------------
AAAA = explore_herbivore_equilibrium(
    1;
    NPP=Float64(npp_DA_relative_to_1000[idx[1][1], idx[1][2]]),
    mu_range=range(0.0, stop=1.0, length=50),
    M_mean=0.1, symmetrical_competition=true, mean_m_alpha=0.0000001, epsilon_val=1.0,
    mu_predation=0.01, artificial_pi=false, alpha=0.25,
    include_predators = false
)

println("Parameter exploration results:")
println(results_table)
