using NLsolve, LinearAlgebra, ForwardDiff, DataFrames, Statistics

#--------------------------------------------------------------------
# Define the equilibrium system for herbivores only.
# We assume the equilibrium condition:
#
#    F(H) = H + M_mod * H - H0 = 0,
#
# where H0 is the effective observed herbivore abundance and M_mod is the modified competition matrix.
#--------------------------------------------------------------------
function herbivore_equilibrium_system!(F, H, params)
    M_mod, H0 = params  # M_mod is S×S, H0 is vector of length S
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
# Define a function to check the stability of the herbivore equilibrium from an ODE solution.
# We will consider the system practically stable if the maximum relative change in the herbivore 
# densities over the last window_fraction of the simulation is below tol.
#--------------------------------------------------------------------
function is_stable(sol; window_fraction=0.1, tol=0.01)
    t_final = sol.t[end]
    window_start = t_final * (1 - window_fraction)
    idx_window = findall(t -> t >= window_start, sol.t)
    # Determine the number of herbivore species from the first state vector.
    S = length(sol.u[1])
    # Build a matrix whose columns are the herbivore states at the time indices in idx_window.
    H_mat = hcat(sol.u...)[1:S, idx_window]  # each column is H at a time point in the window
    # For each herbivore species, compute the maximum relative difference between successive time points.
    rel_diffs = [ maximum(abs.(diff(H_mat[i, :]) ./ H_mat[i, 1:end-1])) for i in 1:S ]
    return maximum(rel_diffs) < tol, maximum(rel_diffs)
end

function is_stable_by_nico(sol)
    final_abundances = sol[1:end, end]
    if any(final_abundances .> 100.0)
        return false
    else
        return true
    end
end

#--------------------------------------------------------------------
# Define a function to explore the stability of herbivore-only equilibria
# for a given cell, sweeping over a range of mu values.
#
# For herbivore-only runs, we exclude predators.
#--------------------------------------------------------------------
function explore_stability(
    cell::Int; 
    NPP::Float64 = 1000.0,
    mu_range = range(0.0, stop=1.0, length=20),
    M_mean::Float64 = 0.1,
    symmetrical_competition::Bool = true,
    mean_m_alpha::Float64 = 0.1,    # not used if predators are excluded
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.0,    # no predation
    artificial_pi = false,
    alpha::Float64 = 0.25,
    include_predators = false,
    plot = false # Beware plotting will be very slow for more than a few iterations
)
    results_table = DataFrame(mu = Float64[], f_converged = Bool[], x_converged = Bool[],
                              stable = Bool[], max_rel_change = Float64[], stable_by_nico = Bool[],
                              survival_rate = Float64[],
                              mu_predation = Float64[], epsilon_val = Float64[])
    
    # Get cell indices and species names.
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    if include_predators
        # Exclude predators for herbivore-only runs.
        filter!(name -> !(name in predator_names), sp_nm)
    end
    
    for mu_val in mu_range
        # Run the community setup with current mu (others fixed)
        results = new_attempt_setup_community(
            local_i, local_j,
            mu_val, mu_predation, epsilon_val, symmetrical_competition, mean_m_alpha;
            localNPP = NPP,
            species_names = sp_nm,
            artificial_pi = artificial_pi,
            alpha = alpha
        )
        if results === nothing
            continue
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
        # Extract equilibrium:
        H_eq = sol.zero
        # Compute numerical Jacobian and its eigenvalues.
        J = compute_jacobian(equilibrium_system!, H_eq, params)
        max_real = maximum(real.(eigvals(J)))
        stable_flag = max_real < 0.0
        # Now, run a full simulation of the herbivore-only ODE (using herbivore_run) with these parameters.
        sim_results, sol_full = herbivore_run(
            cell, mu_val, mu_predation, epsilon_val, symmetrical_competition, mean_m_alpha;
            include_predators=include_predators, time_end=500.0, plot=plot,
            NPP=NPP, artificial_pi=artificial_pi, alpha=alpha,
            ignore_inf_error = true
        )
        if mu_val == 0.0
            println(sol_full[1:end, end])
        end

        if sim_results === nothing
            continue
        end
        sr = sim_results.survival_rate[1]
        # Use our stability check on the full ODE simulation:
        stable_full, max_rel = is_stable(sol_full)
        stable_by_nico = is_stable_by_nico(sol_full)
        push!(results_table, (
            mu_val, sol.f_converged, sol.x_converged, stable_full, max_rel, stable_by_nico, sr,
            mu_predation, epsilon_val 
            )
        )
    end
    if any(results_table.stable_by_nico)
        @info "Stable by Nico at mu = $(results_table.mu[findfirst(results_table.stable_by_nico)])"
    end
    return results_table
end

# Example call:
AAAA = explore_stability(1; 
    NPP=Float64(npp_DA_relative_to_1000[idx[1][1], idx[1][2]]),
    mu_range=range(0.0, stop=1.0, length=50),
    symmetrical_competition=true,
    mean_m_alpha=0.1, epsilon_val=0.0,
    mu_predation=0.0, artificial_pi=false,
    alpha=0.25,
    plot = false, # Beware plotting will be very slow for more than a few iterations (i.e. more than a few mu values)
    include_predators = true
)

h_run = herbivore_run(
    1, 
    0.05069, 0.01448, 0.0888886, true, 0.022613;
    include_predators=true, time_end=500.0, plot=true,
    NPP=Float64(npp_DA_relative_to_1000[idx[1][1], idx[1][2]]), artificial_pi=false,
    alpha=0.25,
    ignore_inf_error = true,
    hollingII = true, h = 0.1
);

println("Stability exploration results:")
println(results_table)
