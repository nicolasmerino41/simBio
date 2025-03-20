using LinearAlgebra, Random, Statistics, NLsolve, ForwardDiff, DataFrames, DifferentialEquations

# --- Helper functions for generating random matrices with given connectance ---
function random_matrix(rows::Int, cols::Int, value::Float64, connectance::Float64)
    M = zeros(rows, cols)
    for i in 1:rows, j in 1:cols
        if rand() < connectance
            M[i, j] = value  # you could also draw from a distribution if desired
        end
    end
    return M
end

# --- Abstract community parametrisation ---
function abstract_parametrise_community(
    S::Int,            # number of herbivore species
    R::Int;            # number of predator species
    NPP::Float64 = 1000.0,
    mu::Float64 = 0.5,           # baseline interspecific competition strength
    mu_pred::Float64 = 0.01,     # baseline attack rate
    epsilon::Float64 = 1.0,      # trophic efficiency
    c::Float64 = 0.5,            # connectance of competition matrix
    f::Float64 = 0.3,            # connectance of trophic interactions
    M_mean::Float64 = 1.0,       # mean body mass (can set to one for simplicity)
    alpha::Float64 = 0.25        # allometric exponent
)
    # Set seed for reproducibility if desired
    Random.seed!(1234)
    
    ## HERBIVORES
    # Generate abstract herbivore abundances (H0) – e.g., positive random numbers.
    H0 = rand(S) .+ 0.5

    # Mortality rates (m): use a positive random variable around M_mean.
    m = abs.(randn(S) .* 0.1 .+ M_mean)
    
    # Competition matrix: diagonal elements represent intraspecific competition.
    comp_matrix = random_matrix(S, S, mu, c)
    for i in 1:S
        comp_matrix[i, i] = 1.0  # fix intraspecific competition at 1.0
    end

    # Body masses: for simplicity, set them all equal to M_mean.
    bodyMasses = fill(M_mean, S)
    
    # Allometric growth rates:
    norm_sq = sum((bodyMasses[i]^(-alpha) * H0[i])^2 for i in 1:S)
    raw_g = [ (NPP * bodyMasses[i]^(-2*alpha) * H0[i]) / norm_sq for i in 1:S ]
    beta = raw_g ./ m .- 1
    g = [ raw_g[i] * (beta[i] / (1 + beta[i])) for i in 1:S ]
    
    ## PREDATORS
    # Create trophic (attack) matrix of dimensions R×S.
    trophic_matrix = random_matrix(R, S, mu_pred, f)
    
    # Predator mortality rates: fixed at 0.1.
    m_alpha = fill(0.1, R)
    
    # Predator self-regulation / interaction matrix:
    A = -10.0 * Matrix{Float64}(I, R, R) + random_matrix(R, R, mu_pred, f)
    A_inv = inv(A)
    
    # Compute nondimensional attack rates for herbivores: A_star (S×R).
    d = m ./ H0
    A_star = zeros(S, R)
    for i in 1:S, α in 1:R
        A_star[i, α] = trophic_matrix[α, i] / d[i]
    end

    # Compute predator-mediated competition (C) and release (G) terms.
    C = zeros(S, S)
    G = zeros(S)
    for i in 1:S
        for j in 1:S
            val = 0.0
            for α in 1:R, β in 1:R
                val += A_star[i, α] * A_inv[α, β] * trophic_matrix[β, j]
            end
            C[i, j] = val
        end
        valG = 0.0
        for α in 1:R, β in 1:R
            valG += A_star[i, α] * A_inv[α, β] * m_alpha[β]
        end
        G[i] = valG
    end

    # Update effective observed herbivore abundance incorporating competition and predation.
    H0_eff = zeros(S)
    for i in 1:S
        comp_sum = sum(i != j ? comp_matrix[i, j] * H0[j] : 0.0 for j in 1:S)
        pred_contrib = sum(A_star[i, :] .* ones(R))
        H0_eff[i] = H0[i] + comp_sum + pred_contrib
    end

    # Return the assembled abstract community parameters.
    return (
        S = S,
        R = R,
        H0 = H0,           # base herbivore abundances
        m = m,             # herbivore mortalities
        comp_matrix = comp_matrix,
        bodyMasses = bodyMasses,
        raw_g = raw_g,     # raw growth rates
        beta = beta,       # niche parameters
        g = g,             # adjusted growth rates
        trophic_matrix = trophic_matrix,  # predator attack matrix (R×S)
        m_alpha = m_alpha, # predator mortalities
        A = A,             # predator self-regulation matrix (R×R)
        A_star = A_star,   # nondimensional attack rates (S×R)
        C = C,             # predator-mediated competition coefficients
        G = G,             # predator release term
        H0_eff = H0_eff    # effective observed herbivore abundance
    )
end

# --- Non-abstract equilibrium system (for herbivores and predators) ---
function abstract_equilibrium_system!(F, u, params)
    # Unpack parameters in the same order as the non_abstract pipeline:
    # S, R, H0_eff, m, g, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha
    S, R, H0_eff, m, g, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha = params
    H = u[1:S]
    P = u[S+1:S+R]
    
    # Herbivore equilibrium equations:
    for i in 1:S
        comp_term = sum(comp_matrix[i, j] * H[j] for j in 1:S)
        linear_pred_term = (R > 0 ? sum(A_star[i, α] * P[α] for α in 1:R) : 0.0)
        F[i] = H[i] + comp_term + linear_pred_term - H0_eff[i]
    end
    
    # Predator equilibrium equations:
    for α in 1:R
        attack_sum = sum(A_pred[α, i] * H[i] for i in 1:S)
        interact_sum = sum(A[α, β] * P[β] for β in 1:R)
        F[S + α] = attack_sum - m_alpha[α] - interact_sum
    end
end

# Compute the Jacobian numerically using ForwardDiff.
function abstract_compute_jacobian(f!, u, params)
    return ForwardDiff.jacobian(x -> begin
        F = similar(x)
        f!(F, x, params)
        F
    end, u)
end

# Check stability using time series from the ODE solution.
function abstract_is_stable(sol; window_fraction=0.1, tol=0.01)
    t_final = sol.t[end]
    window_start = t_final * (1 - window_fraction)
    idx_window = findall(t -> t >= window_start, sol.t)
    S = length(sol.u[1])  # number of state variables (herbivores only if predators are not present)
    H_mat = hcat(sol.u...)[1:S, idx_window]
    rel_diffs = [ maximum(abs.(diff(H_mat[i, :]) ./ H_mat[i, 1:end-1])) for i in 1:S ]
    return maximum(rel_diffs) < tol, maximum(rel_diffs)
end

function abstract_is_stable_by_nico(sol)
    final_abundances = sol.u[end]
    return all(x -> x ≤ 100.0, final_abundances)
end

# --- Stability exploration for abstract communities ---
function abstract_explore_stability(;
    S::Int = 2,                         # number of herbivore species
    R::Int = 0,                         # number of predator species
    NPP::Float64 = 1000.0,
    mu_range = range(0.0, stop=1.0, length=20),
    M_mean::Float64 = 0.1,
    c::Float64 = 0.4,                   # connectance for competition matrix
    alpha::Float64 = 0.25,
    f::Float64 = 0.2,                    # trophic connectance
    mu_pred::Float64 = 0.0,             # baseline attack rate
    epsilon::Float64 = 1.0
)
    results_table = DataFrame(mu = Float64[], f_converged = Bool[], x_converged = Bool[],
                              stable = Bool[], max_rel_change = Float64[],
                              stable_by_nico = Bool[], H_eq = Vector{Float64}[],
                              P_eq = Vector{Float64}[],
                              H_abund = Vector{Float64}[], P_abund = Vector{Float64}[])
    
    # For herbivore-only runs (R == 0) we use:
    f_val = R > 0 ? f : 0.0    # no trophic interactions if no predators
    mu_pred = R > 0 ? mu_pred : 0.0
    epsilon_val = R > 0 ? epsilon : 0.0
    
    for mu_val in mu_range
        # Generate an abstract community.
        results = abstract_parametrise_community(S, R;
                        NPP = NPP, mu = mu_val, mu_pred = mu_pred, epsilon = epsilon_val,
                        c = c, f = f_val, M_mean = M_mean, alpha = alpha)
        
        # Destructure the returned NamedTuple.
        S2       = results.S
        R2       = results.R
        H_i0     = results.H0_eff
        m_i      = results.m
        g_i      = results.g
        beta     = results.beta
        M_mod    = results.comp_matrix    # Modified competition matrix (S×S)
        # For predators: use the trophic_matrix directly as A_pred.
        A_pred   = (R2 > 0 ? results.trophic_matrix : zeros(0, S2))
        A_star   = results.A_star
        A        = results.A
        m_alpha  = results.m_alpha
        # Set predator overpopulation thresholds.
        P0       = m_alpha
        
        # Build the parameter tuple (same order as in abstract_equilibrium_system!).
        params = (S2, R2, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha)
        # Use H_i0 and a guess for predator abundances.
        x0 = float.(vcat(H_i0, (R2 > 0 ? H_i0[1:R2].*0.1 : Float64[])))
        
        sol = nlsolve((F, x) -> abstract_equilibrium_system!(F, x, params), x0)
        if !sol.f_converged
            @warn "nlsolve did not converge for mu = $mu_val"
            continue
        end
        
        H_eq = sol.zero[1:S2]
        P_eq = (R2 > 0 ? sol.zero[S2+1:end] : [])
        
        # Compute the Jacobian eigenvalues.
        state_eq = (R2 > 0 ? sol.zero : H_eq)
        J = abstract_compute_jacobian(abstract_equilibrium_system!, state_eq, params)
        max_real = maximum(real.(eigvals(J)))
        stable_flag = max_real < 0.0
        
        # Run the full simulation.
        sim_results, sol_full = abstract_herbivore_run(
            S, R, mu_val, mu_pred, epsilon_val, M_mean;
            time_end = 500.0, plot = false, NPP = NPP, alpha = alpha,
            do_you_want_sol = true, include_predators = (R2 > 0),
            H_init = H_eq, P_init = P_eq
        )
        if sim_results === nothing
            continue
        end
        
        stable_full, max_rel = abstract_is_stable(sol_full)
        stable_by_nico_flag = abstract_is_stable_by_nico(sol_full)
        
        push!(results_table, (
            mu_val, sol.f_converged, sol.x_converged, stable_full, max_rel,
            stable_by_nico_flag, H_eq, (R2 > 0 ? P_eq : []), sol_full[1:S2, end], (R2 > 0 ? sol_full[S2+1:end, end] : [])
            )
        )
    end
    if any(results_table.stable_by_nico)
        @info "Stable by Nico at mu = $(results_table.mu[findfirst(results_table.stable_by_nico)])"
    end
    return results_table
end

# --- Dynamics simulation for the abstract community ---
function abstract_herbivore_run(
    S::Int, R::Int, mu_val, mu_pred_val, eps_val, M_mean;
    time_end = 500.0, 
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators = false,
    plot = false,
    NPP = 1000.0,
    alpha = 0.25,
    H_init = nothing,
    P_init = nothing,
    hollingII = false,
    h = 0.1
)
    if !include_predators
        R = 0
        f_val = 0.0
    else
        f_val = 0.3
    end
    
    # Generate abstract community parameters.
    results = abstract_parametrise_community(S, R;
                    NPP = NPP, mu = mu_val, mu_pred = mu_pred_val, epsilon = eps_val,
                    c = 0.4, f = f_val, M_mean = M_mean, alpha = alpha)
    
    S2         = results.S
    R2         = results.R
    H0_eff     = results.H0_eff
    m          = results.m
    g          = results.g
    beta       = results.beta
    comp_matrix = results.comp_matrix
    
    # Predator-related parameters.
    A_star         = results.A_star
    A_pred         = (R2 > 0 ? results.trophic_matrix : zeros(0, S2))
    m_alpha        = results.m_alpha
    A              = results.A
    P0             = m_alpha
    
    if H_init === nothing
        H_init = copy(H0_eff)
    end
    if R2 > 0
        if P_init === nothing
            P_init = ones(R2) .* (H_init[1:min(R2, length(H_init))] ./ 10.0)
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    
    # Package parameters for the ODE dynamics.
    if R2 > 0
        if !hollingII
            params = (S2, R2, H0_eff, m, g, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha)
        else
            params = (S2, R2, H0_eff, m, g, beta, comp_matrix, A_star, A_pred, P0, A, m_alpha, h)
        end
    else
        if !hollingII
            params = (S2, R2, H0_eff, m, g, beta, comp_matrix, zeros(S2, 0), zeros(0, S2), zeros(0), zeros(0, 0), zeros(0))
        else
            params = (S2, R2, H0_eff, m, g, beta, comp_matrix, zeros(S2, 0), zeros(0, S2), zeros(0), zeros(0, 0), zeros(0), h)
        end
    end
    
    prob = ODEProblem(new_dynamics!, u0, (0.0, time_end), params)
    sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
    
    if plot    
        fig = Figure(; size = (600, 500))
        ax = Axis(
            fig[1, 1],
            # xlabel="Time",
            # ylabel="Biomass",
            # title="Dynamics for cell $cell",
            # yscale = log ? log10 : identity
            )
        times_combined = sol.t
        
        # @info "1"
        # Plot herbivore dynamics (indices 1:S2) as blue solid lines.
        for i in 1:S2
            lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
        end
        # @info "2"
        # Plot predator dynamics (indices S2+1:S2+R2) as red dashed lines.
        for α in 1:R2
            lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
        end
        # @info "3"
    
        display(fig)
    end
    
    if R2 > 0
        H_end = sol.u[end][1:S2]
        P_end = sol.u[end][S2+1:end]
    else
        H_end = sol.u[end]
        P_end = zeros(0)
    end
    H_end = map(x -> x < 1e-6 ? 0.0 : x, H_end)
    P_end = map(x -> x < 1e-6 ? 0.0 : x, P_end)
    
    survived_herb = count(x -> x > 1e-6, H_end)
    survived_pred = count(x -> x > 1e-6, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species
    
    single_run_results = DataFrame(
        survival_rate = survival_rate,
        H_biomass = sum(H_end),
        P_biomass = sum(P_end),
        total_species = total_species,
        survived_herbivores = survived_herb,
        survived_predators = survived_pred
    )
    
    if do_you_want_params && do_you_want_sol
        return single_run_results, params, sol
    elseif do_you_want_params || do_you_want_sol
        return single_run_results, do_you_want_params ? params : sol
    else
        return single_run_results
    end
end

# Example usage:
AAAA = abstract_explore_stability(;
    S = 5, R = 2,
    epsilon = 1.0, f = 0.2, M_mean = 0.1, alpha = 0.25, mu_pred = 0.0
)

abstract_herbivore_run(
    5, 2, 0.0, 0.01, 0.5, 0.1;
    time_end = 500.0, 
    do_you_want_params = false,
    do_you_want_sol = false,
    include_predators =  true,
    plot = false,
    NPP = 1000.0,
    alpha = 0.25,
    H_init = AAAA[1, :H_eq],
    P_init =fill(0.0, 2)
)

println("P_eq for first entry:")
println(AAAA[1, :P_eq])
