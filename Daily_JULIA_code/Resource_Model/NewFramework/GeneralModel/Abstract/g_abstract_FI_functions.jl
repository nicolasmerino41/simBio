function g_abstract_parametrise_the_community(
    S::Int, O::Int, R::Int;
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    conn = 0.2,                # target connectance (probability an interaction exists)
    cell_abundance_h::Vector{Float64} = Float64[],  # default: not provided → ones
    cell_abundance_p::Vector{Float64} = Float64[],  # default: not provided → ones
    delta_nu = 0.05,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 1.0,
    randomise_attack_rates::Bool = false,
    attack_rate_sd::Float64 = 0.2
)
    # Convert parameters to Float64.
    mu             = Float64(mu)
    epsilon_val  = Float64(epsilon_val)
    mean_m_alpha = Float64(mean_m_alpha)
    delta_nu     = Float64(delta_nu)
    d_alpha      = Float64(d_alpha)
    d_i          = Float64(d_i)
    r_omni_proportion = Float64(r_omni_proportion)
    conn         = Float64(conn)

    # Total number in the herbivore compartment.
    S_total = S + O

    # Define synthetic species names.
    herbivore_list = vcat(["Herbivore $i" for i in 1:S],
                          ["Omnivore $i" for i in 1:O])
    predator_list  = ["Predator $i" for i in 1:R]
    sp_nm = vcat(herbivore_list, predator_list)

    # Process equilibrium abundances:
    H_eq = (S_total > 0) ? (length(cell_abundance_h) == S_total ? copy(cell_abundance_h) : ones(S_total)) : Float64[]
    P_eq = (R > 0) ? (length(cell_abundance_p) == R ? copy(cell_abundance_p) : ones(R)) : Float64[]
    H_star = copy(H_eq)

    # Intrinsic growth rates for herbivores:
    r_vec = vcat(ones(S), ones(O) * r_omni_proportion)  # pure herbivores grow at rate 1, omnivores at r_omni_proportion

    # For predators, assign mortality and conversion efficiency:
    if R > 0
        m_alpha = fill(mean_m_alpha, R)
        # Conversion efficiency scaled by (d_i/d_alpha)
        epsilon_vec = fill(epsilon_val, R) * (d_i / d_alpha)
    else
        m_alpha = Float64[]
        epsilon_vec = Float64[]
    end

    # --- Build the merged interaction matrix A ---
    # A is of size n_total x n_total, where n_total = S_total + R.
    n_total = S_total + R
    A = zeros(Float64, n_total, n_total)
    # We loop over all pairs (i,j). For each potential interaction, we decide (with probability conn)
    # whether an interaction exists. If yes, then we assign:
    #   - The effect on the target (row i) is +1.0 (representing a loss in biomass when attacked).
    #   - The effect on the consumer (row j) is +epsilon_val (a benefit from consumption).
    #
    # Here we distinguish by the types of the species:
    #   - If the target is herbivore (i ≤ S_total) and the consumer is predator (j > S_total), or
    #   - If the target is herbivore and the consumer is omnivore (j ≤ S_total and j > S),
    #   - If both are predators, or if a herbivore (even a pure herbivore) consumes a predator,
    # we use the same rule.
    for i in 1:n_total
        for j in 1:n_total
            if rand() < conn  # an interaction exists
                # Determine species types.
                target_is_herb = (i ≤ S_total)
                consumer_is_pred = (j > S_total)
                consumer_is_omnivore = (j ≤ S_total) && (j > S)  # indices S+1..S_total are omnivores
                consumer_is_herb = (j ≤ S_total) && (j ≤ S)
                # We assign values only if the entry is not already set (to avoid double counting).
                if target_is_herb && consumer_is_pred && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = 1.0         # Herbivore loses biomass.
                    A[j,i] = epsilon_val   # Predator gains biomass.
                elseif target_is_herb && consumer_is_omnivore && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0         # Herbivore loses biomass.
                    A[j,i] = -epsilon_val   # Omnivore gains biomass.
                elseif (!target_is_herb) && consumer_is_pred && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0
                    A[j,i] = epsilon_val
                elseif (!target_is_herb) && (consumer_is_herb || consumer_is_omnivore) && A[i,j] == 0.0 && A[j,i] == 0.0
                    A[i,j] = -1.0
                    A[j,i] = -epsilon_val
                end
            end
        end
    end
    # println("Conn is $(conn) but in reality is $(sum(A .!= 0.0) / (n_total^2))")
    # Scale each row of A:
    for i in 1:S_total
        A[i, :] .= A[i, :] ./ d_i
    end
    for i in S_total+1:n_total
        A[i, :] .= A[i, :] ./ d_alpha
    end

    # --- Build the competition matrix C (for herbivores only) ---
    C = fill(mu, S_total, S_total)
    for i in 1:S_total
        C[i, i] = 0.0
    end
    C .= C ./ d_i

    # --- Compute the effective predation scaling factor ν ---
    # Let X = [H_star; P_eq]
    X = vcat(H_star, P_eq)
    nu_candidates = Float64[]
    # For each predator (indices S_total+alpha), compute:
    for alpha in 1:R
        idx_pred = S_total + alpha
        gain = sum(A[idx_pred, :] .* X)
        if gain > 0
            push!(nu_candidates, (P_eq[alpha] + m_alpha[alpha]) / gain)
        end
    end
    nu = isempty(nu_candidates) ? 0.0 : maximum(nu_candidates)
    nu *= (1.0 + delta_nu)  # apply safety margin

    # --- Embed ν into A (i.e. multiply all entries by ν) ---
    A .= nu .* A

    # --- Optionally randomise attack rates ---
    if randomise_attack_rates
        # Use a truncated Normal distribution with mean 1 and sd attack_rate_sd.
        multiplier_dist = Truncated(Normal(1.0, attack_rate_sd), 0, Inf)
        for i in axes(A,1), j in axes(A,2)
            A[i,j] *= rand(multiplier_dist)
        end
    end

    # --- Recalculate carrying capacities ---
    # For herbivores:
    K_i = zeros(Float64, S_total)
    total_H = sum(H_star)
    for i in 1:S_total
        # Compute trophic effect from all interactions.
        trophic_effect = sum(A[i, :] .* X)
        comp_effect = sum(C[i, :] .* H_star)
        K_i[i] = (1 - mu) * H_star[i] + trophic_effect + comp_effect
    end
    # For predators:
    K_alpha = zeros(Float64, R)
    for alpha in 1:R
        idx_pred = S_total + alpha
        trophic_effect_pred = sum(A[idx_pred, :] .* X)
        K_alpha[alpha] = trophic_effect_pred - P_eq[alpha]
    end

    # (Optional) Compute derived ratios (not used here but could be returned)
    new_di = r_vec ./ K_i
    new_da = m_alpha ./ K_alpha

    return (
        S = S_total, R = R,
        H_eq = H_eq, P_eq = P_eq,
        r_i = r_vec, K_i = K_i,
        mu = mu, nu = nu,
        A_matrix = A,
        C_matrix = C,
        epsilon = epsilon_vec, m_alpha = m_alpha, K_alpha = K_alpha,
        herbivore_list = herbivore_list, predator_list = predator_list,
        species_names = sp_nm,
        H_star = H_star, P_star = P_eq
    )
end

"""
    reverse_abstract_parametrise_the_community(S, O, R; mu, epsilon_val, mean_m_alpha, nu, conn, 
        d_alpha, d_i, r_omni_proportion, randomise_attack_rates, attack_rate_sd, max_sim_time)

Constructs an abstract community (without real data) using the new framework. Instead of specifying equilibrium abundances,
you provide a desired effective predation scaling factor `nu` (which will be “embedded” into A_matrix) and we then determine the equilibrium abundances by simulating the dynamics until near‐steady state. In our new formulation only two matrices are used:
  - A_matrix: the merged interaction matrix (with negative entries on targets and positive on consumers)
  - C_matrix: the herbivore competition matrix.

Arguments:
- S::Int: number of pure herbivores.
- O::Int: number of omnivores.
- R::Int: number of predators.
- mu: competition coefficient among herbivores.
- epsilon_val: baseline conversion efficiency.
- mean_m_alpha: mean predator mortality.
- nu: effective predation scaling factor to embed in the interactions.
- conn: target probability that any given potential interaction exists.
- d_alpha, d_i: self‐regulation (scaling) parameters for predators and herbivores, respectively.
- r_omni_proportion: intrinsic growth rate for omnivores (relative to 1 for pure herbivores).
- randomise_attack_rates (Bool): if true, each entry in A_matrix is multiplied by a random factor.
- attack_rate_sd: standard deviation for the random multiplier (the multiplier is drawn from a truncated Normal with mean 1).
- max_sim_time: maximum simulation time for the ODE solver to (approximately) reach equilibrium.

Returns a NamedTuple containing:
  - S: total number of herbivores (S + O)
  - R: number of predators
  - H_eq: the equilibrium herbivore abundances (obtained by simulation)
  - P_eq: the equilibrium predator abundances (obtained by simulation)
  - r_i: intrinsic herbivore growth rates (vector; omnivores get r_omni_proportion)
  - K_i: carrying capacities computed from the herbivore equilibrium conditions
  - mu: competition coefficient (as input)
  - nu: effective predation rate (as provided, after safety margin)
  - A_matrix: the merged interaction matrix (with ν embedded and optionally randomized)
  - C_matrix: the competition matrix among herbivores
  - epsilon: conversion efficiencies (for predators)
  - m_alpha: predator mortality rates
  - K_alpha: carrying capacities for predators computed from equilibrium conditions
  - herbivore_list: synthetic names for herbivores/omnivores
  - predator_list: synthetic names for predators
  - species_names: the full list of synthetic species names
  - H_star: the equilibrium herbivore abundances (same as H_eq)
  - P_star: the equilibrium predator abundances (same as P_eq)
"""
function reverse_abstract_parametrise_the_community(
    S::Int, O::Int, R::Int;
    mu = 0.5,
    epsilon_val = 1.0,
    mean_m_alpha = 0.1,
    nu = 0.1,
    conn = 0.2,
    d_alpha = 1.0,
    d_i = 1.0,
    r_omni_proportion = 1.0,
    randomise_attack_rates::Bool = false,
    attack_rate_sd::Float64 = 0.2,
    max_sim_time = 1000.0
)
    # Convert parameters to Float64.
    mu              = Float64(mu)
    epsilon_val     = Float64(epsilon_val)
    mean_m_alpha    = Float64(mean_m_alpha)
    nu              = Float64(nu)
    conn            = Float64(conn)
    d_alpha         = Float64(d_alpha)
    d_i             = Float64(d_i)
    r_omni_proportion = Float64(r_omni_proportion)

    # Total number in the herbivore compartment.
    S_total = S + O

    # Generate synthetic species names.
    herbivore_list = vcat(["Herbivore $i" for i in 1:S],
                          ["Omnivore $i" for i in 1:O])
    predator_list  = ["Predator $i" for i in 1:R]
    sp_nm = vcat(herbivore_list, predator_list)

    # Set baseline initial abundances (guess): ones.
    H_init = ones(S_total)
    P_init = ones(R)
    H_star_guess = copy(H_init)

    # Define intrinsic growth rates for herbivores.
    r_i = vcat(ones(S), ones(O) * r_omni_proportion)

    # For predators, assign mortality and conversion efficiency.
    m_alpha = R > 0 ? fill(mean_m_alpha, R) : Float64[]
    epsilon_vec = R > 0 ? fill(epsilon_val, R) * (d_i / d_alpha) : Float64[]

    # --- Build merged interaction matrix A_matrix ---
    n_total = S_total + R
    A = zeros(Float64, n_total, n_total)
    # We now decide independently for each possible (ordered) pair (i,j) whether an interaction exists.
    # The rule is as follows (based on our new framework):
    #   For an interaction where the target is a herbivore (i ≤ S_total) and the consumer is
    #     - a predator (j > S_total) or an omnivore (S < j ≤ S_total): 
    #         The target loses biomass: assign -1.
    #         The consumer gains biomass: assign +ε_val.
    #   For interactions among predators or herbivore–herbivore competition, we use the same rule.
    for i in 1:n_total, j in 1:n_total
        if rand() < conn
            target_is_herb = (i ≤ S_total)
            consumer_is_pred = (j > S_total)
            consumer_is_omnivore = (j ≤ S_total) && (j > S)  # indices S+1..S_total
            consumer_is_herb = (j ≤ S_total) && (j ≤ S)
            # We assign only if not previously set (to avoid double counting):
            if target_is_herb && (consumer_is_pred || consumer_is_omnivore)
                A[i,j] = -1.0     # target loses biomass
                A[j,i] = epsilon_val  # consumer gains biomass
            else
                # For interactions among predators:
                if !target_is_herb && consumer_is_pred
                    A[i,j] = -1.0
                    A[j,i] = epsilon_val
                else
                    # For herbivore-herbivore competitive interactions:
                    A[i,j] = -1.0
                    A[j,i] = 0.0
                end
            end
        end
    end

    # Scale rows of A_matrix: herbivores by d_i, predators by d_alpha.
    for i in 1:S_total
        A[i, :] ./= d_i
    end
    for i in S_total+1:n_total
        A[i, :] ./= d_alpha
    end

    # --- Build the competition matrix C_matrix for herbivores ---
    C = fill(mu, S_total, S_total)
    for i in 1:S_total
        C[i,i] = 0.0
    end
    C ./= d_i

    # --- Embed nu into A_matrix ---
    # Multiply the entire merged interaction matrix by the user‐supplied nu.
    A .= nu .* A

    # --- Optionally randomise the attack rates in A_matrix ---
    if randomise_attack_rates
        multiplier_dist = truncated(Normal(1.0, attack_rate_sd), 0, Inf)
        for i in axes(A,1), j in axes(A,2)
            A[i,j] *= rand(multiplier_dist)
        end
    end

    # --- Obtain equilibrium abundances by simulating the ODE --- 
    # We use the following ODE dynamics:
    #   For herbivores:
    #     dH/dt = H .* r_i .* ( 1 - [ (1-mu)H + C*H + A[1:S_total, :]*X ] ./ K_i )
    #   For predators:
    #     dP/dt = P .* m_alpha .* ( ([A[S_total+1:end, :]*X] - P) ./ K_alpha - 1 )
    # Here we do not supply K_i and K_alpha a priori; instead we choose initial guesses and later compute them based on the equilibrium.
    #
    # We define a temporary ODE function (general_dynamics!) that expects K_i and K_alpha as parameters.
    # In our scheme, we will set K_i and K_alpha according to the equilibrium conditions:
    #   For herbivores: K_i = (1-mu) * H_i* + (A[1:S_total, :]*X)_i + (C*H)_i.
    #   For predators:   K_alpha = (A[S_total+1:end, :]*X)_α - P_α*.
    # Our strategy is to run a simulation for a long time (max_sim_time) with initial guess X0 = ones(n_total)
    # and then use the final state X* as the equilibrium abundances.
    X0 = ones(n_total)
    # Here we create a temporary K vector just to allow simulation; we set them to the values computed from X0.
    K_i_temp = [(1 - mu) * 1.0 + sum(A[i, :] .* X0) + sum(C[i, :] .* ones(S_total)) for i in 1:S_total]
    K_alpha_temp = [sum(A[S_total+α, :] .* X0) - 1.0 for α in 1:R]
    p_dyn = (S_total, R, K_i_temp, r_i, mu, nu, A, C, m_alpha, K_alpha_temp)

    # Run the simulation for max_sim_time to approximate the equilibrium.
    prob = ODEProblem(general_dynamics!, X0, (0.0, max_sim_time), p_dyn)
    sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    X_star = sol.u[end]
    H_eq_found = X_star[1:S_total]
    P_eq_found = X_star[S_total+1:end]

    # --- Recalculate carrying capacities using the found equilibrium abundances ---
    X_star_full = X_star  # = [H_eq_found; P_eq_found]
    K_i = zeros(S_total)
    for i in 1:S_total
        K_i[i] = (1 - mu) * H_eq_found[i] + sum(A[i, :] .* X_star_full) + sum(C[i, :] .* H_eq_found)
    end
    K_alpha = zeros(R)
    for α in 1:R
        idx_pred = S_total + α
        K_alpha[α] = sum(A[idx_pred, :] .* X_star_full) - P_eq_found[α]
    end
    if any(sol[:, end] .< EXTINCTION_THRESHOLD)
        println("Some species went extinct during the simulation.")
    end

    return (
        S = S_total,
        R = R,
        H_eq = H_eq_found,
        P_eq = P_eq_found,
        r_i = r_i,
        K_i = K_i,
        mu = mu,
        nu = nu,
        A_matrix = A,
        C_matrix = C,
        epsilon = epsilon_vec,
        m_alpha = m_alpha,
        K_alpha = K_alpha,
        herbivore_list = herbivore_list,
        predator_list = predator_list,
        species_names = sp_nm,
        H_star = H_eq_found,
        P_star = P_eq_found
    )
end

for i in 0.0:0.1:1.0
    A = reverse_abstract_parametrise_the_community(
        10, 5, 5;
        mu = 0.1,
        epsilon_val = 0.1,
        mean_m_alpha = 0.05,
        nu = i
    );
    println("For nu $i the eq abunda are $(vcat(A.H_star, A.P_star))")
end

A.H_star
A.P_star

function reverse_abstract_run(
    S::Int, O::Int, R::Int,
    mu_val,
    eps_val,
    mean_m_alpha,
    nu_val;
    d_alpha = 1.0, d_i = 1.0,
    time_end = 500.0,
    do_you_want_params = false,
    do_you_want_sol = false,
    plot = false,
    H_init = nothing,
    P_init = nothing,
    ignore_inf_error = false,
    log = false,
    initial_deviation = 0.0,
    extinguishing_species = 0,
    nu_omni_proportion = 1.0,
    nu_b_proportion = 1.0,
    r_omni_proportion = 0.01,
    force_nu_to = nothing,
    use_cb = false,
    perturbation_halfway = false,
    species_to_perturb = 0,
    removal_fraction = 0.0,
    conn = 0.2,  # target connectance for synthetic interaction matrices
)
    AAAA = DataFrame()

    # Build the synthetic community using the new abstract function.
    params_setup = reverse_abstract_parametrise_the_community(
        S, O, R;
        mu = mu_val,
        epsilon_val = eps_val,
        mean_m_alpha = mean_m_alpha,
        nu = nu_val,
        conn = conn,
        d_alpha = d_alpha,
        d_i = d_i,
        r_omni_proportion = r_omni_proportion
    )

    if isnothing(params_setup)
        @error "Error: params_setup is nothing"
        return nothing
    end

    # Destructure the returned parameters.
    S_total    = params_setup.S       # Total number of herbivores+omnivores
    R          = params_setup.R
    H_eq       = params_setup.H_eq      # Baseline herbivore (and omnivore) abundances
    P_eq       = params_setup.P_eq      # Baseline predator abundances
    r_i        = params_setup.r_i       # Herbivore intrinsic growth rates
    K_i        = params_setup.K_i       # Herbivore carrying capacities
    nu_val     = params_setup.nu        # Effective predation rate
    A_matrix   = params_setup.A_matrix   # Trophic interaction matrix
    C_matrix   = params_setup.C_matrix   # Competition interaction matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    herbivore_list = params_setup.herbivore_list
    predator_list  = params_setup.predator_list

    # println("Herbivore list length: ", length(herbivore_list))
    
    # Initial conditions: if not provided, use the baseline equilibrium.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        if extinguishing_species == 0
            u0 = abs.(vcat(H_init, P_init) .+ randn(S_total+R) .* initial_deviation)
        else
            if extinguishing_species > S_total+R
                @error "Error: extinguishing_species > total species"
            end
            u0 = vcat(H_init, P_init) .* [i == extinguishing_species ? 0.0 : 1.0 for i in 1:(S_total+R)]
        end
    else
        u0 = H_init
    end

    # Build the parameter tuple for the ODE.
    if !isnothing(force_nu_to)
        nu_val = force_nu_to
    end
    params = (S+O, R, K_i, r_i, mu_val, nu_val, A_matrix, C_matrix, m_alpha, K_alpha)

    # Solve the ODE using omnivore_dynamics!
    prob = ODEProblem(general_dynamics!, u0, (0.0, time_end), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if use_cb
            solve(prob, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
        else
            solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)
        end
    end

    if !ignore_inf_error
        if sol.t[end] < time_end || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @error "Error: solution did not finish properly"
            return nothing
        end
    end

    # --- PERTURBATION HALFWAY (optional) ---
    if perturbation_halfway
        if species_to_perturb > 0 && species_to_perturb <= S_total+R
            u0_perturbed = sol[:, end]
            u0_perturbed[species_to_perturb] -= u0_perturbed[species_to_perturb] * removal_fraction
            prob2 = ODEProblem(general_dynamics!, u0_perturbed, (time_end, time_end+time_end), params)
            if use_cb
                new_sol = solve(prob2, Tsit5(); callback = cb_no_trigger, abstol = 1e-8, reltol = 1e-6)
            else
                new_sol = solve(prob2, Tsit5(); abstol = 1e-8, reltol = 1e-6)
            end
        elseif species_to_perturb > S_total+R
            error("Error: species_to_perturb > total species")
        else
            error("Error: perturbation_halfway activated but species_to_perturb is not specified")
        end
    end

    # --- Plotting using Makie ---
    if plot    
        fig = Figure(; size = (900, 600))
        ax = Axis(fig[1, 1],
            xlabel = "Time",
            ylabel = "Biomass",
            title = "Dynamics for abstract community",
            yscale = log ? log10 : identity)
        times = sol.t
        for i in 1:S_total
            # For simplicity, use blue for herbivores and green for omnivores.
            color = (herbivore_list[i] in herbivore_list && !occursin("Omnivore", herbivore_list[i])) ? :blue : :green
            lines!(ax, times, sol[i, :], label = "$(herbivore_list[i])", color = color)
        end
        for alpha in 1:R
            lines!(ax, times, sol[S_total+alpha, :], label = "$(predator_list[alpha])", linestyle = :dash, color = :red)
        end

        if perturbation_halfway
            ax2 = Axis(fig[1, 2],
                xlabel = "Time",
                ylabel = "Biomass",
                title = "After perturbation",
                yscale = log ? log10 : identity)
            times2 = new_sol.t
            for i in 1:S_total
                color = (herbivore_list[i] in herbivore_list && !occursin("Omnivore", herbivore_list[i])) ? :blue : :green
                lines!(ax2, times2, new_sol[i, :], label = "$(herbivore_list[i])", color = color)
            end
            for alpha in 1:R
                lines!(ax2, times2, new_sol[S_total+alpha, :], label = "$(predator_list[alpha])", linestyle = :dash, color = :red)
            end

            # # Plot total biomass before and after perturbation (excluding data between time_end and time_end+20)
            # ax3 = Axis(fig[2, 1],
            #     xlabel = "Time",
            #     ylabel = "Total Biomass",
            #     title = "Total Biomass Dynamics")
            # filtered_old = sol[:, sol.t .< time_end]
            # filtered_new = new_sol[:, new_sol.t .> time_end + 20]
            # times_old = sol.t[sol.t .< time_end]
            # times_new = new_sol.t[new_sol.t .> time_end + 20]
            # combined_times = vcat(times_old, times_new)
            # combined_sol = hcat(filtered_old, filtered_new)
            # total_biomass = vec(sum(combined_sol, dims=1))
            # lines!(ax3, combined_times, total_biomass, label = "Total Biomass", color = :black)
        end
        display(fig)
    end

    # --- Post-Processing: Extract final densities and compute survival rates ---
    H_end = sol[1:S_total, end]
    P_end = sol[S_total+1:S_total+R, end]
    H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, H_end)
    P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, P_end)
    survived_herb = count(x -> x > EXTINCTION_THRESHOLD, H_end)
    survived_pred = count(x -> x > EXTINCTION_THRESHOLD, P_end)
    total_surv = survived_herb + survived_pred
    total_species = S_total + R
    survival_rate = total_surv / total_species
    herbivore_survival_rate = survived_herb / S_total
    predator_survival_rate = (R > 0 ? survived_pred / R : 0.0)
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0 ? NaN : (P_biomass / H_biomass))
    println("Estimated nu = $nu_val")

    single_run_results = DataFrame(
        survival_rate           = survival_rate,
        h_survival              = "$(survived_herb)/$S_total",
        p_survival              = "$(survived_pred)/$R",
        biomass_at_the_end      = biomass_at_the_end,
        H_vector                = [H_end],
        P_vector                = [P_end],
        H_eq                    = [H_eq],
        P_eq                    = [P_eq],
        survived_herbivores     = survived_herb,
        survived_predators      = survived_pred,
        H_biomass               = H_biomass,
        P_biomass               = P_biomass,
        herb_pred_ratio         = ratio,
        total_survivors         = total_surv,
        total_species           = total_species,
        herbivore_survival_rate = herbivore_survival_rate,
        predator_survival_rate  = predator_survival_rate,
        mu                      = mu_val,
        nu                      = nu_val,
        epsilon_val             = eps_val,
        mean_m_alpha            = mean_m_alpha
    )
    
    append!(AAAA, single_run_results)
    
    if perturbation_halfway
        new_H_end = new_sol[1:S_total, end]
        new_P_end = new_sol[S_total+1:S_total+R, end]
        new_H_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_H_end)
        new_P_end = map(x -> x < EXTINCTION_THRESHOLD ? 0.0 : x, new_P_end)
        new_survived_herb = count(x -> x > EXTINCTION_THRESHOLD, new_H_end)
        new_survived_pred = count(x -> x > EXTINCTION_THRESHOLD, new_P_end)
        new_total_surv = new_survived_herb + new_survived_pred
        new_total_species = S_total + R
        new_survival_rate = new_total_surv / new_total_species
        new_herbivore_survival_rate = new_survived_herb / S_total 
        new_predator_survival_rate = (R > 0 ? new_survived_pred / R : 0.0)
        new_H_biomass = sum(new_H_end)
        new_P_biomass = sum(new_P_end)
        new_biomass_at_the_end = new_H_biomass + new_P_biomass
        new_ratio = (new_H_biomass == 0.0 ? NaN : (new_P_biomass / new_H_biomass))

        new_single_run_results = DataFrame(
            survival_rate           = new_survival_rate,
            h_survival              = "$(new_survived_herb)/$S_total",
            p_survival              = "$(new_survived_pred)/$R",
            biomass_at_the_end      = new_biomass_at_the_end,
            H_vector                = [new_H_end],
            P_vector                = [new_P_end],
            H_eq                    = [H_eq],
            P_eq                    = [P_eq],
            survived_herbivores     = new_survived_herb,    
            survived_predators      = new_survived_pred,
            H_biomass               = new_H_biomass,
            P_biomass               = new_P_biomass,
            herb_pred_ratio         = new_ratio,
            total_survivors         = new_total_surv,
            total_species           = new_total_species,
            herbivore_survival_rate = new_herbivore_survival_rate,
            predator_survival_rate  = new_predator_survival_rate,
            mu                      = mu_val,
            nu                      = nu_val,
            epsilon_val             = eps_val,
            mean_m_alpha            = mean_m_alpha
        )
        
        append!(AAAA, new_single_run_results)
    end

    initial_biomass = sum(u0)
    final_biomass = sum(sol[:, end])
    println("Initial biomass was $initial_biomass and final biomass is $final_biomass")

    if perturbation_halfway
        diff = []
        for i in 2:size(new_sol, 2)
            push!(diff, abs(sum(new_sol[:, i]) - biomass_at_the_end))
        end
        println("The maximum absolute difference is ", maximum(diff), " and it was found at time ", new_sol.t[argmax(diff)])
    end

    if do_you_want_params && do_you_want_sol
        return AAAA, params, sol
    elseif do_you_want_params || do_you_want_sol
        return AAAA, do_you_want_params ? params : sol
    else
        return AAAA
    end

end

begin
    S, O, R = 10, 5, 5
    # cb_no_trigger, cb_trigger = build_callbacks(S+O, R, EXTINCTION_THRESHOLD, T_ext, 1)
    @time A_run = reverse_abstract_run(
        S, O, R,    
        0.8, 1.0, 0.5, 1.0; # mu, epsilon, m_alpha, nu
        d_alpha = 1.0, d_i = 1.0,
        time_end = 500.0,
        do_you_want_params = false,
        do_you_want_sol = false,
        plot = true,
        H_init = nothing,
        P_init = nothing,
        ignore_inf_error = true,
        log = false,
        initial_deviation = 0.0,
        extinguishing_species = 0,
        nu_omni_proportion = 1.0,
        nu_b_proportion = 1.0,
        r_omni_proportion = 1.0,
        force_nu_to = nothing,
        use_cb = false,
        perturbation_halfway = false,
        species_to_perturb = 10,
        removal_fraction = 0.1,
        conn = 0.1  # target connectance for synthetic interaction matrices
    )
end