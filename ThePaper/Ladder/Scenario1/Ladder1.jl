################ INFORMATION #######################
# a) Ladder1.jl is for all the functions
# b) Ladder2.jl is for running the simulations
# c) Ladder3.jl is for post-processing and plotting
#############################################################################
#############################################################################
################## TROPHIC DYNAMICS ODE MODEL ###############################
#############################################################################
#############################################################################
function trophic_ode!(du, u, p, t)
    # Unpack parameters from tuple
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    
    # u is arranged as:
    # u[1:R]          ? resources
    # u[R+1 : R+C]  ? consumers
    
    # Resources dynamics (indices 1:R)
    for i in 1:R
        # Sum over predators of resource i (these are consumers, indices: R+1 to R+C)
        pred_sum = 0.0
        for j in R+1:R+C
            pred_sum += A[i, j] * u[j]
        end
        # Resource i dynamics
        du[i] = u[i] * d_res[i] * ( (r_res[i] / d_res[i]) - u[i] + pred_sum )
    end

    # Consumers dynamics (indices R+1 to R+C)
    for k in 1:C
        i = R + k  # global index for consumer k
        
        # Gains from prey (only A_used[i,j]>0)
        # println("A[i,j] = ", typeof(A))
        # println("epsilon[i,j] = ", typeof(epsilon))
        sum_prey = sum(
            epsilon[i,j] * A[i,j] * u[j]
            for j in 1:R+C if A[i,j] > 0.0;
            init = 0.0
        )

        # Losses to predators (only A_used[i,j]<0)
        sum_pred = sum(
            A[i,j] * u[j]
            for j in 1:R+C if A[i,j] < 0.0;
            init = 0.0
        )
        
        du[i] = u[i] * (m_cons[k] / xi_cons[k]) * ( -xi_cons[k] - u[i] + sum_prey + sum_pred )
    end
end

#############################################################################
#############################################################################
################## CALIBRATE PARAMETRISATION ################################
#############################################################################
#############################################################################
function calibrate_params(
    R_eq::Vector{Float64},
    C_eq::Vector{Float64},
    p_calib;
    xi_threshold = 0.3
)
    R, C, m_cons, d_res, epsilon, A_used = p_calib
    total = R + C
    B_eq = vcat(R_eq, C_eq)

    new_xi_cons = zeros(C)
    new_r_res   = zeros(R)

    # 1) Consumer thresholds
    for k in 1:C
        i = R + k

        prey_sum = sum(
            epsilon[i,j] * A_used[i,j] * B_eq[j]
            for j in 1:total if A_used[i,j] > 0.0;
            init = 0.0
        )

        pred_sum = sum(
            A_used[i,j] * B_eq[j]
            for j in 1:total if A_used[i,j] < 0.0;
            init = 0.0
        )

        xi = -C_eq[k] + prey_sum - pred_sum
        new_xi_cons[k] = xi > 0 ? xi : NaN

        xi = prey_sum - pred_sum - C_eq[k]
        if xi > 0 && (C_eq[k]/xi) < xi_threshold
            new_xi_cons[k] = xi
        else
            new_xi_cons[k] = NaN
        end
    end

    # 2) Resource growth rates
    for i in 1:R
        # Consumer rows are R+1:total; losses where those entries <0
        pred_sum = sum(
            A_used[i,j] * B_eq[j]
            for j in 1:R+C if A_used[i,j] < 0.0;
            init = 0.0
        )
        r = d_res[i] * (R_eq[i] + pred_sum)
        new_r_res[i] = r > 0 ? r : NaN
    end

    return new_xi_cons, new_r_res
end
#############################################################################
#############################################################################
################## JACOBIAN, RESILIENCE AND RESISTANCE ######################
#############################################################################
#############################################################################
# Revised Jacobian computation
function compute_jacobian(B, p)
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    # println("B is ", typeof(B))
    # println("r_res is ", typeof(r_res))
    total = R + C
    J = zeros(total, total)
    # Resources: indices 1:R
    for i in 1:R
        # Diagonal: diffR_i/diffR_i
        J[i, i] = - B[i] * d_res[i]
        # For consumers (indices R+1 to R+C)
        for j in (R+1):total
            J[i, j] = - B[i] * d_res[i] * A[j, i]
        end
        # Off-diagonals with other resources remain 0.
    end
    # Consumers: indices R+1 to total.
    for k in 1:C
        i = R + k  # global index for consumer k
        alpha = m_cons[k] / xi_cons[k]
        # Diagonal (diffC_k/diffC_k)
        J[i, i] = - alpha * B[i]
        # Derivatives with respect to resources (indices 1:R)
        for j in 1:R
            J[i, j] = B[i] * alpha * epsilon[i, j] * A[i, j]
        end
        # Derivatives with respect to other consumers (indices R+1 to total)
        for j in (R+1):total
            if j != i
                J[i, j] = - B[i] * alpha * A[j, i]
            end
        end
    end
    return J
end

# Resilience: negative of the largest real part of the Jacobian eigenvalues.
function compute_resilience(B, p)
    J = compute_jacobian(B, p)
    ev = eigvals(J)
    return maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of the symmetric part of the Jacobian.
function compute_reactivity(B, p)
    J = compute_jacobian(B, p)
    J_sym = (J + J') / 2
    ev_sym = eigvals(J_sym)
    return maximum(real.(ev_sym))
end

#############################################################################
#############################################################################
################## SIMULATE PRESS PERTURBATIONS #############################
#############################################################################
#############################################################################
function simulate_press_perturbation(
    u0, p, tspan, t_perturb, delta;
    solver=Tsit5(),
    plot=false,
    show_warnings=true,
    full_or_simple=true,
    cb = cb
)
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p

    # --- Phase 1: run up to t_perturb ---
    prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], t_perturb), p)
    sol1 = show_warnings ? solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8) :
                           with_logger(logger) do
                               solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8)
                           end
    pre_state = sol1.u[end]
    before_persistence = count(x -> x > EXTINCTION_THRESHOLD, pre_state) / length(pre_state)
    # pre_state[pre_state .> EXTINCTION_THRESHOLD] .= 0.0
    # --- Phase 2: apply press (reduce thresholds by delta) ---
    xi_press = xi_cons .* (1 .- delta)
    p_press  = (R, C, m_cons, xi_press, r_res, d_res, epsilon, A)
    prob2    = ODEProblem(trophic_ode!, pre_state, (t_perturb, tspan[2]), p_press)
    sol2     = solve(prob2, solver; callback = cb, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)

    # --- Return Times ---
    return_times = fill(NaN, n)
    for i in 1:n
        target = new_equil[i]
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target) / (abs(target) + 1e-8) < 0.10
                return_times[i] = t - t_perturb
                break
            end
        end
    end

    # --- Maximum Overshoot & Integrated Recovery Error ---
    overshoot = zeros(n)
    ire       = zeros(n)
    for i in 1:n
        # compute relative errors at each timepoint
        rel_errors = [abs(state[i] - new_equil[i]) / (abs(new_equil[i]) + 1e-8)
                      for state in sol2.u]
        overshoot[i] = maximum(rel_errors)
        ire[i]       = mean(rel_errors)
    end

    # --- Optional plotting ---
    if plot && sol1.t[end] == t_perturb && all(isfinite, pre_state)
        fig = Figure(; size=(1600,800))
        ax1 = Axis(fig[1,1];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "Before Press (Full Model)" :
                                           "Before Press (Simplified)")
        ax2 = Axis(fig[1,2];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "After Press (Full Model)" :
                                           "After Press (Simplified)")
        for i in 1:R
            lines!(ax1, sol1.t[1:end-1], sol1[i, 1:end-1], color=:blue)
            lines!(ax2, sol2.t, sol2[i, :], color=:blue)
        end
        for i in R+1:R+C
            lines!(ax1, sol1.t[1:end-1], sol1[i, 1:end-1], color=:red)
            lines!(ax2, sol2.t, sol2[i, :], color=:red)
        end
        display(fig)
    end

    return return_times, overshoot, ire, before_persistence, new_equil
end

# -------------------------------
# Example parameter definitions
if false
    # Define target equilibrium abundances
    R_eq = fill(5.0, R)  # Equilibrium biomass for each resource
    C_eq = fill(1.0, C)  # Equilibrium biomass for each consumer

    # Compute the new parameters
    new_xi_cons, new_r_res = calibrate_parametrization(R_eq, C_eq, p_tuple)

    println("New consumer thresholds (xi_cons): ", new_xi_cons)
    println("New resource growth rates (r_res): ", new_r_res)

    # Numbers of resources and consumers
    R = 2      # e.g., two resources
    C = 3      # e.g., three consumers

    # Consumer parameters (vectors of length C)
    m_cons  = ones(C)             # Mortality rates (example: all ones)
    xi_cons = ones(C)             # Overpopulation thresholds

    # Resource parameters (vectors of length R)
    r_res = ones(R)               # Intrinsic growth rates
    d_res = ones(R)               # Scaling factors for resources

    # Conversion efficiency
    epsilon = 0.5

    # Global interaction matrix A (size: (R+C) x (R+C))
    # Here we create a random matrix for demonstration.
    A = randn(R+C, R+C)

    # Bundle parameters into a tuple in the specified order
    p_tuple = (R, C, m_cons, new_xi_cons, new_r_res, d_res, epsilon, A)
    # -------------------------------
    # Initial conditions:
    # Resources start at 5.0, consumers at 1.0
    u0 = [ones(R) .* 5.0; ones(C) .* 1.0]
    tspan = (0.0, 50.0)

    # Set up and solve the ODE problem
    prob = ODEProblem(trophic_ode!, u0, tspan, p_tuple)
    sol = solve(prob, Tsit5())

    # Plot the solution
    # Labels: first R are resources, next C are consumers.
    begin
        fig = Figure(; size=(800, 600), title="Trophic Dynamics")
        ax = Axis(fig[1, 1], xlabel="Time", ylabel="Population Size", title="Trophic Dynamics")
        for i in 1:(R)
            MK.lines!(ax, sol.t, sol[i, :], label="Resource $i", color=:blue, linewidth=1)
        end
        for i in 1:(C)
            MK.lines!(ax, sol.t, sol[R+i, :], label="Consumer $i", color=:red, linewidth=1)
        end
        display(fig)
    end
end
#############################################################################
#############################################################################
################## Build Jacobian & Analytical V ############################
#############################################################################
#############################################################################
function build_jacobian(B_eq, p)
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    total = R + C
    J = zeros(total, total)

    # Resource block
    for i in 1:R
        # diff(du_i)/diffB_i = - d_res[i] * B_eq[i]
        J[i,i] = -d_res[i]*B_eq[i]
        # diff(du_i)/diffB_j (j a consumer)
        for j in (R+1):total
            J[i,j] =  d_res[i]*B_eq[i]*A[i,j]
        end
    end

    # Consumer block
    for k in 1:C
        i = R + k
        psi = B_eq[i] / xi_cons[k]
        prefactor = m_cons[k] * psi
        for j in 1:total
            # -delta_{ij} + epsilon[i,j]*A[i,j]  - A[j,i]
            J[i,j] = prefactor * ( (i==j ? -1 : 0) + epsilon[i,j]*A[i,j] - A[j,i] )
        end
    end

    return J
end

"""
    compute_analytical_V(J, R, C, m_cons, xi_cons)

Return V = J^{-1} * D, with D_{ii}=m_i/xi_i for consumers (i>R), D_{ii}=0 for resources.
"""
function compute_analytical_V(J, R, C, m_cons, xi_cons)
    total = size(J,1)
    D = zeros(total, total)
    for k in 1:C
        D[R+k, R+k] = m_cons[k] / xi_cons[k]
    end
    # Solve J * V = D  ?  V = J \ D
    return J \ D
end

##########################################################################
##########################################################################
################### Ladder-of-simplification transforms ##################
##########################################################################
##########################################################################
function transform_for_ladder_step(step, A_adj, epsilon_full)
    total = size(A_adj,1)

    if step == 1 # Full
        return A_adj, epsilon_full
    elseif step == 2 # S2 epsilon mean per row
        epsilon2 = similar(epsilon_full)
        for i in 1:total
            epsilon2[i,:] .= mean(epsilon_full[i,:])
        end
        return A_adj, epsilon2
    elseif step == 3 # S3 epsilon global mean
        return A_adj, fill(mean(epsilon_full), total, total)
    elseif step == 4 # S4 epsilon re-randomised
        epsilon4 = similar(epsilon_full)
        for i in 1:total, j in 1:total
            epsilon4[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
        end
        return A_adj, epsilon4
    elseif 5 <= step <= 8 
        
        pos, neg = mean(A_adj[A_adj.>0]), mean(A_adj[A_adj.<0])
        
        A5 = map(x-> x>0 ? pos : x<0 ? neg : 0, A_adj)
        if step==5 # epsilon is full
            return A5, epsilon_full
        elseif step==6 # epsilon is mean per row
            epsilon6 = similar(epsilon_full)
            for i in 1:total 
                epsilon6[i,:] .= mean(epsilon_full[i,:])
            end
            return A5, epsilon6
        elseif step==7 # epsilon is global mean
            epsilon7 = fill(mean(epsilon_full), total, total)
            return A5, epsilon7
        elseif step==8 # epsilon is re-randomised
            epsilon8 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon8[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A5, epsilon8
        end
    elseif 9 <= step <= 12 # A is averaged across all non-zero values
        
        m = mean(abs.(A_adj[A_adj .!= 0]))
        A6 = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)
        
        if step==9 # epsilon is full
            return A6, epsilon_full
        elseif step==10 # epsilon is mean per row
            epsilon10 = similar(epsilon_full)
            for i in 1:total
                epsilon10[i,:] .= mean(epsilon_full[i,:])
            end
            return A6, epsilon10
        elseif step==11 # epsilon is global mean
            return A6, fill(mean(epsilon_full), total, total)
        elseif step==12 # epsilon is re-randomised
            epsilon12 = similar(epsilon_full)
            for i in 1:total*total
                epsilon12[i] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A6, epsilon12
        end
    elseif 13 <= step <= 16 # A is fully randomised with same sparsity pattern as A_adj
        # Fully randomized A matrix using Normal(mean, std) of abs non-zero A entries
        A_vals = abs.(A_adj[A_adj .!= 0])
        A_mean, A_std = mean(A_vals), std(A_vals)
    
        A_rand = similar(A_adj)
        for i in 1:total, j in 1:total
            if A_adj[i, j] != 0
                A_rand[i, j] = rand(Normal(A_mean, A_std)) * sign(A_adj[i, j])
            else
                A_rand[i, j] = 0.0
            end
        end
    
        if step == 13 # epsilon is full
            return A_rand, epsilon_full
        elseif step == 14 # epsilon is mean per row
            epsilon14 = similar(epsilon_full)
            for i in 1:total
                epsilon14[i, :] .= mean(epsilon_full[i, :])
            end
            return A_rand, epsilon14
        elseif step == 15 # epsilon is global mean
            return A_rand, fill(mean(epsilon_full), total, total)
        elseif step == 16 # epsilon is re-randomised
            epsilon16 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon16[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A_rand, epsilon16
        end    
    else
        error("Unknown ladder step $step")
    end
end

#########################################################################
#########################################################################
######################### BUILD CALLBACKS ###############################
#########################################################################
#########################################################################
function build_callbacks(S::Int, EXTINCTION_THRESHOLD::Float64)
    # A) Continuous threshold-based extinctions
    callbacks = []

    # 1) Always ensure positivity
    push!(callbacks, PositiveDomain())

    # 2) Herbivores: set to zero if below EXTINCTION_THRESHOLD
    for x in 1:S
        function threshold_condition(u, t, integrator)
            # event if u[x] < EXTINCTION_THRESHOLD
            return u[x] - EXTINCTION_THRESHOLD
        end
        function threshold_affect!(integrator)
            integrator.u[x] = 0.0
        end
        push!(callbacks, ContinuousCallback(threshold_condition, threshold_affect!))
    end

    # Build callback set => no forced extinction
    cb_no_trigger = CallbackSet(callbacks...)

    return cb_no_trigger
end

# -----------------------------------------------------------
# Type II Trophic ODE
# -----------------------------------------------------------
function trophic_two!(du, u, p, t)
    R, C, m_cons, xi_cons, r_res, d_res, ε, A, h = p
    total = R + C

    # resources (same as before)
    for i in 1:R
      pred_sum = sum(A[i,j]*u[j] for j in R+1:total)
      du[i] = u[i]*d_res[i]*((r_res[i]/d_res[i]) - u[i] + pred_sum)
    end

    # consumers with Type II
    for k in 1:C
      i = R + k

      # handling‐denominator D_i = 1 + ∑ₖ h_{i,k}·A_{i,k}·u[k]  (only positive A links)
      D = 1 + sum(h[i,j]*A[i,j]*u[j] for j in 1:total if A[i,j] > 0; init = 0.0)

      # saturated prey‐gain
      gain = sum(ε[i,j]*A[i,j]*u[j] for j in 1:total if A[i,j] > 0; init = 0.0) / D

      # linear predation‐loss
      loss = sum(A[i,j]*u[j] for j in 1:total if A[i,j] < 0; init = 0.0)

      du[i] = u[i] * (m_cons[k]/xi_cons[k]) * (-xi_cons[k] - u[i] + gain + loss)
    end
end

# -----------------------------------------------------------
# Calibration for Type II
# -----------------------------------------------------------
function calibrate_two(
    R_eq::Vector{Float64},
    C_eq::Vector{Float64},
    p_calib;
    xi_threshold = 0.3
)
    R, C, m_cons, d_res, ε, A, h = p_calib
    total = R + C
    B_eq = vcat(R_eq, C_eq)

    new_xi = similar(C_eq)
    new_r  = similar(R_eq)

    # 1) consumers: xi solves
    #    0 = -xi - B_eq[i] + ∑ ε A B / (1+∑ h A B)  - ∑ (loss)   ⇒ xi = -B_eq + gain - loss
    for k in 1:C
      i = R + k
      D = 1 + sum(h[i,j]*A[i,j]*B_eq[j] for j in 1:total if A[i,j]>0; init=0.0)
      gain = sum(ε[i,j]*A[i,j]*B_eq[j] for j in 1:total if A[i,j]>0; init=0.0) / D
      loss = sum(A[i,j]*B_eq[j] for j in 1:total if A[i,j]<0; init=0.0)
      xi  = -C_eq[k] + gain - loss
      new_xi[k] = (xi>0 && C_eq[k]/xi < xi_threshold) ? xi : NaN
    end

    # 2) resources: same as before
    for i in 1:R
      loss = sum(A[i,j]*B_eq[j] for j in R+1:R+C if A[i,j]<0; init=0.0)
      new_r[i] = max( d_res[i] * (R_eq[i] + loss), 0 )
    end

    return new_xi, new_r
end