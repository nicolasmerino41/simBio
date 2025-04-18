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
    # u[1:R]          → resources
    # u[R+1 : R+C]  → consumers
    
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
# --- Revised Calibration Function ---
# p_calib = (R, C, m_cons, d_res, ε, A)
function calibrate_params(
    R_eq::Vector{Float64},
    C_eq::Vector{Float64},
    p_calib
)
    R, C, m_cons, d_res, epsilon, A_adj = p_calib
    total = R + C
    B_eq = vcat(R_eq, C_eq)

    viable = false
    count  = 0

    while !viable
        # 1) Randomize strengths on the +1/–1 adjacency
        A_used = zeros(total, total)
        for i in 1:total, j in 1:total
            if A_adj[i,j] == +1
                A_used[i,j] =  rand()     # positive link strength
            elseif A_adj[i,j] == -1
                A_used[i,j] = -rand()     # negative link strength
            end
        end

        # 2) Prepare output vectors
        total = R + C
        B_eq = vcat(R_eq, C_eq)

        local new_xi_cons = zeros(C)
        local new_r_res   = zeros(R)

        # 1) Consumer thresholds
        for k in 1:C
            i = R + k

            # Gains from prey (only A_used[i,j]>0)
            prey_sum = sum(
                epsilon[i,j] * A_used[i,j] * B_eq[j]
                for j in 1:total if A_used[i,j] > 0.0;
                init = 0.0
            )

            # Losses to predators (only A_used[i,j]<0)
            pred_sum = sum(
                A_used[i,j] * B_eq[j]
                for j in 1:total if A_used[i,j] < 0.0;
                init = 0.0
            )

            ξ = -C_eq[k] + prey_sum + pred_sum
            new_xi_cons[k] = ξ > 0 ? ξ : NaN
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

        # 5) Viability check
        count += 1
        if (!any(isnan, new_xi_cons) && !any(isnan, new_r_res)) || count ≥ 10
            viable = true
            if count ≥ 10 && (any(isnan, new_xi_cons) || any(isnan, new_r_res))
                @warn "calibrate_params: gave up after 10 tries"
            end
        end
        if viable
            return new_xi_cons, new_r_res, A_used
        end
    end
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
        # Diagonal: ∂R_i/∂R_i
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
        α = m_cons[k] / xi_cons[k]
        # Diagonal (∂C_k/∂C_k)
        J[i, i] = - α * B[i]
        # Derivatives with respect to resources (indices 1:R)
        for j in 1:R
            J[i, j] = B[i] * α * epsilon[i, j] * A[i, j]
        end
        # Derivatives with respect to other consumers (indices R+1 to total)
        for j in (R+1):total
            if j != i
                J[i, j] = - B[i] * α * A[j, i]
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
function simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=false, show_warnings=true, full_or_simple = true)
    # Unpack p: (R, C, m_cons, xi_cons, r_res, d_res, epsilon, A)
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p

    # Phase 1: simulate from tspan[1] to t_perturb.
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(trophic_ode!, u0, tspan1, p)
    if show_warnings
        sol1 = solve(prob1, solver; reltol=1e-8, abstol=1e-8)
    else
        sol1 = with_logger(logger) do
            solve(prob1, solver; reltol=1e-8, abstol=1e-8)
        end
    end
    pre_state = sol1.u[end]
    before_persistence = sum(x -> x > EXTINCTION_THRESHOLD, pre_state) / length(pre_state)
    # Phase 2: apply the press perturbation (increase thresholds by delta)
    xi_press = copy(xi_cons)
    xi_press .= xi_press .* (1 - delta)
    # Build updated parameter tuple.
    p_press = (R, C, m_cons, xi_press, r_res, d_res, epsilon, A)

    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(trophic_ode!, pre_state, tspan2, p_press)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)

    # Compute return times: time for each species to come within 10% of new equilibrium.
    return_times = zeros(n)
    for i in 1:n
        target = new_equil[i]
        recovered = false
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target) / (abs(target) + 1e-8) < 0.1
                return_times[i] = t - t_perturb
                recovered = true
                break
            end
        end
        if !recovered
            return_times[i] = NaN
        end
    end
    
    if plot && sol1.t[end] == t_perturb && all(!isnan, sol1.u[end]) && all(!isinf, sol1.u[end])
        fig = Figure(; size =(1600,800))
        ax1 = Axis(
            fig[1,1], xlabel="Time", ylabel="Abundance",
            title= full_or_simple ? "Community Response Before Press Perturbation (Full Model)" : "Community Response Before Press Perturbation (Simplified Model)"
        )
        ax2 = Axis(
            fig[1,2], xlabel="Time", ylabel="Abundance",
            title= full_or_simple ? "Community Response After Press Perturbation (Full Model)" : "Community Response After Press Perturbation (Simplified Model)"
        )
        for i in 1:R
            MK.lines!(ax1, sol1.t, sol1[i, :], label="Resource $i", color=:blue)
            MK.lines!(ax2, sol2.t, sol2[i, :], label="Resource $i", color=:blue)
        end
        for i in (R+1):(R+C)
            MK.lines!(ax1, sol1.t, sol1[i, :], label="Consumer $(i-R)", color=:red)
            MK.lines!(ax2, sol2.t, sol2[i, :], label="Consumer $(i-R)", color=:red)
        end
        # axislegend(ax, position=:rb)
        display(fig)
    end

    return return_times, before_persistence, new_equil
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