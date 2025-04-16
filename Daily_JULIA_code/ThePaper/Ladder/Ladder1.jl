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
        for j in R+1 : R+C
            pred_sum += A[j, i] * u[j]
        end
        # Resource i dynamics
        du[i] = u[i] * d_res[i] * ( (r_res[i] / d_res[i]) - u[i] - pred_sum )
    end

    # Consumers dynamics (indices R+1 to R+C)
    for k in 1:C
        i = R + k  # global index for consumer k
        # Prey term: sum over resources (indices 1 to R)
        sum_prey = 0.0
        for j in 1:R
            sum_prey += epsilon * A[i, j] * u[j]
        end
        # Predator term: sum over consumers (indices R+1 to R+C)
        sum_pred = 0.0
        for j in R+1 : R+C
            sum_pred += A[j, i] * u[j]
        end
        du[i] = u[i] * (m_cons[k] / xi_cons[k]) * ( -xi_cons[k] - u[i] + sum_prey - sum_pred )
    end
end

#############################################################################
#############################################################################
################## CALIBRATE PARAMETRISATION ################################
#############################################################################
#############################################################################
# --- Revised Calibration Function ---
# p_calib = (R, C, m_cons, d_res, ε, A)
function calibrate_params(R_eq::Vector{Float64}, C_eq::Vector{Float64}, p_calib)
    R, C, m_cons, d_res, epsilon, A_adjacency = p_calib

    viable = false
    count = 0
    while !viable
        # Create a local copy of the adjacency matrix.
        A_local = copy(A_adjacency)
        new_xi_cons = zeros(C)
        new_r_res = zeros(R)
        for i in 1:(R+C), j in 1:(R+C)
            A_local[i,j] = A_adjacency[i,j] * rand()  # Randomize the strength
        end

        # Consumer calibration:
        for k in 1:C
            i = R + k
            prey_sum = sum(epsilon * A_local[i, j] * R_eq[j] for j in 1:R)
            pred_sum = sum(A_local[j, i] * C_eq[j - R] for j in (R+1):(R+C))
            ξ_val = - C_eq[k] + prey_sum - pred_sum
            new_xi_cons[k] = (ξ_val > 0) ? ξ_val : NaN
        end

        # Resource calibration:
        for i in 1:R
            pred_sum = sum(A_local[j, i] * C_eq[j - R] for j in (R+1):(R+C))
            r_val = d_res[i] * (R_eq[i] + pred_sum)
            new_r_res[i] = (r_val > 0) ? r_val : NaN
        end

        count += 1
        if (any(isnan, new_xi_cons) || any(isnan, new_r_res)) && count < 10
            viable = false
        elseif !any(isnan, new_xi_cons) && !any(isnan, new_r_res)
            viable = true
        elseif (any(isnan, new_xi_cons) || any(isnan, new_r_res)) && count >= 10
            # println("Tried 10 times and still not viable; breaking loop")
            viable = true
        end
        if viable
            return new_xi_cons, new_r_res, A_local
        end
    end
end


# Define target equilibrium abundances
R_eq = fill(5.0, R)  # Equilibrium biomass for each resource
C_eq = fill(1.0, C)  # Equilibrium biomass for each consumer

# Compute the new parameters
new_xi_cons, new_r_res = calibrate_parametrization(R_eq, C_eq, p_tuple)

println("New consumer thresholds (xi_cons): ", new_xi_cons)
println("New resource growth rates (r_res): ", new_r_res)

# -------------------------------
# Example parameter definitions

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
p_tuple = (R, C, m_cons, xi_cons, r_res, d_res, epsilon, A)
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