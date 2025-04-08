using DifferentialEquations, Random, Statistics, Plots

# ----------------------------
# Step 1: Define the Full Model
# ----------------------------
# The full model is defined by:
#   dx_i/dt = x_i * ( r_i + Σ_j a_ij * x_j )
# where r is a vector of intrinsic growth rates and A is the interaction matrix.

# Set seed for reproducibility
Random.seed!(1234)

# Number of species
n = 20

# Define intrinsic growth rates: sample uniformly between 0.5 and 1.5
r_full = rand(n) .+ 0.5

# Construct full interaction matrix A:
#   - Set self-interaction (diagonal) to -1.0 (self-limitation)
#   - Off-diagonals are drawn from a normal distribution (mean 0, sd 0.1)
A_full = zeros(n, n)
for i in 1:n, j in 1:n
    if i == j
        A_full[i, j] = -1.0
    else
        A_full[i, j] = randn() * 0.1
    end
end

# Define the full model ODE function
function full_model!(du, u, p, t)
    r, A = p
    for i in 1:length(u)
        s = 0.0
        for j in 1:length(u)
            s += A[i, j] * u[j]
        end
        du[i] = u[i] * (r[i] + s)
    end
end

# Initial conditions: start with abundance 1 for each species
u0 = ones(n)
# Time span for simulation
tspan = (0.0, 50.0)

# Define parameters as a tuple and set up the ODE problem
p_full = (r_full, A_full)
prob_full = ODEProblem(full_model!, u0, tspan, p_full)

# Solve the full model
sol_full = solve(prob_full, Tsit5())

# Extract time-series data for later analysis
data_full = Array(sol_full)

# ----------------------------
# Step 2: Compute Metrics for the Full Model
# ----------------------------

# Define a threshold for persistence (e.g., abundance > 0.1)
threshold = 0.1
persistence_full = sum(data_full[end, :] .> threshold) / n

# Compute fluctuations: average standard deviation over the simulation for each species
fluctuations_full = mean([std(data_full[:, i]) for i in 1:n])

println("Full Model Metrics:")
println("  Persistence (fraction above threshold): ", persistence_full)
println("  Average fluctuation (std) of species abundances: ", fluctuations_full)

# Function to simulate a perturbation and estimate return times
function simulate_perturbation(u0, p, tspan, perturb_time, perturbation_factor)
    # Simulate until the perturbation time
    tspan1 = (tspan[1], perturb_time)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, Tsit5())
    
    # Apply a perturbation: multiply species abundances by perturbation_factor
    u_perturb = sol1.u[end] .* perturbation_factor

    # Continue simulation after perturbation
    tspan2 = (perturb_time, tspan[2])
    prob2 = ODEProblem(full_model!, u_perturb, tspan2, p)
    sol2 = solve(prob2, Tsit5())
    
    # Define recovery as reaching within 10% of pre-perturbation abundance
    target = sol1.u[end]
    recovery_times = zeros(n)
    for i in 1:n
        recovered = false
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target[i]) / target[i] < 0.1
                recovery_times[i] = t - perturb_time
                recovered = true
                break
            end
        end
        if !recovered
            recovery_times[i] = NaN
        end
    end
    return recovery_times
end

# Simulate a perturbation for the full model:
#   - Perturb at time 25.0 by reducing abundances by 50%
perturb_time = 25.0
perturbation_factor = 0.5
recovery_times_full = simulate_perturbation(u0, p_full, tspan, perturb_time, perturbation_factor)
mean_return_time_full = mean(skipmissing(recovery_times_full))
println("  Mean Return Time after perturbation: ", mean_return_time_full)

# ----------------------------
# Step 3: Simplification of the Model
# ----------------------------
# Construct a simplified interaction matrix that removes detailed pairwise information.
# Here we replace all off-diagonal elements with their overall mean.
offdiag_vals = [A_full[i, j] for i in 1:n, j in 1:n if i != j]
mean_off = mean(offdiag_vals)

A_simpl = fill(mean_off, n, n)
for i in 1:n
    A_simpl[i, i] = -1.0  # Preserve self-limitation in the diagonal
end

# Use the same growth rate vector
p_simpl = (r_full, A_simpl)

# ----------------------------
# Step 4: Simulate and Compute Metrics for the Simplified Model
# ----------------------------
prob_simpl = ODEProblem(full_model!, u0, tspan, p_simpl)
sol_simpl = solve(prob_simpl, Tsit5())
data_simpl = Array(sol_simpl)

# Compute persistence for the simplified model
persistence_simpl = sum(data_simpl[end, :] .> threshold) / n
# Compute average fluctuations
fluctuations_simpl = mean([std(data_simpl[:, i]) for i in 1:n])

println("\nSimplified Model Metrics:")
println("  Persistence (fraction above threshold): ", persistence_simpl)
println("  Average fluctuation (std) of species abundances: ", fluctuations_simpl)

# Simulate the perturbation for the simplified model
recovery_times_simpl = simulate_perturbation(u0, p_simpl, tspan, perturb_time, perturbation_factor)
mean_return_time_simpl = mean(skipmissing(recovery_times_simpl))
println("  Mean Return Time after perturbation: ", mean_return_time_simpl)

# ----------------------------
# Step 5: Compare Full vs. Simplified Models
# ----------------------------
println("\nComparison Summary:")
println("  Persistence (Full vs. Simplified): ", persistence_full, " vs ", persistence_simpl)
println("  Average Fluctuation (Full vs. Simplified): ", fluctuations_full, " vs ", fluctuations_simpl)
println("  Mean Return Time (Full vs. Simplified): ", mean_return_time_full, " vs ", mean_return_time_simpl)

# ----------------------------
# Optional: Plotting the Time Series
# ----------------------------
# Plot the time series of species abundances for the full model
plt1 = plot(sol_full, title="Full Model Time Series", xlabel="Time", ylabel="Species Abundance")
# Plot the time series for the simplified model
plt2 = plot(sol_simpl, title="Simplified Model Time Series", xlabel="Time", ylabel="Species Abundance")
plot(plt1, plt2, layout=(2,1))
