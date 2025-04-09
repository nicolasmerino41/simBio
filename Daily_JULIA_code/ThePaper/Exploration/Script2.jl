#############################
# Julia Pipeline for Incorporating Abundance Distributions into Stability Predictions
# Inspired by Barbier et al. (2018)
#############################
# ========= Step 1: Simulate the Full Complex Model ===========
# We use a generalized Lotka–Volterra model:
#   d𝑥[i]/dt = 𝑥[i] * ( r[i] + sum_{j=1}^n A[i,j] * 𝑥[j] )
# where r[i] can be interpreted as intrinsic growth rate (or carrying capacity when D = 1),
# and A[i,j] are the interaction coefficients.
#
# For our study, we choose:
# - r[i] drawn from a LogNormal distribution (to yield fat-tailed abundances)
# - Self-regulation is imposed by setting diagonal of A to -1
# - Off-diagonals of A are drawn from a Normal distribution with a small variance.

Random.seed!(1234)  # For reproducibility
n = 50  # Number of species

# Generate intrinsic growth rates r[i] from a lognormal distribution.
μ_r = 0.0           # mean parameter for log-normal (in log space)
σ_r = 0.5           # standard deviation in log space
r_full = rand(LogNormal(μ_r, σ_r), n)

# Self-regulation parameters D[i] (set to 1 for simplicity)
D_full = ones(n)

# Construct the interaction matrix A_full:
# - Diagonals: set to -1 (self-limitation)
# - Off-diagonals: drawn from Normal(0, 0.1)
A_full = zeros(n, n)
for i in 1:n
    for j in 1:n
        if i == j
            A_full[i, j] = -1.0
        else
            A_full[i, j] = rand(Normal(0, 0.1))
        end
    end
end

# Define the GLV ODE for the full model.
function full_model!(du, u, p, t)
    r, A = p
    n = length(u)
    for i in 1:n
        # Sum interaction from all species j
        interaction = sum(A[i, :] .* u)
        du[i] = u[i] * (r[i] + interaction)
    end
end

# Set initial conditions; here we start at half the carrying capacity r_full.
u0_full = 0.5 .* r_full
tspan = (0.0, 100.0)
p_full = (r_full, A_full)

# Setup and solve the ODE problem.
prob_full = ODEProblem(full_model!, u0_full, tspan, p_full)
sol_full = solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)
equilibrium_full = sol_full.u[end]  # Equilibrium abundances at final time

println("=== Full Model Equilibrium Abundances ===")
println(equilibrium_full)

# ========= Step 2: Compute Community Metrics from the Full Model ===========
threshold = 1e-3  # Abundance threshold to consider a species as surviving.
survivors_idx = findall(x -> x > threshold, equilibrium_full)
φ_full = length(survivors_idx) / n  # Persistence: fraction of surviving species
total_biomass_full = sum(equilibrium_full)
mean_abundance_full = mean(equilibrium_full[survivors_idx])
std_abundance_full = std(equilibrium_full[survivors_idx])

println("\n--- Full Model Community Metrics ---")
println("Persistence (φ): ", φ_full)
println("Total Biomass: ", total_biomass_full)
println("Mean Abundance (survivors): ", mean_abundance_full)
println("Std Abundance (survivors): ", std_abundance_full)

# ========= Step 3: Fit the Abundance Distribution ===========
# We now take the surviving species’ abundances and fit a lognormal distribution.
survivor_abundances = equilibrium_full[survivors_idx]

# Fit the lognormal to survivors using maximum likelihood estimation.
lognormal_fit = fit(LogNormal, survivor_abundances)
μ_fit = lognormal_fit.μ       # Mean of log(abundance)
σ_fit = lognormal_fit.σ       # Std dev of log(abundance)

# Compute the moments of the fitted lognormal distribution.
fitted_mean = exp(μ_fit + σ_fit^2 / 2)   # E[B] for abundance B
fitted_var = exp(2*μ_fit + σ_fit^2) * (exp(σ_fit^2) - 1)

println("\n--- Lognormal Fit of Survivors ---")
println("Fitted μ (log-space): ", μ_fit)
println("Fitted σ (log-space): ", σ_fit)
println("Fitted Mean Abundance: ", fitted_mean)
println("Fitted Variance: ", fitted_var)

# ========= Step 4: Compute Emergent Parameters ===========
# In the reduced (disordered) model, the key emergent parameters include:
# ζ = std(K) where K = r (since D=1 in our simulation)
# μ_eff = n * mean(off-diagonal elements of A_full)
# σ_eff = sqrt(n) * std(off-diagonals of A_full)
# γ_eff = correlation between A[i,j] and A[j,i] (for i ≠ j)

K_full = r_full  # Since D_full=1, effective carrying capacities K = r.
ζ = std(K_full)

# Extract off-diagonal values.
offdiag_vals = [A_full[i,j] for i in 1:n for j in 1:n if i != j]
μ_eff = n * mean(offdiag_vals)
σ_eff = sqrt(n) * std(offdiag_vals)

# Compute reciprocity: correlation between A[i,j] and A[j,i] for i < j.
pairs = [(A_full[i,j], A_full[j,i]) for i in 1:n for j in 1:n if i < j]
vals1 = [p[1] for p in pairs]
vals2 = [p[2] for p in pairs]
γ_eff = cor(vals1, vals2)

println("\n--- Emergent Parameters (Reference Model) ---")
println("ζ (std of K): ", ζ)
println("μ_eff: ", μ_eff)
println("σ_eff: ", σ_eff)
println("γ_eff: ", γ_eff)
println("Persistence (φ): ", φ_full)

# ========= Step 5: Conduct a Perturbation Experiment ===========
# We define a function to simulate a perturbation at a specified time, then measure the recovery time
# (return time) for each species – defined here as the time to recover to within 10% of its pre-perturbation value.

function simulate_perturbation(u0, p, tspan, perturb_time, perturb_factor)
    # Integrate to the perturbation time.
    tspan1 = (tspan[1], perturb_time)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8)
    state_pre = sol1.u[end]
    
    # Apply the perturbation (e.g. reduce abundances by perturb_factor).
    state_perturbed = state_pre .* perturb_factor
    
    # Integrate from the perturbation time to tspan end.
    tspan2 = (perturb_time, tspan[2])
    prob2 = ODEProblem(full_model!, state_perturbed, tspan2, p)
    sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    # For each species, record the time until its abundance recovers
    # to within 10% of its value before the perturbation.
    return_times = zeros(n)
    for i in 1:n
        target = state_pre[i]
        recovered = false
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target) / (abs(target) + 1e-8) < 0.1
                return_times[i] = t - perturb_time
                recovered = true
                break
            end
        end
        if !recovered
            return_times[i] = NaN
        end
    end
    return return_times
end

perturb_time = 50.0
perturb_factor = 0.5  # Reduce each abundance by 50% at perturb_time.
rt_full = simulate_perturbation(u0_full, p_full, tspan, perturb_time, perturb_factor)
mean_return_time_full = mean(filter(!isnan, rt_full))
println("\nMean Return Time after Perturbation (Full Model): ", mean_return_time_full)

# ========= Step 6: Construct the Simplified (Reference) Model ===========
# The idea behind the disordered (reference) model is to remove detailed structure.
# Here, we replace all off-diagonal interactions with their overall mean.
A_simpl = zeros(n, n)
for i in 1:n
    for j in 1:n
        if i == j
            A_simpl[i,j] = -1.0  # Preserve self-limitation.
        else
            A_simpl[i,j] = mean(offdiag_vals)
        end
    end
end

p_simpl = (r_full, A_simpl)
u0_simpl = u0_full  # Use the same initial condition.
prob_simpl = ODEProblem(full_model!, u0_simpl, tspan, p_simpl)
sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
equilibrium_simpl = sol_simpl.u[end]

println("\n=== Simplified Model Equilibrium Abundances ===")
println(equilibrium_simpl)

# Compute simplified model community metrics.
survivors_idx_simpl = findall(x -> x > threshold, equilibrium_simpl)
φ_simpl = length(survivors_idx_simpl) / n
total_biomass_simpl = sum(equilibrium_simpl)
mean_abundance_simpl = mean(equilibrium_simpl[survivors_idx_simpl])
println("\n--- Simplified Model Community Metrics ---")
println("Persistence (φ): ", φ_simpl)
println("Total Biomass: ", total_biomass_simpl)
println("Mean Abundance (survivors): ", mean_abundance_simpl)

# Perturbation for simplified model.
rt_simpl = simulate_perturbation(u0_simpl, p_simpl, tspan, perturb_time, perturb_factor)
mean_return_time_simpl = mean(filter(!isnan, rt_simpl))
println("\nMean Return Time after Perturbation (Simplified Model): ", mean_return_time_simpl)

# ========= Step 7: Incorporate the Abundance Distribution into Stability Predictions ===========
# Our goal is now to modify a stability (return time) prediction function using the
# fitted abundance distribution. The fitted lognormal provides effective moments that
# may differ from the Gaussian assumption.
#
# Here, we compute the relative variance of abundances (variance divided by mean²) from the lognormal fit.
rel_variance = fitted_var / fitted_mean^2
# We now define a hypothetical stability prediction that uses φ, γ_eff and the relative variance.
# (For example, the predicted return time might be inversely proportional to [1 - φ * γ_eff * rel_variance].)
# This formula is illustrative and can be refined as needed.
function predicted_return_time(φ, γ, rel_var)
    denom = 1 - φ * γ * rel_var
    return denom > 0 ? 1 / denom : Inf
end

T_pred = predicted_return_time(φ_full, γ_eff, rel_variance)
println("\nPredicted Return Time from Modified Reference Model (using abundance distribution): ", T_pred)

# ========= Step 8: Compare Stability Measurements ===========
println("\n=== Comparison of Return Times ===")
println("Full Model Mean Return Time: ", mean_return_time_full)
println("Simplified Model Mean Return Time: ", mean_return_time_simpl)
println("Predicted Return Time (Reference Model): ", T_pred)

# ========= Step 9: Visualization ===========
# Plot the histogram of survivor abundances and overlay the fitted lognormal PDF.
MK.hist(survivor_abundances, bins=20,
        # normalize=true,
        #   title="Abundance Distribution of Survivors",
        #   xlabel="Abundance", ylabel="Density", label="Data"
)
x_vals = range(minimum(survivor_abundances), stop=maximum(survivor_abundances), length=100)
pdf_vals = pdf.(lognormal_fit, x_vals)
begin
    fig = Figure(; size = (600, 400))
    ax = Axis(fig[1, 1], title="Abundance Distribution of Survivors",
              xlabel="Abundance", ylabel="Density")
    MK.lines!(ax, x_vals, pdf_vals)
    display(fig)
end

# Plot time series (optional) for full and simplified models.
begin
    fig = Figure(; size = (600, 400))
    ax1 = Axis(fig[1, 1], title="Time Series: Full Model", xlabel="Time", ylabel="Species Abundance")
    ax2 = Axis(fig[2, 1], title="Time Series: Simplified Model", xlabel="Time", ylabel="Species Abundance")
    for i in 1:n
        lines!(ax1, sol_full.t, sol_full[i, :], label="Species $i")
    end
    for i in 1:n
        lines!(ax2, sol_simpl.t, sol_simpl[i, :], label="Species $i")
    end
    display(fig)
end