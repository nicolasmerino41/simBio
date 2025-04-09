#############################
# Julia Pipeline for Model Reduction and Stability Prediction
# Inspired by Barbier et al. (2018)
#############################

# --- Step 0. Set random seed for reproducibility ---
Random.seed!(1234)

# --- Step 1. Simulate the Full (Complex) Model ---
# We use a generalized Lotka–Volterra model:
#   dx[i]/dt = x[i] * ( r[i] - x[i] + sum_{j≠i} A[i,j] * x[j] )
#
# Here we generate r[i] from a lognormal distribution to help produce fat-tailed abundance distributions,
# and we assume a constant self-limitation (D = 1). Off-diagonal interactions A[i,j] are drawn from a Normal.
# (See Barbier et al. for the idea of using a “disordered limit” where only a few aggregate parameters matter.)

n = 50  # number of species

# Generate intrinsic growth rates (also viewed as carrying capacities since D=1)
mu_r = 0.0
sigma_r = 0.5
r_full = rand(LogNormal(mu_r, sigma_r), n)
D_full = ones(n)  # self-regulation parameters (set to 1)

# Construct interaction matrix A:
# Diagonal elements correspond to self-limitation, here set to -1.
# Off-diagonals are drawn from Normal(0, 0.1)
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

# Define the ODE for the full model
function full_model!(du, u, p, t)
    r, A = p
    n = length(u)
    for i in 1:n
        interaction = 0.0
        for j in 1:n
            interaction += A[i, j] * u[j]
        end
        du[i] = u[i] * (r[i] + interaction)
    end
end

# Set initial conditions as a fraction (e.g. 0.5) of the carrying capacities
u0_full = 0.5 .* r_full
tspan = (0.0, 100.0)
p_full = (r_full, A_full)
prob_full = ODEProblem(full_model!, u0_full, tspan, p_full)
sol_full = solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)

# Extract the equilibrium state (assumed to be reached at t = 100)
equilibrium_full = sol_full.u[end]
println("=== Full Model Equilibrium Abundances ===")
println(equilibrium_full)

# --- Step 2. Compute Community Metrics from the Full Model ---
threshold = 1e-3  # species are considered surviving if abundance > threshold
survivors_idx = findall(x -> x > threshold, equilibrium_full)
φ_full = length(survivors_idx) / n  # persistence = fraction of survivors
total_biomass_full = sum(equilibrium_full)
mean_abundance_full = mean(equilibrium_full[survivors_idx])
std_abundance_full = std(equilibrium_full[survivors_idx])

println("\n--- Full Model Community Metrics ---")
println("Persistence (φ): ", φ_full)
println("Total Biomass: ", total_biomass_full)
println("Mean Abundance (survivors): ", mean_abundance_full)
println("Std Abundance (survivors): ", std_abundance_full)

# --- Step 3. Compute Emergent Parameters ---
# Following the framework in Barbier et al., we derive:
#   - Effective carrying capacities: K[i] = r[i]/D[i] (with D = 1, so K = r)
#   - Interaction coefficients: α[i,j] = A[i,j] (for i ≠ j)
#
# The key emergent parameters are:
#   ζ = std(K)       (spread of carrying capacities)
#   μ_eff = n * mean(α_offdiag)
#   σ_eff = sqrt(n)*std(α_offdiag)
#   γ_eff = correlation between A[i,j] and A[j,i] (for i ≠ j)

K_full = r_full
ζ = std(K_full)

# Off-diagonal elements:
offdiag_vals = [A_full[i, j] for i in 1:n, j in 1:n if i != j]
μ_eff = n * mean(offdiag_vals)
σ_eff = sqrt(n) * std(offdiag_vals)

# Compute reciprocity (γ_eff) using pairs (i,j) with i < j
pairs = [(A_full[i,j], A_full[j,i]) for i in 1:n for j in 1:n if i < j]
vals1 = [p[1] for p in pairs]
vals2 = [p[2] for p in pairs]
γ_eff = cor(vals1, vals2)

println("\n--- Emergent Parameters (Reference Model) ---")
println("ζ (std of carrying capacities): ", ζ)
println("μ_eff: ", μ_eff)
println("σ_eff: ", σ_eff)
println("γ_eff: ", γ_eff)
println("Persistence (φ): ", φ_full)

# --- Step 4. Analyze the Abundance Distribution ---
# We now consider the distribution of equilibrium abundances for survivors.
survivor_abundances = equilibrium_full[survivors_idx]
log_abundances = log.(survivor_abundances)
mu_log = mean(log_abundances)
sigma_log = std(log_abundances)

println("\n--- Abundance Distribution Fit (Lognormal) ---")
println("mu_log: ", mu_log, ", sigma_log: ", sigma_log)
# One might also extract the tail exponent (if the distribution is fat-tailed) from a log-log plot of the CCDF.

# --- Step 5. Measure Stability via a Perturbation Experiment ---
# We define a function that perturbs the system (e.g. halves species abundances) at a given time,
# and then computes the recovery (return) time for each species (i.e. time to come within 10% of pre-perturbation value).

function simulate_perturbation(u0, p, tspan, perturb_time, perturb_factor)
    # Integrate up to the perturbation time.
    tspan1 = (tspan[1], perturb_time)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8)
    state_pre = sol1.u[end]
    
    # Apply perturbation (e.g., multiply abundances by perturb_factor)
    state_perturbed = state_pre .* perturb_factor
    
    # Integrate from perturbation time to end.
    tspan2 = (perturb_time, tspan[2])
    prob2 = ODEProblem(full_model!, state_perturbed, tspan2, p)
    sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    # For each species, record the time until its abundance recovers within 10% of the pre-perturbation value.
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
perturb_factor = 0.5  # for instance, reduce abundances by 50%
rt_full = simulate_perturbation(u0_full, p_full, tspan, perturb_time, perturb_factor)
mean_return_time_full = mean(filter(!isnan, rt_full))
println("\nMean Return Time after perturbation (Full Model): ", mean_return_time_full)

# --- Step 6. Construct the Simplified (Reference) Model ---
# To obtain the disordered limit, we “remove” detailed information by replacing each off-diagonal A[i,j]
# with its overall mean value. This yields a model parameterized solely by aggregate statistics.
A_simpl = zeros(n, n)
for i in 1:n
    for j in 1:n
        if i == j
            A_simpl[i, j] = -1.0
        else
            A_simpl[i, j] = mean(offdiag_vals)
        end
    end
end

p_simpl = (r_full, A_simpl)
u0_simpl = u0_full  # same initial conditions as before
prob_simpl = ODEProblem(full_model!, u0_simpl, tspan, p_simpl)
sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
equilibrium_simpl = sol_simpl.u[end]
println("\n=== Simplified Model Equilibrium Abundances ===")
println(equilibrium_simpl)

# Compute metrics for the simplified model
survivors_idx_simpl = findall(x -> x > threshold, equilibrium_simpl)
φ_simpl = length(survivors_idx_simpl) / n
total_biomass_simpl = sum(equilibrium_simpl)
mean_abundance_simpl = mean(equilibrium_simpl[survivors_idx_simpl])
std_abundance_simpl = std(equilibrium_simpl[survivors_idx_simpl])
println("\n--- Simplified Model Community Metrics ---")
println("Persistence (φ): ", φ_simpl)
println("Total Biomass: ", total_biomass_simpl)
println("Mean Abundance (survivors): ", mean_abundance_simpl)
println("Std Abundance (survivors): ", std_abundance_simpl)

# Perturbation experiment for simplified model
rt_simpl = simulate_perturbation(u0_simpl, p_simpl, tspan, perturb_time, perturb_factor)
mean_return_time_simpl = mean(skipmissing(rt_simpl))
println("\nMean Return Time after perturbation (Simplified Model): ", mean_return_time_simpl)

# --- Step 7. Predict Stability Attributes from Emergent Parameters & Abundance Distribution ---
# In the disordered (reference) model framework, many global properties (diversity, biomass, stability)
# are functions of a small set of aggregate parameters. Here we define a hypothetical function to predict
# the (return) time as a proxy for stability based on φ, μ_eff, σ_eff and γ_eff.
#
# (Note: This prediction function is illustrative. In Barbier et al., analytical formulas are derived via the cavity method.)
function predicted_return_time(φ, μ, σ, γ)
    denom = 1 - φ * γ * σ^2
    return denom > 0 ? 1 / denom : Inf
end

T_pred = predicted_return_time(φ_full, μ_eff, σ_eff, γ_eff)
println("\nPredicted Return Time from Reference Model (Hypothetical): ", T_pred)

# --- Step 8. Compare Stability Measurements ---
println("\n=== Comparison of Return Times ===")
println("Full Model Mean Return Time: ", mean_return_time_full)
println("Simplified Model Mean Return Time: ", mean_return_time_simpl)
println("Predicted Return Time (Reference Model): ", T_pred)

# --- Step 9. Visualize Abundance Distributions ---
MK.hist(
    survivor_abundances, bins=20,
    # title="Abundance Distribution (Full Model Survivors)",
    # xlabel="Abundance", ylabel="Frequency", label="Simulated",
    # legend=:topright
)

