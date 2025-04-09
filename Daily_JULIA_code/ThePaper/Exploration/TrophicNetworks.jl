#############################
# Julia Pipeline: Fixed-Abundance Trophic Community Study
#############################
# ---------------------------
# Function Definitions
# ---------------------------

# Full GLV predator-prey model.
# Here, we assume D[i] = 1 so that r[i] = k[i].
# The equation is:
#   dB[i]/dt = B[i] * ( r[i] - B[i] + sum_{j ≠ i} A[i,j] * B[j] )
function full_model!(du, B, p, t)
    r, A = p
    n = length(B)
    for i in 1:n
        du[i] = B[i] * ( r[i] - B[i] + sum(A[i, :] .* B) )
    end
end

# Function to simulate a perturbation:
# At time t_perturb, each species abundance is multiplied by perturb_factor.
# Return times are measured as the time required for each species to recover within 10% of pre-perturbation level.
function simulate_perturbation(u0, p, tspan, t_perturb, perturb_factor)
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]
    
    perturbed_state = pre_state .* perturb_factor
    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(full_model!, perturbed_state, tspan2, p)
    sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    n = length(pre_state)
    return_times = zeros(n)
    for i in 1:n
        target = pre_state[i]
        recovered = false
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i]-target) / (abs(target)+1e-8) < 0.1
                return_times[i] = t - t_perturb
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

# A hypothetical predicted return time function.
# Combines persistence, network reciprocity, and the relative variance of abundances.
function predicted_return_time(φ, γ_eff, relVar)
    denom = 1 - φ * γ_eff * relVar
    return denom > 0 ? 1 / denom : Inf
end

# ---------------------------
# Pipeline Parameters
# ---------------------------
species_scenarios = [20, 40, 60]   # Different numbers of species
Niter = 10                         # Iterations per scenario
tspan = (0.0, 100.0)               # Integration time span
t_perturb = 50.0                   # Time of perturbation
perturb_factor = 0.5               # Factor to multiply abundances at perturbation

# Containers for results
results = Dict{Int, Dict{Symbol, Vector{Float64}}}()

# ---------------------------
# Main Loop: Over Species Scenarios and Iterations
# ---------------------------
for n in species_scenarios
    persistence_arr = Float64[]
    rt_full_arr = Float64[]
    rt_simpl_arr = Float64[]
    pred_rt_arr = Float64[]
    relvar_arr = Float64[]
    
    for iter in 1:Niter
        #############################
        # Step 1: Generate Fixed Equilibrium Abundances
        #############################
        # Assume the observed equilibrium abundances B* come from a lognormal distribution.
        # These are fixed, empirical values.
        μ_log_obs = 0.5
        σ_log_obs = 1.0
        fixed_B = rand(LogNormal(μ_log_obs, σ_log_obs), n)
        
        #############################
        # Step 2: Construct a Trophic (Predator–Prey) Network
        #############################
        # Generate a directed interaction matrix A based on a connectance parameter.
        connectance = 0.2
        A = zeros(n, n)
        for i in 1:n
            for j in 1:n
                if i != j && rand() < connectance
                    # Randomly choose the direction: 50% chance species i preys on j.
                    if rand() < 0.5
                        β = rand(Exponential(1.0))       # interaction magnitude
                        A[i, j] = β                      # prey is harmed, predator gains
                        A[j, i] = -0.3 * β               # conversion efficiency (epsilon=0.3)
                    else
                        β = rand(Exponential(1.0))
                        A[j, i] = β
                        A[i, j] = -0.3 * β
                    end
                end
            end
        end
        
        #############################
        # Step 3: Calibrate Model with Fixed Abundances
        #############################
        # For each species i, adjust the carrying capacity to force the equilibrium:
        #   B_i* = k_i - sum_{j ≠ i} A[i,j] * B_j*
        # Thus set:
        #   k_i = B_i* + sum_{j ≠ i} A[i,j] * B_j*
        k = zeros(n)
        for i in 1:n
            k[i] = fixed_B[i] + sum(A[i, :] .* fixed_B)  # summing over j≠i (diagonal of A assumed to be 0 here)
        end
        # Define intrinsic growth rates: r[i] = k[i] (and assume D[i] = 1)
        r = k
        
        #############################
        # Step 4: Simulate the Full Model
        #############################
        # Our GLV model becomes:
        #   dB[i]/dt = B[i]*( r[i] - B[i] + sum_{j ≠ i} A[i,j]*B[j] )
        u0 = fixed_B  # Start at the fixed abundances (should be equilibrium if the system is stable)
        p = (r, A)
        prob = ODEProblem(full_model!, u0, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
        B_eq_full = sol.u[end]
        
        # Persistence: fraction of species with B > threshold
        thresh = 1e-3
        φ = sum(B_eq_full .> thresh) / n
        
        #############################
        # Step 5: Perturbation Experiment and Return Time Measurement
        #############################
        rt_full = simulate_perturbation(u0, p, tspan, t_perturb, perturb_factor)
        mean_rt_full = mean(skipmissing(rt_full))
        
        #############################
        # Step 6: Construct the Simplified Model
        #############################
        # Remove detailed trophic structure by replacing all off-diagonal entries with their mean.
        offdiag_vals = [A[i,j] for i in 1:n for j in 1:n if i != j]
        avg_off = mean(offdiag_vals)
        A_simpl = zeros(n, n)
        for i in 1:n
            for j in 1:n
                if i == j
                    A_simpl[i,j] = -1.0
                else
                    A_simpl[i,j] = avg_off
                end
            end
        end
        p_simpl = (r, A_simpl)
        prob_simpl = ODEProblem(full_model!, u0, tspan, p_simpl)
        sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
        B_eq_simpl = sol_simpl.u[end]
        rt_simpl = simulate_perturbation(u0, p_simpl, tspan, t_perturb, perturb_factor)
        mean_rt_simpl = mean(skipmissing(rt_simpl))
        
        #############################
        # Step 7: Abundance Distribution Metrics
        #############################
        # Since the fixed abundances are our “observed” equilibrium,
        # we fit them to a lognormal distribution.
        fit_ln = fit(LogNormal, fixed_B)
        # Compute the mean and variance of the fitted lognormal distribution:
        mean_ln = exp(fit_ln.μ + (fit_ln.σ^2)/2)
        var_ln = exp(2*fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2)-1)
        RelVar = var_ln / mean_ln^2  # Relative variance
        
        #############################
        # Step 8: Compute Emergent Network Parameters
        #############################
        # For the full interaction matrix A, we compute emergent parameters:
        # Effective interaction strength, variability, and reciprocity.
        μ_eff = n * mean(offdiag_vals)
        σ_eff = sqrt(n) * std(offdiag_vals)
        pairs = [(A[i,j], A[j,i]) for i in 1:n for j in 1:n if i < j]
        if length(pairs) > 1
            vals1 = [p[1] for p in pairs]
            vals2 = [p[2] for p in pairs]
            γ_eff = cor(vals1, vals2)
        else
            γ_eff = 0.0
        end
        
        #############################
        # Step 9: Predict Return Time Using Emergent & Abundance Metrics
        #############################
        T_pred = predicted_return_time(φ, γ_eff, RelVar)
        
        #############################
        # Save Metrics for the Current Iteration
        #############################
        push!(persistence_arr, φ)
        push!(rt_full_arr, mean_rt_full)
        push!(rt_simpl_arr, mean_rt_simpl)
        push!(pred_rt_arr, T_pred)
        push!(relvar_arr, RelVar)
    end
    
    results[n] = Dict(:persistence => persistence_arr,
                      :return_full => rt_full_arr,
                      :return_simpl => rt_simpl_arr,
                      :return_pred => pred_rt_arr,
                      :relvar => relvar_arr)
end

# ---------------------------
# Plotting Results
# ---------------------------
# Boxplots for persistence, return times, and relative variance across scenarios.

# Define numerical categories
category_ids = 1:length(species_scenarios)
category_labels = string.(species_scenarios)

# Helper function to flatten values and repeat category index
function flatten_data(data_by_group)
    values = reduce(vcat, data_by_group)
    categories = reduce(vcat, [fill(i, length(data_by_group[i])) for i in 1:length(data_by_group)])
    return categories, values
end

# === Prepare data ===
persistence_groups = [results[s][:persistence] for s in species_scenarios]
return_full_groups = [results[s][:return_full] for s in species_scenarios]
return_simpl_groups = [results[s][:return_simpl] for s in species_scenarios]
return_pred_groups = [results[s][:return_pred] for s in species_scenarios]
relvar_groups = [results[s][:relvar] for s in species_scenarios]

# === Plot 1: Persistence ===
fig1 = Figure()
ax1 = Axis(fig1[1, 1], xlabel = "Species Count", ylabel = "Persistence (φ)", 
           title = "Persistence Across Scenarios", xticks = (category_ids, category_labels))
cats, vals = flatten_data(persistence_groups)
MK.boxplot!(ax1, cats, vals)
display(fig1)

# === Plot 2: Return Times ===
fig2 = Figure(resolution = (1000, 300))

# Full Model
ax2 = Axis(fig2[1, 1], xlabel = "Species Count", ylabel = "Mean Return Time", 
           title = "Return Time (Full Model)", xticks = (category_ids, category_labels))
cats, vals = flatten_data(return_full_groups)
MK.boxplot!(ax2, cats, vals)

# Simplified Model
ax3 = Axis(fig2[1, 2], xlabel = "Species Count", ylabel = "Mean Return Time", 
           title = "Return Time (Simplified Model)", xticks = (category_ids, category_labels))
cats, vals = flatten_data(return_simpl_groups)
MK.boxplot!(ax3, cats, vals)

# Predicted
ax4 = Axis(fig2[1, 3], xlabel = "Species Count", ylabel = "Predicted Return Time", 
           title = "Predicted Return Time", xticks = (category_ids, category_labels))
cats, vals = flatten_data(return_pred_groups)
MK.boxplot!(ax4, cats, vals)

display(fig2)

# === Plot 3: Relative Variance ===
fig3 = Figure()
ax5 = Axis(fig3[1, 1], xlabel = "Species Count", ylabel = "Relative Variance", 
           title = "RelVar of Fixed Abundances", xticks = (category_ids, category_labels))
cats, vals = flatten_data(relvar_groups)
MK.boxplot!(ax5, cats, vals)
display(fig3)

