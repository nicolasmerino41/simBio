#############################
# Julia Pipeline: Fixed-Abundance Trophic Community Study
# Revised to address stiff dynamics (e.g., dt forced below epsilon)
#############################
# ---------------------------
# Function Definitions
# ---------------------------
# Full GLV predator-prey model.
# We assume D[i] = 1 so that r[i] = k[i].
# Model equation:
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
# Return times are measured as the time required for each species to recover within 10%
# of its pre-perturbation value.
function simulate_perturbation(u0, p, tspan, t_perturb, perturb_factor; solver=Rodas5())
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, solver, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]
    
    perturbed_state = pre_state .* perturb_factor
    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(full_model!, perturbed_state, tspan2, p)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)
    
    n = length(pre_state)
    return_times = zeros(n)
    for i in 1:n
        target = pre_state[i]
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
perturb_factor = 0.9              # Use a gentler perturbation factor (e.g., 0.8) to avoid extreme stiffness

# Container for results (per species scenario)
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
        try
            #############################
            # Step 1: Generate Fixed Equilibrium Abundances
            #############################
            # Assume the observed equilibrium abundances B* come from a lognormal distribution.
            μ_log_obs = 0.5
            σ_log_obs = 1.0
            fixed_B = rand(LogNormal(μ_log_obs, σ_log_obs), n)
            
            #############################
            # Step 2: Construct a Trophic (Predator–Prey) Network
            #############################
            connectance = 0.2
            A = zeros(n, n)
            for i in 1:n
                for j in 1:n
                    if i != j && rand() < connectance
                        if rand() < 0.5
                            β = rand(Exponential(1.0))
                            A[i, j] = β
                            A[j, i] = -0.3 * β
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
            # Force the equilibrium by setting:
            #   k_i = B_i* + sum_{j ≠ i} A[i,j] * B_j*
            k = zeros(n)
            for i in 1:n
                k[i] = fixed_B[i] + sum(A[i, :] .* fixed_B)
            end
            r = k  # Let r[i] = k[i], with D[i]=1
            
            #############################
            # Step 4: Simulate the Full Model
            #############################
            u0 = fixed_B  # Initial conditions are the fixed abundances.
            p = (r, A)
            prob = ODEProblem(full_model!, u0, tspan, p)
            sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)
            B_eq_full = sol.u[end]
            
            thresh = 1e-3
            φ = sum(B_eq_full .> thresh) / n
            
            #############################
            # Step 5: Perturbation Experiment and Return Time
            #############################
            rt_full = simulate_perturbation(u0, p, tspan, t_perturb, perturb_factor; solver=Rodas5())
            mean_rt_full = mean(skipmissing(rt_full))
            
            #############################
            # Step 6: Simplified Model
            #############################
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
            sol_simpl = solve(prob_simpl, Rodas5(), reltol=1e-8, abstol=1e-8)
            B_eq_simpl = sol_simpl.u[end]
            rt_simpl = simulate_perturbation(u0, p_simpl, tspan, t_perturb, perturb_factor; solver=Rodas5())
            mean_rt_simpl = mean(skipmissing(rt_simpl))
            
            #############################
            # Step 7: Abundance Distribution Metrics
            #############################
            fit_ln = fit(LogNormal, fixed_B)
            mean_ln = exp(fit_ln.μ + (fit_ln.σ^2)/2)
            var_ln = exp(2*fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2)-1)
            RelVar = var_ln / mean_ln^2
            
            #############################
            # Step 8: Emergent Network Parameters
            #############################
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
            # Step 9: Predict Return Time
            #############################
            T_pred = predicted_return_time(φ, γ_eff, RelVar)
            
            #############################
            # Save iteration metrics
            #############################
            push!(persistence_arr, φ)
            push!(rt_full_arr, mean_rt_full)
            push!(rt_simpl_arr, mean_rt_simpl)
            push!(pred_rt_arr, T_pred)
            push!(relvar_arr, RelVar)
        catch e
            @warn "Iteration skipped due to error" exception = e
        end
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
begin

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
end
# === Plot 3: Relative Variance ===
fig3 = Figure()
ax5 = Axis(fig3[1, 1], xlabel = "Species Count", ylabel = "Relative Variance", 
           title = "RelVar of Fixed Abundances", xticks = (category_ids, category_labels))
cats, vals = flatten_data(relvar_groups)
MK.boxplot!(ax5, cats, vals)
display(fig3)

