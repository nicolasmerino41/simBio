#############################
# Julia Pipeline: Effective Model Study with Abundance Distribution Integration
# Scenarios with 20, 40, and 60 species; multiple iterations to capture randomness.
# We simulate a full GLV model, extract community metrics, conduct perturbation experiments,
# build a simplified version of the model, and compute analytical (predicted) return times using
# parameters from a lognormal fit to the survivors’ abundances.
#############################
# --- Define our basic functions ---

# Full generalized Lotka–Volterra ODE model.
function full_model!(du, u, p, t)
    r, A = p
    n = length(u)
    for i in 1:n
        du[i] = u[i] * (r[i] + sum(A[i, :] .* u))
    end
end

# Perturbation experiment: integrate until perturbation, then apply a factor, and measure return time.
function simulate_perturbation(u0, p, tspan, perturb_time, perturb_factor)
    # Integrate up to the perturbation
    tspan1 = (tspan[1], perturb_time)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8)
    state_pre = sol1.u[end]
    # Apply perturbation: multiply abundances by perturb_factor
    state_perturbed = state_pre .* perturb_factor
    # Integrate from perturbation time to end
    tspan2 = (perturb_time, tspan[2])
    prob2 = ODEProblem(full_model!, state_perturbed, tspan2, p)
    sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
    n = length(state_pre)
    return_times = zeros(n)
    for i in 1:n
        target = state_pre[i]
        recovered = false
        for (t, state) in zip(sol2.t, sol2.u)
            # Recovery criterion: within 10% of pre-perturbation value
            if abs(state[i] - target) / (abs(target)+1e-8) < 0.1
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

# A simple hypothetical prediction function for return time:
# Here, T_pred is assumed to be inversely related to 1 - φ * γ * (relative variance)
function predicted_return_time(φ, γ, rel_var)
    denom = 1 - φ * γ * rel_var
    return denom > 0 ? 1 / denom : Inf
end

# --- Simulation parameters --- 
species_scenarios = [20, 40, 60]
Niter = 100         # number of iterations per species count
tspan = (0.0, 100.0)
perturb_time = 50.0
perturb_factor = 0.5  # reduce abundances by 50% at perturbation

# Containers to store results for each species scenario.
results = Dict{Int,Dict{Symbol, Vector{Float64}}}()

for n in species_scenarios
    # Create arrays for each metric for the current species count.
    persistence_full_arr = Float64[]
    return_full_arr = Float64[]
    return_simpl_arr = Float64[]
    pred_return_arr = Float64[]
    rel_variance_arr = Float64[]
    
    for iter in 1:Niter
        # --- Step 1: Generate Full Model Parameters ---
        # r is drawn from a LogNormal to promote a fat-tailed abundance distribution.
        μ_r = 0.0
        σ_r = 0.5
        r_full = rand(LogNormal(μ_r, σ_r), n)
        D_full = ones(n)
        # Build interaction matrix A_full:
        A_full = zeros(n, n)
        for i in 1:n
            for j in 1:n
                if i == j
                    A_full[i,j] = -1.0
                else
                    A_full[i,j] = rand(Normal(0, 0.1))
                end
            end
        end
        
        # --- Step 2: Solve Full Model ---
        u0_full = 0.5 .* r_full  # initial conditions
        p_full = (r_full, A_full)
        prob_full = ODEProblem(full_model!, u0_full, tspan, p_full)
        sol_full = solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)
        equilibrium_full = sol_full.u[end]
        
        # Define surviving species (abundance above a threshold)
        threshold = 1e-3
        survivors = equilibrium_full .> threshold
        φ_full = sum(survivors) / n
        
        # --- Step 3: Fit Lognormal to Survivors ---
        if sum(survivors) > 0
            survivors_abund = equilibrium_full[survivors]
            lognormal_fit = fit(LogNormal, survivors_abund)
            fitted_mean = exp(lognormal_fit.μ + (lognormal_fit.σ^2)/2)
            fitted_var = exp(2*lognormal_fit.μ + lognormal_fit.σ^2) * (exp(lognormal_fit.σ^2)-1)
            rel_var = fitted_var / fitted_mean^2  # relative variance (coefficient of variation^2)
        else
            rel_var = NaN
        end
        
        # --- Step 4: Perturbation Experiments ---
        rt_full = simulate_perturbation(u0_full, p_full, tspan, perturb_time, perturb_factor)
        mean_rt_full = mean(skipmissing(rt_full))
        
        # --- Step 5: Construct the Simplified Model ---
        # Replace all off-diagonals with their overall mean.
        offdiag_vals = [] 
        for i in 1:n
            for j in 1:n
                if i != j
                    push!(offdiag_vals, A_full[i,j])
                end
            end
        end
        
        avg_offdiag = mean(offdiag_vals)
        A_simpl = zeros(n, n)
        for i in 1:n
            for j in 1:n
                if i == j
                    A_simpl[i,j] = -1.0
                else
                    A_simpl[i,j] = avg_offdiag
                end
            end
        end
        p_simpl = (r_full, A_simpl)
        u0_simpl = u0_full
        prob_simpl = ODEProblem(full_model!, u0_simpl, tspan, p_simpl)
        sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
        equilibrium_simpl = sol_simpl.u[end]
        rt_simpl = simulate_perturbation(u0_simpl, p_simpl, tspan, perturb_time, perturb_factor)
        mean_rt_simpl = mean(skipmissing(rt_simpl))
        
        # --- Step 6: Compute Emergent Parameters from Full Model ---
        # Effective carrying capacities K (here simply r since D=1).
        K_full = r_full
        ζ = std(K_full)
        # Off-diagonals of A.
        offdiag_full = [A_full[i,j] for i in 1:n for j in 1:n if i != j]
        μ_eff = n * mean(offdiag_full)
        σ_eff = sqrt(n) * std(offdiag_full)
        # Reciprocity γ_eff computed over pairs (i<j)
        pairs = [(A_full[i,j], A_full[j,i]) for i in 1:n for j in 1:n if i < j]
        vals1 = [p[1] for p in pairs]
        vals2 = [p[2] for p in pairs]
        γ_eff = length(vals1) > 1 ? cor(vals1, vals2) : 0.0
        
        # --- Step 7: Compute Predicted Return Time ---
        T_pred = predicted_return_time(φ_full, γ_eff, rel_var)
        
        # --- Save metrics for this iteration ---
        push!(persistence_full_arr, φ_full)
        push!(return_full_arr, mean_rt_full)
        push!(return_simpl_arr, mean_rt_simpl)
        push!(pred_return_arr, T_pred)
        push!(rel_variance_arr, rel_var)
    end
    
    results[n] = Dict(:persistence => persistence_full_arr,
                       :return_full => return_full_arr,
                       :return_simpl => return_simpl_arr,
                       :return_pred => pred_return_arr,
                       :rel_var => rel_variance_arr)
end

# ========= Step 8: Plotting the Results ==========
# Create box plots summarizing each attribute (persistence, return time, relative variance) across scenarios.
# Persistence Plot:
begin
    # Prepare data
    species_labels = reduce(vcat, [fill(s, length(results[s][:persistence])) for s in species_scenarios])

    persistence_data = reduce(vcat, [results[s][:persistence] for s in species_scenarios])
    return_full_data = reduce(vcat, [results[s][:return_full] for s in species_scenarios])
    return_simpl_data = reduce(vcat, [results[s][:return_simpl] for s in species_scenarios])
    return_pred_data = reduce(vcat, [results[s][:return_pred] for s in species_scenarios])
    rel_var_data = reduce(vcat, [results[s][:rel_var] for s in species_scenarios])

    # Plot 1: Persistence
    fig1 = Figure()
    ax1 = Axis(fig1[1, 1], xlabel="Species Count", ylabel="Persistence (φ)", title="Persistence Comparison")
    MK.boxplot!(ax1, species_labels, persistence_data)
    display(fig1)

    # Plot 2-4: Return Times
    fig2 = Figure(resolution=(1000, 300))
    
    ##### AX2: Full Model Return Time #####
    ax2 = Axis(fig2[1, 1], xlabel="Species Count", ylabel="Mean Return Time", title="Return Time (Full Model)")
    # Remove NaN indexes
    valid_indexes = .!isnan.(return_full_data)
    species_labels_ = species_labels[valid_indexes]
    return_full_data_ = return_full_data[valid_indexes]
    MK.boxplot!(ax2, species_labels_, return_full_data_)

    ##### AX3: Simplified Model Return Time #####
    # Remove NaN indexes
    valid_indexes = .!isnan.(return_simpl_data)
    species_labels_ = species_labels[valid_indexes]
    return_simpl_data_ = return_simpl_data[valid_indexes]
    ax3 = Axis(fig2[1, 2], xlabel="Species Count", ylabel="Mean Return Time", title="Return Time (Simplified Model)")
    MK.boxplot!(ax3, species_labels_, return_simpl_data_)

    ##### AX4: Predicted Return Time #####
    # Remove NaN indexes
    valid_indexes = .!isnan.(return_pred_data)
    species_labels_ = species_labels[valid_indexes]
    return_pred_data_ = return_pred_data[valid_indexes]
    ax4 = Axis(fig2[1, 3], xlabel="Species Count", ylabel="Predicted Return Time", title="Predicted Return Time")
    MK.boxplot!(ax4, species_labels_, return_pred_data_)

    display(fig2)

    # Plot 5: Relative Variance
    # Prepare data
    valid_indexes = .!isnan.(rel_var_data)
    species_labels_ = species_labels[valid_indexes]
    rel_var_data_ = rel_var_data[valid_indexes]
    fig3 = Figure()
    ax5 = Axis(fig3[1, 1], xlabel="Species Count", ylabel="Relative Variance", title="Abundance Distribution: Relative Variance")
    MK.boxplot!(ax5, species_labels_, rel_var_data_)
    display(fig3)
end
