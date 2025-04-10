#############################
# Julia Pipeline: Fixed-Abundance Trophic Community Study with Press Perturbation
#############################
# ---------------------------
# Function Definitions
# ---------------------------

# Full GLV predator-prey model.
# We assume D[i] = 1 so that r[i] = k[i].
# Model:  dB[i]/dt = B[i] * (r[i] - B[i] + sum_{j≠i} A[i,j] * B[j])
function full_model!(du, B, p, t)
    r, A = p
    n = length(B)
    for i in 1:n
        du[i] = B[i] * (r[i] - B[i] + sum(A[i, :] .* B))
    end
end

# Function to simulate a "press perturbation" by changing carrying capacities.
# At t_perturb, r is updated to r_press = (1 - δ)*r.
# The new equilibrium (target) is approximated by integrating until t_end with the new r.
# The return time for species i is defined as the first time after t_perturb when its abundance
# is within 10% of the new equilibrium abundance.
function simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=false)
    # Integrate with original parameters to t_perturb.
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, solver, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]

    # Modify the intrinsic growth rates (carrying capacities) for a press perturbation.
    r, A = p
    r_press = (1 - delta) .* r
    p_press = (r_press, A)

    # Integrate the system with new parameters from t_perturb to end.
    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(full_model!, pre_state, tspan2, p_press)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)

    # Approximate the new equilibrium as the final state.
    new_equil = sol2.u[end]
    n = length(new_equil)
    return_times = zeros(n)
    # For each species, find the first time after t_perturb with error < 10%.
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
    if plot
        fig = Figure(; size=(800, 600))
        ax = Axis(fig[1, 1], xlabel="Time", ylabel="Abundance", title="Community Response to a Perturbation")
        # Plot time series for each species 
        for i in 1:n
            # Extract the time series for species i across all time steps.
            lines!(ax, sol2.t, sol2[i, :], label="Species $i")
        end
        display(fig)
    end
    return return_times
end

# A hypothetical predicted return time function.
# Here, the prediction combines persistence, network reciprocity (γ_eff) and the abundance distribution heterogeneity (RelVar).
function predicted_return_time(φ, γ_eff, relVar)
    # For press perturbations with fixed abundances (φ may be 1), the formula simplifies to:
    # T_pred ≈ 1 / (1 - γ_eff * relVar)
    denom = 1 - φ * γ_eff * relVar 
    return denom > 0 ? 1 / denom : Inf
end

# ---------------------------
# Pipeline Parameters
# ---------------------------
species_scenarios = [10, 20, 30]   # Different number of species
Niter = 1                      # Iterations per scenario
tspan = (0.0, 100.0)               # Total integration time
t_perturb = 50.0                   # Time of press perturbation
delta = 0.2                        # Relative reduction in carrying capacity (20% reduction)

# Containers for results (per species scenario)
results_shit = Dict{Int,Dict{Symbol,Vector{Float64}}}()

# ---------------------------
# Main Loop: Over Species Scenarios and Iterations
# ---------------------------
# lock_result = ReentrantLock()
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
            # These fixed abundances represent the empirical data.
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
                            A[i, j] = 0.3 * β     # i preys on j
                            A[j, i] = -β     # effect on prey
                        else
                            β = rand(Exponential(1.0))
                            A[j, i] = 0.3 * β # j preys on i
                            A[i, j] = -β     # effect on prey
                        end
                    end
                end
            end

            #############################
            # Step 3: Calibrate Model with Fixed Abundances
            #############################
            # For each species i, force the equilibrium by adjusting the carrying capacity:
            #   k_i = B_i* + sum_{j ≠ i} A[i,j] * B_j*
            k = zeros(n)
            for i in 1:n
                # Calculate the interaction sum, excluding the diagonal (since A[i,i] is assumed 0)
                s = 0.0
                for j in 1:n
                    if i != j
                        s += A[i, j] * fixed_B[j]
                    end
                end
                k[i] = fixed_B[i] - s
            end
            r = k

            #############################
            # Step 4: Simulate the Full Model under Fixed Conditions
            #############################
            # Our GLV model becomes:
            #   dB[i]/dt = B[i]*(r[i] - B[i] + sum_{j ≠ i} A[i,j]*B[j])
            u0 = fixed_B  # Start exactly at the fixed abundances
            p = (r, A)
            prob = ODEProblem(full_model!, u0, tspan, p)
            sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
            B_eq_full = sol[:, end]
            thresh = 1e-3
            φ = sum(B_eq_full .> thresh) / n   # With fixed abundances, φ should be 1 (or near 1)
            # println("φ = $φ")

            if true
                fig = Figure(; size=(800, 600))
                ax = Axis(fig[1, 1], xlabel="Time", ylabel="Abundance", title="Before a Perturbation Full Model")
                # Plot time series for each species
                for i in 1:n
                    # Extract the time series for species i across all time steps.
                    species_ts = sol[i, :]
                    lines!(ax, sol.t, species_ts, label="Species $i")
                end
                display(fig)
            end
            #############################
            # Step 5: Press Perturbation via Carrying Capacity Change (Full Model)
            #############################
            rt_full = simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=true)
            mean_rt_full = mean(skipmissing(rt_full))

            #############################
            # Step 6: Construct the Simplified Model
            #############################
            offdiag_vals = [A[i, j] for i in 1:n for j in 1:n if i != j && A[i, j] < 0.0]
            avg_off = abs(mean(offdiag_vals))
            A_simpl = zeros(n, n)

            for i in 1:n
                for j in 1:n
                    if i == j
                        A_simpl[i, j] = 0.0
                    elseif i != j && A[i, j] > 0.0
                        A_simpl[i, j] = avg_off*0.3
                    elseif i != j && A[i, j] < 0.0
                        A_simpl[i, j] = -avg_off
                    end
                end
            end
            p_simpl = (r, A_simpl)
            prob_simpl = ODEProblem(full_model!, u0, tspan, p_simpl)
            sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
            B_eq_simpl = sol_simpl.u[end]

            if true
                fig = Figure(; size=(800, 600))
                ax = Axis(fig[1, 1], xlabel="Time", ylabel="Abundance", title="Before a Perturbation Simplified Model")
                # Plot time series for each species
                for i in 1:n
                    # Extract the time series for species i across all time steps.
                    species_ts = sol_simpl[i, :]
                    lines!(ax, sol_simpl.t, species_ts, label="Species $i")
                end
                display(fig)
            end

            rt_simpl = simulate_press_perturbation(u0, p_simpl, tspan, t_perturb, delta; solver=Tsit5(), plot=true)
            mean_rt_simpl = mean(skipmissing(rt_simpl))

            #############################
            # Step 7: Abundance Distribution Metrics
            #############################
            fit_ln = fit(LogNormal, fixed_B)
            mean_ln = exp(fit_ln.μ + (fit_ln.σ^2) / 2)
            var_ln = exp(2 * fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2) - 1)
            RelVar = var_ln / mean_ln^2

            #############################
            # Step 8: Emergent Network Parameters
            #############################
            μ_eff = n * mean(offdiag_vals)
            σ_eff = sqrt(n) * std(offdiag_vals)
            pairs = [(A[i, j], A[j, i]) for i in 1:n for j in 1:n if i < j]
            if length(pairs) > 1
                vals1 = [p[1] for p in pairs]
                vals2 = [p[2] for p in pairs]
                γ_eff = cor(vals1, vals2)
            else
                γ_eff = 0.0
            end

            #############################
            # Step 9: Predict Return Time Using the Emergent & Abundance Metrics
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
        catch e
            @warn "Iteration skipped due to error" exception = e
        end
    end
    results_shit[n] = Dict(
        :persistence => persistence_arr,
        :return_full => rt_full_arr,
        :return_simpl => rt_simpl_arr,
        :return_pred => pred_rt_arr,
        :relvar => relvar_arr
    )
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
persistence_groups = [results_shit[s][:persistence] for s in species_scenarios]
return_full_groups = [results_shit[s][:return_full] for s in species_scenarios]
return_simpl_groups = [results_shit[s][:return_simpl] for s in species_scenarios]
return_pred_groups = [results_shit[s][:return_pred] for s in species_scenarios]
relvar_groups = [results_shit[s][:relvar] for s in species_scenarios]
# === Plot 1: Persistence ===
begin
    fig1 = Figure()
    ax1 = Axis(fig1[1, 1], xlabel="Species Count", ylabel="Persistence (φ)",
        title="Persistence Across Scenarios", xticks=(category_ids, category_labels))
    cats, vals = flatten_data(persistence_groups)
    MK.boxplot!(ax1, cats, vals)
    display(fig1)
end
# === Plot 2: Return Times ===
begin
    fig2 = Figure(resolution=(1000, 300))

    # Full Model
    ax2 = Axis(
        fig2[1, 1], xlabel="Species Count", ylabel="Mean Return Time",
        title="Return Time (Full Model)", xticks=(category_ids, category_labels),
        # yticks = 0:0.1:1.0
    )
    cats, vals = flatten_data(return_full_groups)
    MK.boxplot!(ax2, cats, vals)

    # Simplified Model
    ax3 = Axis(
        fig2[1, 2], xlabel="Species Count", ylabel="Mean Return Time",
        title="Return Time (Simplified Model)", xticks=(category_ids, category_labels),
        # yticks = 0:0.1:1.0
    )
    cats, vals = flatten_data(return_simpl_groups)
    MK.boxplot!(ax3, cats, vals)

    # Predicted
    ax4 = Axis(
        fig2[1, 3], xlabel="Species Count", ylabel="Predicted Return Time",
        title="Predicted Return Time", xticks=(category_ids, category_labels),
        # yticks = 0:0.1:1.0
    )
    cats, vals = flatten_data(return_pred_groups)
    MK.boxplot!(ax4, cats, vals)

    display(fig2)
end

# === Plot 4: Relative Variance ===
begin
    fig3 = Figure()
    ax5 = Axis(fig3[1, 1], xlabel="Species Count", ylabel="Relative Variance",
        title="RelVar of Fixed Abundances", xticks=(category_ids, category_labels))
    cats, vals = flatten_data(relvar_groups)
    MK.boxplot!(ax5, cats, vals)
    display(fig3)
end

begin
    fig2 = Figure(resolution=(1000, 300))

    function mean_sd_plot!(ax, data_groups, title)
        means = [mean(g) for g in data_groups]
        sds = [std(g) for g in data_groups]
        x = 1:length(data_groups)

        # Plot error bars
        # errorbars!(ax, x, means, sds; linewidth=2)
        # Plot means as dots
        scatter!(ax, x, means; color=:black, markersize=8)

        ax.xticks = (x, category_labels)
        ax.title = title
        ax.xlabel = "Species Count"
        ax.ylabel = "Mean Return Time"
        # ax.yticks = 0:0.1:3.0
    end

    # Full Model
    ax2 = Axis(fig2[1, 1])
    mean_sd_plot!(ax2, return_full_groups, "Return Time (Full Model)")

    # Simplified Model
    ax3 = Axis(fig2[1, 2])
    mean_sd_plot!(ax3, return_simpl_groups, "Return Time (Simplified Model)")

    # Predicted Return Time
    ax4 = Axis(fig2[1, 3])
    mean_sd_plot!(ax4, return_pred_groups, "Predicted Return Time")

    display(fig2)
end