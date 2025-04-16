# Revised Jacobian computation
function compute_jacobian(B, p)
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
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
            J[i, j] = B[i] * α * epsilon * A[i, j]
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
function simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=false, show_warnings=true)
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
    
    if plot
        fig = Figure(; size =(1600,800))
        ax1 = Axis(fig[1,1], xlabel="Time", ylabel="Abundance",
                  title="Community Response Before Press Perturbation")
        ax2 = Axis(fig[1,2], xlabel="Time", ylabel="Abundance",
                    title="Community Response After Press Perturbation")
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

    return return_times, new_equil
end

# ---------------------------
# Master pipeline function.
"""
    run_simulation_pipeline(R_vals, C_vals, conn_vals, delta_vals; Niter=10, tspan=(0.0,50.0), t_perturb=25.0, plot=false)

For each combination of resource count, consumer count, connectance, and press perturbation intensity,
this function:
  - Generates a target equilibrium state from a LogNormal.
  - Constructs a global interaction matrix A (with resources not preying).
  - Defines model parameters and calibrates them (using calibrate_params, assumed defined).
  - Simulates the full model to obtain equilibrium, computes persistence, resilience, and reactivity.
  - Simulates a press perturbation for both the full and a simplified model.
  - Returns a DataFrame of results.
"""
function run_simulation_pipeline(
    R_vals::Vector{Int}, C_vals::Vector{Int},
    conn_vals::AbstractRange{Float64}, delta_vals::AbstractRange{Float64};
    Niter::Int=10,
    tspan::Tuple{Float64,Float64}=(0.0,50.0),
    t_perturb::Float64=25.0,
    plot::Bool=false, 
    show_warnings::Bool=true
)
    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()
    for R in R_vals
        for C in C_vals
            total_species = R + C
            for conn in conn_vals
                for delta in delta_vals
                    for iter in 1:Niter
                        try
                            # Step 1: Generate fixed equilibrium abundances.
                            fixed_state = abs.(rand(LogNormal(5.0, 10.0), total_species))
                            
                            # Step 2: Construct global interaction matrix A.
                            A = zeros(total_species, total_species)
                            for i in (R+1):total_species
                                for j in 1:total_species
                                    if i != j && rand() < conn && iszero(A[i,j])
                                        A[i,j] = 1.0
                                        A[j,i] = -1.0
                                    end
                                end
                            end
                            
                            # Step 3: Define remaining parameters.
                            m_cons = rand(Uniform(0.1, 1.0), C)    # Mortality rates for consumers.
                            xi_cons = ones(C)                      # Initial consumer thresholds.
                            r_res = ones(R)                        # Intrinsic growth rates for resources.
                            d_res = ones(R)                        # Scaling factors for resources.
                            epsilon = 0.5
                            p_tuple = (R, C, m_cons, d_res, epsilon, A)
                            
                            # Step 4: Calibrate thresholds so that fixed_state is an equilibrium.
                            R_eq = fixed_state[1:R]
                            C_eq = fixed_state[R+1:total_species]
                            new_xi_cons, new_r_res, A = calibrate_params(R_eq, C_eq, p_tuple)
                            count_calib = 0
                            max_calib = 10
                            while (any(isnan, new_xi_cons) || any(isnan, new_r_res)) && (count_calib < max_calib)
                                fixed_state = abs.(rand(Normal(0.5, 1.0), total_species))
                                R_eq = fixed_state[1:R]
                                C_eq = fixed_state[R+1:total_species]
                                new_xi_cons, new_r_res, A = calibrate_params(R_eq, C_eq, p_tuple)
                                count_calib += 1
                            end

                            if any(isnan, new_xi_cons) || any(isnan, new_r_res)
                                error("Calibration failed after $max_calib attempts; skipping this iteration.")
                            end

                            # Update the parameter tuple with the calibrated thresholds and growth rates.
                            p_tuple = (R, C, m_cons, new_xi_cons, new_r_res, d_res, epsilon, A)
                            println("New consumer thresholds (xi_cons): ", new_xi_cons)
                            println("New resource growth rates (r_res): ", new_r_res)
                            
                            # Step 5: Simulate full model.
                            prob_full = ODEProblem(trophic_ode!, fixed_state, tspan, p_tuple)
                            sol_full = solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)
                            B_eq_full = sol_full.u[end]
                            
                            if sol_full.t[end] < t_perturb || any(isnan, sol_full.u[end]) || any(isinf, sol_full.u[end])
                                error("Error: solution did not finish properly")
                            end

                            persistence = sum(x -> x > EXTINCTION_THRESHOLD, B_eq_full) / total_species
                            resilience = compute_resilience(B_eq_full, p_tuple)
                            reactivity = compute_reactivity(B_eq_full, p_tuple)
                            
                            # Step 6: Full model press perturbation.
                            rt_full, new_equil_full = simulate_press_perturbation(fixed_state, p_tuple, tspan, t_perturb, delta;
                                                                solver=Tsit5(), plot=plot)
                            mean_rt_full = mean(filter(!isnan, rt_full))
                            after_persitence_full = sum(x -> x > EXTINCTION_THRESHOLD, new_equil_full) / total_species
                            
                            # Step 7: Simplified model version.
                            offdiag_values = [A[i,j] for i in 1:total_species, j in 1:total_species if i != j && A[i,j] > 0.0]
                            avg_off = isempty(offdiag_values) ? 0.0 : mean(offdiag_values)
                            A_simpl = copy(A)
                            A_simpl[A .> 0.0] .= avg_off
                            A_simpl[A .< 0.0] .= -avg_off*epsilon
                            p_simpl = (R, C, m_cons, new_xi_cons, new_r_res, d_res, epsilon, A_simpl)
                            rt_simpl, new_equil_simpl = simulate_press_perturbation(fixed_state, p_simpl, tspan, t_perturb, delta;
                                                                    solver=Tsit5(), plot=plot)
                            mean_rt_simpl = mean(filter(!isnan, rt_simpl))
                            after_persitence_simpl = sum(x -> x > EXTINCTION_THRESHOLD, new_equil_simpl) / total_species
                            # Step 8: Additional metrics: e.g., relative variance from LogNormal fit.
                            fit_ln = fit(LogNormal, fixed_state)
                            mean_ln = exp(fit_ln.μ + (fit_ln.σ^2)/2)
                            var_ln = exp(2 * fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2) - 1)
                            relVar = var_ln / mean_ln^2
                            
                            rec = (; 
                                fully_stable = all(B_eq_full .> EXTINCTION_THRESHOLD),
                                resource_count = R,
                                consumer_count = C,
                                total_species = total_species,
                                connectance = conn,
                                perturbation_delta = delta,
                                iteration = iter,
                                persistence = persistence,
                                after_persistence_full = after_persitence_full,
                                after_persistence_simpl = after_persitence_simpl,
                                return_time_full = mean_rt_full,
                                return_time_simplified = mean_rt_simpl,
                                resilience = resilience,
                                reactivity = reactivity,
                                relVar = relVar
                            )

                            # Step 9: Store results in a thread-safe manner.
                            lock(results_lock) do
                                push!(results, rec)
                            end
                        catch e
                            # Handle exceptions gracefully, e.g., log or print a message.
                            println("Error in iteration $iter for R=$R, C=$C, conn=$conn, delta=$delta: $e")
                            continue  # Skip to the next iteration.
                        end
                    end
                end
            end
        end
    end
    return DataFrame(results)
end

###############################################################################
# 3. Run Simulations and Gather Results
###############################################################################
# Specify parameter ranges:
R_vals = [20]            # e.g., try 4 and 6 resources
C_vals = [5]            # e.g., try 5 and 8 consumers
connectance_list = 0.1:0.05:0.5
delta_list = 0.01:0.01:0.02

# Run the simulation pipeline.
df_results = run_simulation_pipeline(R_vals, C_vals, connectance_list, delta_list;
                                     Niter=10, tspan=(0.0,500.0), t_perturb=250.0, plot=false)
println("First 10 rows of simulation results:")
println(first(df_results, 10))

###############################################################################
# 4. Post-processing Analysis and Plotting
###############################################################################
using GLM
# (A) Scatter plot of resilience (predictor) vs. full-model mean return time.
begin
    fig1 = Figure(; size=(1200,600))
    ax1 = Axis(fig1[1,1], xlabel="Resilience", ylabel="Mean Return Time (Full Model)",
            title="Return Time vs. Resilience")
    # df_filtered = filter(row -> row.resilience > -2.0 && row.resilience < 0.0, df_results)
    df_filtered = filter(row -> !iszero(row.return_time_full), df_results)
    df_filtered = filter(row -> row.fully_stable, df_filtered)
    MK.scatter!(ax1, df_filtered.resilience, df_filtered.return_time_full, markersize=8, color=:blue)
    
    # Fit a linear model for return_time_full ~ resilience.
    model1 = lm(@formula(return_time_full ~ resilience), df_filtered)
    res_lin = sort(df_filtered.resilience)
    predicted = predict(model1)
    # For plotting a regression line, we sort the data.
    sorted_idx = sortperm(df_filtered.resilience)
    MK.lines!(ax1, df_filtered.resilience[sorted_idx], predicted[sorted_idx], color=:red, linewidth=2)
    # Label(fig1, "Pearson cor = $(round(cor(df_filtered.resilience, df_filtered.return_time_full), digits=2))", position=(20,20))
    display(fig1)
end
begin
    # (B) Scatter plot comparing full model vs. simplified model return times.
    fig2 = Figure(; size=(1200,600))
    ax2 = Axis(fig2[1,1], xlabel="Return Time (Full Model)", ylabel="Return Time (Simplified Model)",
            title="Full vs. Simplified Return Times")
    MK.scatter!(ax2, df_filtered.return_time_full, df_filtered.return_time_simplified, markersize=8, color=:green)
    # Plot a 1:1 line.
    min_val = minimum(vcat(df_filtered.return_time_full, df_filtered.return_time_simplified))
    max_val = maximum(vcat(df_filtered.return_time_full, df_filtered.return_time_simplified))
    lines!(ax2, [min_val, max_val], [min_val, max_val], color=:black, linestyle=:dash, linewidth=2)
    display(fig2)
end

begin
    # (C) Analyze persistence and relative variance vs. connectance.
    fig3 = Figure(; size=(1500,600))
    ax3 = Axis(fig3[1,1], xlabel="Connectance", ylabel="Persistence (Initial)",
            title="Initial Persistence vs. Connectance")
    ax4 = Axis(fig3[1,2], xlabel="Connectance", ylabel="Persistence Full(Post-Perturbation)",
            title="Post-Perturbation Persistence Full vs. Connectance")
    ax5 = Axis(fig3[1,3], xlabel="Connectance", ylabel="Persistence Simplified(Post-Perturbation)",
            title="Post-Perturbation Persistence Simplified vs. Connectance")
    ax6 = Axis(fig3[1,4], xlabel="Persistence Full(Post-Perturbation)", ylabel="Persistence Simplified(Post-Perturbation)")
    
    MK.scatter!(ax3, df_filtered.connectance, df_filtered.persistence, markersize=8, color=:purple)
    MK.scatter!(ax4, df_filtered.connectance, df_filtered.after_persistence_full, markersize=8, color=:orange)
    MK.scatter!(ax5, df_filtered.connectance, df_filtered.after_persistence_simpl, markersize=8, color=:brown)
    MK.scatter!(ax6, df_filtered.after_persistence_full, df_filtered.after_persistence_simpl, markersize=8, color=:pink)

    display(fig3)
end

# Print simple correlation statistics:
println("\nCorrelation between resilience and full-model return time: ",
        round(cor(df_filtered.resilience, df_filtered.return_time_full), digits=3))
println("Correlation between full and simplified return times: ",
        round(cor(df_filtered.return_time_full, df_filtered.return_time_simplified), digits=3))
println("Correlation between relative variance and persistence: ",
        round(cor(df_results.relVar, df_results.persistence), digits=3))