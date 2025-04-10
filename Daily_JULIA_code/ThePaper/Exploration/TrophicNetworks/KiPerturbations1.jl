#############################
# Function Definitions
#############################
# Full GLV predator-prey model.
function full_model!(du, B, p, t)
    r, A = p
    n = length(B)
    for i in 1:n
        du[i] = B[i] * (r[i] - B[i] + sum(A[i, :] .* B))
    end
end

# Simulate a press perturbation.
function simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=false)
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, solver, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]

    # Update intrinsic growth rates.
    r, A = p
    r_press = (1 - delta) .* r
    p_press = (r_press, A)

    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(full_model!, pre_state, tspan2, p_press)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)
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
        fig = Figure(resolution = (800, 600))
        ax = Axis(fig[1, 1], xlabel="Time", ylabel="Abundance", title="Community Response")
        for i in 1:n
            lines!(ax, sol2.t, sol2[i, :], label="Species $i")
        end
        display(fig)
    end
    return return_times
end

# Predicted return time function.
function predicted_return_time(φ, γ_eff, relVar)
    denom = 1 - φ * γ_eff * relVar
    return denom > 0 ? 1 / denom : Inf
end

# Compute the Jacobian matrix at equilibrium.
function compute_jacobian(B, A)
    n = length(B)
    J = zeros(n, n)
    for i in 1:n
        J[i, i] = -B[i]
        for j in 1:n
            if i != j
                J[i, j] = B[i] * A[i, j]
            end
        end
    end
    return J
end

# Resilience: negative of largest real eigenvalue.
function compute_resilience(B, A)
    J = compute_jacobian(B, A)
    ev = eigvals(J)
    return -maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of the symmetric part of J.
function compute_reactivity(B, A)
    J = compute_jacobian(B, A)
    Jsym = (J + J') / 2
    ev_sym = eigvals(Jsym)
    return maximum(real.(ev_sym))
end

#############################
# Pipeline Parameters
#############################
species_scenarios = [10, 20, 30]
Niter = 100
tspan = (0.0, 100.0)
t_perturb = 50.0
connectance_list = 0.05:0.05:0.25
delta_list = [0.2, 0.4, 0.6]

# We'll save all simulation results in an array of NamedTuples.
sim_results = Vector{NamedTuple{(:species_count, :delta, :iteration, :connectance, 
    :persistence, :return_full, :return_simpl, :return_pred, :relvar, :resilience, :reactivity), Tuple{Int, Float64, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()

# A lock for pushing safely from multiple threads.
results_lock = ReentrantLock()

#############################
# Main Simulation Loop (Threads)
#############################
for n in species_scenarios
    for conn in connectance_list
        for δ in delta_list
            Threads.@threads for iter = 1:Niter
                try
                    # Step 1: Generate fixed equilibrium abundances.
                    μ_log_obs = 0.5
                    σ_log_obs = 1.0
                    fixed_B = rand(LogNormal(μ_log_obs, σ_log_obs), n)
                    
                    # Step 2: Construct trophic network (using connectance).
                    A = zeros(n, n)
                    for i in 1:n
                        for j in 1:n
                            if i != j && rand() < conn
                                if rand() < 0.5
                                    β = rand(Exponential(1.0))
                                    A[i, j] = 0.3 * β   # i preys on j
                                    A[j, i] = -β        # effect on prey
                                else
                                    β = rand(Exponential(1.0))
                                    A[j, i] = 0.3 * β   # j preys on i
                                    A[i, j] = -β
                                end
                            end
                        end
                    end
                    
                    # Step 3: Calibrate model.
                    k = [ fixed_B[i] - sum(A[i, j] * fixed_B[j] for j in 1:n if j != i) for i in 1:n ]
                    r = k  # (D[i] = 1)
                    
                    # Step 4: Full model simulation.
                    u0 = fixed_B
                    p = (r, A)
                    prob = ODEProblem(full_model!, u0, tspan, p)
                    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
                    B_eq_full = sol[:, end]
                    thresh = 1e-3
                    φ = sum(B_eq_full .> thresh) / n
                    
                    # Step 5: Press perturbation (full).
                    rt_full = simulate_press_perturbation(u0, p, tspan, t_perturb, δ; solver=Tsit5(), plot=false)
                    mean_rt_full = mean(skipmissing(rt_full))
                    
                    # Step 6: Simplified model (using average off-diagonals).
                    offdiag_vals = [ A[i, j] for i in 1:n for j in 1:n if i != j ]
                    avg_off = mean(offdiag_vals)
                    A_simpl = [ i == j ? 0.0 : avg_off for i in 1:n, j in 1:n ]
                    p_simpl = (r, A_simpl)
                    prob_simpl = ODEProblem(full_model!, u0, tspan, p_simpl)
                    sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
                    rt_simpl = simulate_press_perturbation(u0, p_simpl, tspan, t_perturb, δ; solver=Tsit5())
                    mean_rt_simpl = mean(skipmissing(rt_simpl))
                    
                    # Step 7: Abundance distribution metrics.
                    fit_ln = fit(LogNormal, fixed_B)
                    mean_ln = exp(fit_ln.μ + (fit_ln.σ^2)/2)
                    var_ln = exp(2*fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2) - 1)
                    RelVar = var_ln / mean_ln^2
                    
                    # Step 8: Emergent network parameters.
                    μ_eff = n * mean(offdiag_vals)
                    σ_eff = sqrt(n) * std(offdiag_vals)
                    pairs = [(A[i, j], A[j, i]) for i in 1:n for j in 1:n if i < j]
                    γ_eff = length(pairs) > 1 ? cor(getindex.(pairs, 1), getindex.(pairs, 2)) : 0.0
                    
                    # Step 9: Predicted return time.
                    T_pred = predicted_return_time(φ, γ_eff, RelVar)
                    
                    # Step 10: Additional stability metrics.
                    resilience = compute_resilience(B_eq_full, A)
                    reactivity = compute_reactivity(B_eq_full, A)
                    
                    # Assemble a NamedTuple result for this iteration.
                    record = (species_count=n,
                              delta=δ,
                              iteration=iter,
                              connectance=conn,
                              persistence=φ,
                              return_full=mean_rt_full,
                              return_simpl=mean_rt_simpl,
                              return_pred=T_pred,
                              relvar=RelVar,
                              resilience=resilience,
                              reactivity=reactivity)
                    
                    lock(results_lock) do
                        push!(sim_results, record)
                    end
                catch e
                    @warn "Iteration $iter skipped for n=$n, conn=$conn, δ=$δ" exception = e
                end
            end  # End threads for this parameter combination
        end
    end
end

# Convert simulation results to a DataFrame.
df = DataFrame(sim_results)

# Display the first few rows.
println(first(df, 10))

#############################
# Plotting: Return Time vs. Species Count for Each Perturbation and for Each Metric (Full, Simplified, Predicted)
#############################
# Grouped by Connectance (Side-by-Side Boxplots)
#############################
begin
    
    unique_deltas = sort(unique(df.delta))
    n_deltas = length(unique_deltas)

    # Define the three metrics for the three rows.
    metrics = [:return_full, :return_simpl, :return_pred]
    titles = ["Return Time (Full Model)", "Return Time (Simplified Model)", "Return Time (Analytical Prediction)"]

    fig = Figure(; size = (350*n_deltas, 585))

    # Loop over columns (δ values) and rows (metrics)
    for (row_idx, d) in enumerate(unique_deltas)
        # Filter for the current delta value.
        df_d = filter(row -> row.delta == d, df)
        df_d.return_pred .= df_d.return_pred .* df_d.species_count
        # Get the x values (species count) and grouping variable for connectance.
        species_vals = df_d.species_count
        conn_vals = df_d.connectance
        
        conn_cat_dodge = Dict{Float64, Int64}()
        for (i, conn) in enumerate(unique(vals))
            conn_cat_dodge[conn] = i
        end
        conn_cat_dodge = [conn_cat_dodge[c] for c in conn_vals]
        
        for (col_idx, metric) in enumerate(metrics)
            ax = Axis(fig[row_idx, col_idx],
                xlabel = "Species Count",
                ylabel = string(metric),
                title = "δ = $(d): " * titles[col_idx],
                xticks = 10:10:30)
            # Boxplot: we pass species_vals as categories (numeric here) and dodge by dodge_vals.
            MK.boxplot!(ax, species_vals, df_d[!, metric],
                dodge = conn_cat_dodge, show_notch = true,
                color = conn_vals, gap = 10.0, dodge_gap = 0.1)
        end
    end

    # Optionally, add a colorbar for the grouping (connectance)
    # fig[end+1, 1] = Colorbar(fig, colormap = cgrad(:viridis), label = "Connectance")

    display(fig)

end

#############################
# SAME THING BUT FOR RELATIVE VARIANCE, RESILIENCE, AND REACTIVITY
#############################
begin
    
    unique_deltas = sort(unique(df.delta))
    n_deltas = length(unique_deltas)

    # Define the three metrics for the three rows.
    metrics = [:relvar, :resilience, :reactivity]
    titles = ["Relative Variance", "Resilience", "Reactivity"]

    fig = Figure(; size = (350*n_deltas, 585))

    # Loop over columns (δ values) and rows (metrics)
    for (row_idx, d) in enumerate(unique_deltas)
        # Filter for the current delta value.
        df_d = filter(row -> row.delta == d, df)
        # Get the x values (species count) and grouping variable for connectance.
        species_vals = df_d.species_count
        conn_vals = df_d.connectance
        
        conn_cat_dodge = Dict{Float64, Int64}()
        for (i, conn) in enumerate(unique(conn_cat))
            conn_cat_dodge[conn] = i
        end
        conn_cat_dodge = [conn_cat_dodge[c] for c in conn_cat]
        
        for (col_idx, metric) in enumerate(metrics)
            ax = Axis(fig[row_idx, col_idx],
                xlabel = "Species Count",
                ylabel = string(metric),
                title = "δ = $(d): " * titles[col_idx],
                xticks = 10:10:30)
            # Boxplot: we pass species_vals as categories (numeric here) and dodge by dodge_vals.
            MK.boxplot!(ax, species_vals, df_d[!, metric],
                dodge = conn_cat_dodge, show_notch = true,
                color = conn_vals, gap = 10.0, dodge_gap = 0.1)
        end
    end

    # Optionally, add a colorbar for the grouping (connectance)
    # fig[end+1, 1] = Colorbar(fig, colormap = cgrad(:viridis), label = "Connectance")

    display(fig)

end
