# -----------------------------
# Function Definitions
# -----------------------------
# Full GLV predator–prey model.
function full_model!(du, B, p, t)
    r, A = p
    n = length(B)
    for i in 1:n
        du[i] = B[i] * (r[i] - B[i] + sum(A[i, :] .* B))
    end
end

# Simulate the press perturbation and compute stability metrics:
# Return time (T_R): the time at which a species first recovers within 10% of new equilibrium.
# Maximum overshoot (M_OS): maximum relative deviation from the new equilibrium during recovery.
# Integrated recovery error (I_RE): mean relative error over the recovery period.
function simulate_stability_metrics(u0, p, tspan, t_perturb, delta; solver=Tsit5())
    # Integrate up to the perturbation.
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(full_model!, u0, tspan1, p)
    sol1 = solve(prob1, solver, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]
    
    # Apply the press perturbation.
    r, A = p
    r_press = (1 - delta) .* r
    p_press = (r_press, A)
    
    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(full_model!, pre_state, tspan2, p_press)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)
    new_eq = sol2.u[end]
    n = length(new_eq)
    
    T_R_arr = fill(NaN, n)
    M_OS_arr = zeros(n)
    I_RE_arr = zeros(n)
    counts = zeros(n)
    
    for i in 1:n
        max_err = 0.0
        sum_err = 0.0
        found = false
        for (t, state) in zip(sol2.t, sol2.u)
            err = abs(state[i] - new_eq[i]) / (abs(new_eq[i]) + 1e-8)
            sum_err += err
            counts[i] += 1
            max_err = max(max_err, err)
            if !found && err < 0.1
                T_R_arr[i] = t - t_perturb
                found = true
            end
        end
        M_OS_arr[i] = max_err
        I_RE_arr[i] = counts[i] > 0 ? sum_err / counts[i] : NaN
    end
    
    # Return the average over species.
    return (T_R = mean(skipmissing(T_R_arr)),
            M_OS = mean(M_OS_arr),
            I_RE = mean(I_RE_arr))
end

# Compute effective network parameters from A and calibrated carrying capacities k.
function compute_network_properties(A, k, n)
    offdiag = Float64[]
    for i in 1:n
        for j in 1:n
            if i != j
                push!(offdiag, A[i, j])
            end
        end
    end
    μ_eff = n * mean(offdiag)
    σ_eff = sqrt(n) * std(offdiag)
    pairs = [(A[i, j], A[j, i]) for i in 1:n for j in 1:n if i < j]
    γ_eff = length(pairs) > 1 ? cor(getindex.(pairs, 1), getindex.(pairs, 2)) : 0.0
    ζ2 = var(k)
    return (μ_eff = μ_eff, σ_eff = σ_eff, γ_eff = γ_eff, ζ2 = ζ2)
end

# Compute degree-based metrics (for additional structural information).
function compute_degree_metrics(A, n)
    degrees = zeros(n)
    for i in 1:n
        cnt = 0
        for j in 1:n
            if i != j && abs(A[i, j]) > 1e-8
                cnt += 1
            end
        end
        degrees[i] = cnt
    end
    mean_deg = mean(degrees)
    cv_deg = (mean_deg > 0) ? std(degrees) / mean_deg : 0.0
    return (mean_deg = mean_deg, cv_deg = cv_deg)
end

# Analytical prediction for return time (basic model), using fixed_B, persistence φ, and effective reciprocity γ_eff.
function predicted_return_time_basic(fixed_B, φ, γ_eff)
    meanB = mean(fixed_B)
    varB = var(fixed_B)
    RelVar = varB / (meanB^2 + 1e-8)
    denom = 1 - φ * γ_eff * RelVar
    return denom > 0 ? 1 / denom : Inf
end

# Extended prediction incorporates degree CV.
function predicted_return_time_extended(fixed_B, φ, γ_eff, cv_deg; α_degree=0.5)
    basic = predicted_return_time_basic(fixed_B, φ, γ_eff)
    return basic * (1 + α_degree * cv_deg)
end

# -----------------------------
# Pipeline Parameters & Data Collection
# -----------------------------
species_scenarios = [10, 30, 50]   # You can add more if desired.
Niter = 100
tspan = (0.0, 100.0)
t_perturb = 50.0
connectance_list = 0.05:0.05:0.25
delta_list = [0.2, 0.4, 0.6]

# We'll store each simulation run as a NamedTuple.
results_vec = Vector{NamedTuple{(
    :species_count, :connectance, :delta, :iteration, :persistence,
    :T_sim, :M_OS, :I_RE, :μ_eff, :σ_eff, :γ_eff, :ζ2, :mean_deg, :cv_deg,
    :T_pred_basic, :T_pred_ext, :err_RT_basic, :err_RT_ext, :err_OS, :err_I_RE, :compound_err_basic, :compound_err_ext),
    Tuple{Int, Float64, Float64, Int, Float64,
          Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
          Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()

results_lock = ReentrantLock()

# -----------------------------
# Main Simulation Loop (Threaded)
# -----------------------------
for n in species_scenarios
    for conn in connectance_list
        for δ in delta_list
            Threads.@threads for iter = 1:Niter
                try
                    # Step 1: Generate fixed equilibrium abundances.
                    μ_log_obs = 0.5; σ_log_obs = 1.0
                    fixed_B = rand(LogNormal(μ_log_obs, σ_log_obs), n)
                    
                    # Step 2: Construct a fully random trophic network with the given connectance.
                    A = zeros(n, n)
                    for i in 1:n
                        for j in 1:n
                            if i != j && rand() < conn
                                if rand() < 0.5
                                    β = rand(Exponential(1.0))
                                    A[i, j] = 0.3 * β
                                    A[j, i] = -β
                                else
                                    β = rand(Exponential(1.0))
                                    A[j, i] = 0.3 * β
                                    A[i, j] = -β
                                end
                            end
                        end
                    end
                    
                    # Step 3: Calibrate the model so that fixed_B is an equilibrium.
                    k = zeros(n)
                    for i in 1:n
                        local s = 0.0
                        for j in 1:n
                            if j != i
                                s += A[i, j] * fixed_B[j]
                            end
                        end
                        k[i] = fixed_B[i] - s
                    end
                    r = k  # since D[i] = 1
                                        
                    # Step 4: Simulate the full model.
                    u0 = fixed_B
                    p = (r, A)
                    prob = ODEProblem(full_model!, u0, tspan, p)
                    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
                    B_eq = sol[:, end]
                    thresh = 1e-3
                    φ = sum(B_eq .> thresh) / n
                    
                    # Step 5: Compute the three stability metrics from perturbation.
                    stab = simulate_stability_metrics(u0, p, tspan, t_perturb, δ; solver=Tsit5())
                    T_sim = stab.T_R
                    M_OS_val = stab.M_OS
                    I_RE_val = stab.I_RE
                    
                    # Step 6: Compute effective network parameters.
                    netp = compute_network_properties(A, k, n)
                    μ_eff = netp.μ_eff
                    σ_eff = netp.σ_eff
                    γ_eff = netp.γ_eff
                    ζ2 = netp.ζ2
                    
                    # Step 7: Compute degree-based metrics.
                    degp = compute_degree_metrics(A, n)
                    mean_deg = degp.mean_deg
                    cv_deg = degp.cv_deg
                    
                    # Step 8: Compute predictions.
                    T_pred_basic = predicted_return_time_basic(fixed_B, φ, γ_eff)
                    T_pred_ext = predicted_return_time_extended(fixed_B, φ, γ_eff, cv_deg; α_degree=0.5)
                    
                    # Step 9: Compute relative errors for Return Time.
                    err_RT_basic = abs(T_sim - T_pred_basic) / ((T_sim + T_pred_basic) / 2 + 1e-8)
                    err_RT_ext   = abs(T_sim - T_pred_ext)   / ((T_sim + T_pred_ext)   / 2 + 1e-8)
                    # For overshoot and integrated error, we assume the ideal prediction is 0.
                    err_OS = M_OS_val  # absolute value, as target = 0.
                    err_I_RE = I_RE_val
                    
                    # Compound errors: average of the three relative errors.
                    compound_err_basic = (err_RT_basic + err_OS + err_I_RE) / 3
                    compound_err_ext   = (err_RT_ext   + err_OS + err_I_RE) / 3
                    
                    # Assemble the record.
                    rec = (
                        species_count = n,
                        connectance = conn,
                        delta = δ,
                        iteration = iter,
                        persistence = φ,
                        T_sim = T_sim,
                        M_OS = M_OS_val,
                        I_RE = I_RE_val,
                        μ_eff = μ_eff,
                        σ_eff = σ_eff,
                        γ_eff = γ_eff,
                        ζ2 = ζ2,
                        mean_deg = mean_deg,
                        cv_deg = cv_deg,
                        T_pred_basic = T_pred_basic,
                        T_pred_ext = T_pred_ext,
                        err_RT_basic = err_RT_basic,
                        err_RT_ext = err_RT_ext,
                        err_OS = err_OS,
                        err_I_RE = err_I_RE,
                        compound_err_basic = compound_err_basic,
                        compound_err_ext = compound_err_ext
                    )
                    lock(results_lock) do
                        push!(results_vec, rec)
                    end
                    
                catch e
                    @warn "Iteration $iter skipped for n=$n, conn=$conn, δ=$δ" exception = e
                end
            end  # End @threads loop
        end
    end
end

# Convert results to a DataFrame.
df = DataFrame(results_vec)
df = dropmissing(df)
println("First few rows:")
println(first(df, 10))

# Print overall compound errors averaged over all simulations.
compound_err_basic_overall = mean(df.compound_err_basic)
compound_err_ext_overall   = mean(df.compound_err_ext)
println("Overall Compound Relative Error (Basic): ", compound_err_basic_overall)
println("Overall Compound Relative Error (Extended): ", compound_err_ext_overall)

# -----------------------------
# Plotting Examples
# -----------------------------
# Plotting Example: Compare Simulated Return Time vs. Predicted (Basic & Extended)
begin
    # Example 1: Scatter plot comparing simulated T_sim vs predicted T (basic & extended)
    fig1 = Figure(; size = (1200, 600))
    ax1 = Axis(fig1[1,1], xlabel = "Simulated Return Time", ylabel = "Predicted Return Time (Basic)", 
        title = "Basic Prediction vs Simulation")
    MK.scatter!(ax1, df.T_sim, df.T_pred_basic, markersize=8, color = df.connectance, colormap = cgrad(:viridis))
    Colorbar(fig1[1,2], limits=(minimum(df.connectance), maximum(df.connectance)), colormap=cgrad(:viridis), label="Connectance")

    ax2 = Axis(fig1[2,1], xlabel = "Simulated Return Time", ylabel = "Predicted Return Time (Extended)", 
        title = "Extended Prediction vs Simulation")
    MK.scatter!(ax2, df.T_sim, df.T_pred_ext, markersize=8, color = df.cv_deg, colormap = cgrad(:plasma))
    Colorbar(fig1[2,2], limits=(minimum(df.cv_deg), maximum(df.cv_deg)), colormap=cgrad(:plasma), label="CV of Degree")

    display(fig1)
end

begin
    # Example 2: Boxplots of Compound Errors by Connectance (grouped by δ)
    unique_delta = sort(unique(df.delta))
    fig2 = Figure(; size = (400*length(unique_delta), 600))
    for (col_idx, δ) in enumerate(unique_delta[3])
        df_d = filter(row -> row.delta == δ, df)
        ax = Axis(fig2[1, col_idx], xlabel = "Connectance", ylabel = "Compound Err (Basic)", 
            title = "δ = $(δ)")
        MK.boxplot!(ax, df_d.connectance, df_d.compound_err_basic; show_notch=true)
    end
    display(fig2)
end