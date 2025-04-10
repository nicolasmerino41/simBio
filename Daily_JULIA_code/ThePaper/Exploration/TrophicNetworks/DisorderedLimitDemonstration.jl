#############################
# Function Definitions
#############################
# Fully random GLV predator-prey model.
# We assume D[i] = 1 so that r[i] = k[i].
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

    # Perturbed growth rates.
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
    return return_times
end

# Predicted return time (for illustration) using a simple formula.
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

# Resilience: negative of maximum real part of J's eigenvalues.
function compute_resilience(B, A)
    J = compute_jacobian(B, A)
    ev = eigvals(J)
    return -maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of (J+J')/2.
function compute_reactivity(B, A)
    J = compute_jacobian(B, A)
    Jsym = (J + J') / 2
    ev_sym = eigvals(Jsym)
    return maximum(real.(ev_sym))
end

# Compute variability from the Lyapunov equation: J V + V J' = diag(B)
function compute_variability(B, A)
    J = compute_jacobian(B, A)
    N = Diagonal(B)
    try
        # Using ControlSystems.lyap to solve the Lyapunov equation.
        V_mat = lyap(J, N)
        variability = trace(V_mat) / length(B)
    catch
        variability = NaN
    end
    return variability
end

#############################
# Pipeline Parameters
#############################
species_scenarios = [50]   # You may fix n or sweep over n; here we choose one example.
Niter = 200
tspan = (0.0, 100.0)
t_perturb = 50.0
connectance_list = 0.05:0.05:0.25
delta_list = [0.2, 0.4, 0.6]

# We'll collect simulation results in a vector of NamedTuples.
sim_results = Vector{NamedTuple{(:species_count, :connectance, :delta, :iteration, :persistence, :T, :Simpson_div, :resilience, :reactivity, :μ_eff, :σ_eff, :γ_eff, :ζ2), Tuple{Int, Float64, Float64, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()

results_lock = ReentrantLock()

#############################
# Main Simulation Loop (Threaded)
#############################
for n in species_scenarios
    for conn in connectance_list
        for δ in delta_list
            Threads.@threads for iter = 1:Niter
                try
                    # Step 1: Draw fixed equilibrium abundances B* from a lognormal distribution.
                    μ_log_obs = 0.5; σ_log_obs = 1.0
                    fixed_B = rand(LogNormal(μ_log_obs, σ_log_obs), n)

                    # Step 2: Construct a fully random trophic network with connectance.
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
                        local s = 0.0  # declare a local accumulator
                        for j in 1:n
                            if j != i
                                s += A[i, j] * fixed_B[j]
                            end
                        end
                        k[i] = fixed_B[i] - s
                    end

                    r = k   # As D[i] = 1 by assumption.

                    # Step 4: Simulate the full model (should be at equilibrium).
                    u0 = fixed_B; p = (r, A)
                    prob = ODEProblem(full_model!, u0, tspan, p)
                    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
                    B_eq = sol[:, end]
                    thresh = 1e-3
                    φ = sum(B_eq .> thresh) / n   # persistence

                    # Compute Total Biomass T.
                    T_val = sum(B_eq)

                    # Compute Simpson diversity.
                    rel_abund = B_eq / T_val
                    Simpson_inv = sum(rel_abund.^2)
                    Simpson_div = 1 / Simpson_inv

                    # Compute variability using Lyapunov.
                    # V_val = compute_variability(B_eq, A)

                    # Compute resilience and reactivity.
                    resilience_val = compute_resilience(B_eq, A)
                    reactivity_val = compute_reactivity(B_eq, A)

                    # Effective network parameters.
                    offdiag_vals = [A[i, j] for i in 1:n for j in 1:n if i != j]
                    μ_eff = n * mean(offdiag_vals)
                    σ_eff = sqrt(n) * std(offdiag_vals)
                    pairs = [(A[i, j], A[j, i]) for i in 1:n for j in 1:n if i < j]
                    γ_eff = length(pairs) > 1 ? cor(getindex.(pairs, 1), getindex.(pairs, 2)) : 0.0
                    ζ2 = var(k)

                    # Assemble the result.
                    record = (species_count = n,
                              connectance = conn,
                              delta = δ,
                              iteration = iter,
                              persistence = φ,
                              T = T_val,
                              Simpson_div = Simpson_div,
                            #   variability = V_val,
                              resilience = resilience_val,
                              reactivity = reactivity_val,
                              μ_eff = μ_eff,
                              σ_eff = σ_eff,
                              γ_eff = γ_eff,
                              ζ2 = ζ2)

                    lock(results_lock) do
                        push!(sim_results, record)
                    end

                catch e
                    @warn "Iteration $iter skipped for n=$n, conn=$conn, δ=$δ" exception = e
                end
            end  # End @threads loop
        end
    end
end

# Convert results to a DataFrame.
df = DataFrame(sim_results)
println(first(df, 10))

#############################
# Plotting: Compare Simulation Outcomes vs. Analytical Predictions
#############################
# For demonstration, let’s create a multi-panel (3x?) layout that shows the three return time metrics.
# In this case, we demonstrate the fully simulated properties.
# (Since the analytical solution of the full reference model involves solving self-consistent equations,
# here we simply illustrate the simulation outcomes and how they collapse as a function of the effective parameters.)

# For instance, we can create a layout showing:
# Row 1: Total Biomass T vs. species_count (grouped by connectance and delta)
# Row 2: Simpson diversity vs. species_count
# Row 3: Variability, resilience, and reactivity vs. species_count

begin
    # Assume df is your DataFrame with columns: species_count, delta, T, Simpson_div, resilience, reactivity, connectance.
    unique_conn = sort(unique(df.connectance))
    unique_delta = sort(unique(df.delta))

    n_deltas = length(unique_delta)
    fig = Figure(; size = (900, 600))

    for (col_idx, δ) in enumerate(unique_delta)
        
        # Build dodge mapping for connectance based on unique_conn (from the whole df, as in your example)
        conn_cat_dodge = Dict{Float64, Int64}()
        for (i, conn) in enumerate(unique_conn)
            conn_cat_dodge[conn] = i
        end
        conn_cat_dodge = [conn_cat_dodge[c] for c in df.connectance]
        
        # Filter for the current delta value.
        df_d = filter(row -> row.delta == δ, df)
        
        # Row 1: Total Biomass T.
        ax1 = Axis(fig[1, col_idx],
            xlabel = "Species Count",
            ylabel = "Total Biomass T",
            title = "δ = $(δ) : Total Biomass")
        boxplot!(ax1, df_d.species_count, df_d.T)

        # Row 2: Simpson diversity.
        ax2 = Axis(fig[2, col_idx],
            xlabel = "Species Count",
            ylabel = "Simpson Diversity",
            title = "δ = $(δ) : Simpson Diversity",
            xticks = 10:10:30)
        boxplot!(ax2, df_d.species_count, df_d.Simpson_div; 
            show_notch = true)

        # Row 3: Stability Metrics: combine resilience and reactivity.
        ax3 = Axis(fig[3, col_idx],
            xlabel = "Species Count",
            ylabel = "Stability Metrics",
            title = "δ = $(δ) : Resilience / Reactivity",
            xticks = 10:10:30)

        # Create a dodge mapping for connectance for df_d.
        dodge_mapping = Dict{Float64, Float64}()
        unique_conn_vals = unique(df_d.connectance)
        for (i, conn_val) in enumerate(unique_conn_vals)
            dodge_mapping[conn_val] = 0.1 * (i - 1)  # small offset
        end
        # Build dodge offsets for each metric.
        dodge_resilience = [dodge_mapping[c] for c in df_d.connectance]
        dodge_reactivity = [dodge_mapping[c] + 0.4 for c in df_d.connectance]

        boxplot!(ax3, df_d.species_count, df_d.resilience; label = "Resilience")
        boxplot!(ax3, df_d.species_count, df_d.reactivity; label = "Reactivity")
        # axislegend(ax3; position = :rb)
    end

    # Optionally add a colorbar to indicate Connectance.
    fig[1, end+1] = Colorbar(fig, colormap = cgrad(:viridis), label = "Connectance")
    display(fig)

end


fig[4, 1] = Colorbar(fig, colormap = cgrad(:viridis), label = "Connectance")
display(fig)
