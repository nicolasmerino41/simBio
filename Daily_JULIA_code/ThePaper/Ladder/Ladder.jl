# Full consumer–dynamics model (full model).
# Here, p is a tuple: (m, xi, epsilon, A, preys, predators)
function consumer_dynamics!(dB, B, p, t)
    m, xi, epsilon, A, preys, predators = p
    n = length(B)
    for i in 1:n
        # Sum over preys (if empty, returns 0.0)
        prey_effect = sum((epsilon * A[i,j] * B[j] for j in preys[i]), init=0.0)
        # Sum over predators (if empty, returns 0.0)
        pred_effect = sum((A[j,i] * B[j] for j in predators[i]), init=0.0)
        # Compute rate of change for species i
        dB[i] = B[i] * (m[i] / xi[i]) * (-xi[i] - B[i] + prey_effect - pred_effect)
    end
end

# Calibration function: adjust thresholds (xi) to force equilibrium at B_star.
function calibrate_thresholds(B_star::Vector{Float64}, p)
    m, xi, epsilon, A, preys, predators = p
    n = length(B_star)
    xi_cal = zeros(n)
    for i in 1:n
        prey_term = sum((epsilon * A[i,j] * B_star[j] for j in preys[i]), init=0.0)
        pred_term = sum((A[j,i] * B_star[j] for j in predators[i]), init=0.0)
        xi_cal[i] = B_star[i] + prey_term - pred_term
    end
    return xi_cal
end

# Simulate a press perturbation.
# Here the perturbation is applied by increasing the threshold xi by a fraction delta.
function simulate_press_perturbation(u0, p, tspan, t_perturb, delta; solver=Tsit5(), plot=false)
    m, xi, epsilon, A, preys, predators = p
    # Phase 1: simulate until the time of press perturbation.
    tspan1 = (tspan[1], t_perturb)
    prob1 = ODEProblem(consumer_dynamics!, u0, tspan1, p)
    sol1 = solve(prob1, solver, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]
    
    # Phase 2: apply the press perturbation by modifying xi.
    xi_press = copy(xi)
    xi_press .= xi_press .* (1 + delta)
    p_press = (m, xi_press, epsilon, A, preys, predators)
    
    tspan2 = (t_perturb, tspan[2])
    prob2 = ODEProblem(consumer_dynamics!, pre_state, tspan2, p_press)
    sol2 = solve(prob2, solver, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)
    
    # Compute return times (time required for each species to get within 10% of new equilibrium)
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
        fig = Figure(; size = (800,600))
        ax = Axis(fig[1,1], xlabel="Time", ylabel="Abundance", title="Community Response")
        for i in 1:n
            # Plot time series for species i
            lines!(ax, sol2.t, [state[i] for state in sol2.u], label="Species $i")
        end
        axislegend(ax, position = :rb)
        display(fig)
    end
    return return_times
end

# Compute the Jacobian at equilibrium.
function compute_jacobian(B, p)
    m, xi, epsilon, A, preys, predators = p
    n = length(B)
    J = zeros(n, n)
    psi = B ./ xi  # relative yield term
    for i in 1:n
        for j in 1:n
            δ = (i == j) ? 1.0 : 0.0
            J[i,j] = m[i] * psi[i] * (-δ + epsilon * A[i,j] - A[j,i])
        end
    end
    return J
end

# Resilience: negative of the largest real part of the Jacobian eigenvalues.
function compute_resilience(B, p)
    J = compute_jacobian(B, p)
    ev = eigvals(J)
    return -maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of the symmetric part of the Jacobian.
function compute_reactivity(B, p)
    J = compute_jacobian(B, p)
    Jsym = (J + J') / 2
    ev_sym = eigvals(Jsym)
    return maximum(real.(ev_sym))
end

# Predicted return time (example prediction function; adjust as needed).
function predicted_return_time(phi, gamma_eff, relVar)
    denom = 1 - phi * gamma_eff * relVar
    return denom > 0 ? 1 / denom : Inf
end

#############################
# Pipeline Parameters
#############################

species_scenarios = [10, 30, 50]         # Number of species
Niter = 10                               # Iterations per scenario (adjust as needed)
tspan = (0.0, 100.0)                     # Simulation time
t_perturb = 50.0                         # Press perturbation time
connectance_list = 0.1:0.1:0.3             # Range for network connectance
delta_list = [0.1, 0.2, 0.3]               # Perturbation intensities

# We'll save all simulation results here.
sim_results = Vector{NamedTuple{(
    :species_count, :delta, :iteration, :connectance, 
    :persistence, :return_full, :return_simpl, :return_pred, 
    :relvar, :resilience, :reactivity), Tuple{Int, Float64, Int, Float64, 
    Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()

# Lock to protect shared array writes when using threads.
results_lock = ReentrantLock()

#############################
# Main Simulation Loop (Threads)
#############################
for n in species_scenarios
    for conn in connectance_list
        for delta in delta_list
            Threads.@threads for iter = 1:Niter
                try
                    # Step 1: Generate fixed equilibrium abundances.
                    fixed_B = rand(LogNormal(0.5, 1.0), n)
                    
                    # Step 2: Construct the trophic network.
                    A = zeros(n, n)
                    preys = [Int[] for i in 1:n]
                    predators = [Int[] for i in 1:n]
                    for i in 1:n
                        for j in 1:n
                            if i != j && rand() < conn
                                strength = rand() * 0.3  # interaction strength (adjust scale for stability)
                                if rand() < 0.5
                                    # species i preys on j
                                    A[i,j] = strength
                                    push!(preys[i], j)
                                    push!(predators[j], i)
                                else
                                    # species j preys on i
                                    A[j,i] = strength
                                    push!(preys[j], i)
                                    push!(predators[i], j)
                                end
                            end
                        end
                    end
                    
                    # Step 3: Calibrate the model.
                    # Species-specific mortality (random between 0.1 and 1.0)
                    m = rand(Uniform(0.1, 1.0), n)
                    xi_guess = ones(n)
                    epsilon = 0.1  # Lower conversion efficiency to help stability
                    
                    # Build initial parameters tuple.
                    p_model = (m, xi_guess, epsilon, A, preys, predators)
                    
                    # Calibrate xi so that fixed_B is an equilibrium.
                    xi_cal = calibrate_thresholds(fixed_B, p_model)
                    p_model = (m, xi_cal, epsilon, A, preys, predators)
                    
                    # Step 4: Full model simulation.
                    u0 = fixed_B
                    prob = ODEProblem(consumer_dynamics!, u0, tspan, p_model)
                    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
                    B_eq_full = sol.u[end]
                    thresh = 1e-3
                    phi = sum(B_eq_full .> thresh) / n   # Persistence metric
                    
                    # Step 5: Press perturbation simulation (full model).
                    rt_full = simulate_press_perturbation(u0, p_model, tspan, t_perturb, delta; solver=Tsit5(), plot=true)
                    mean_rt_full = mean(skipmissing(rt_full))
                    
                    # Step 6: Simplified model simulation.
                    # Create a simplified interaction matrix with average off-diagonal strength.
                    offdiag_vals = [ A[i,j] for i in 1:n for j in 1:n if i != j ]
                    avg_off = isempty(offdiag_vals) ? 0.0 : mean(offdiag_vals)
                    A_simpl = [i == j ? 0.0 : avg_off for i in 1:n, j in 1:n]
                    p_simpl = (m, xi_cal, epsilon, A_simpl, preys, predators)
                    
                    # Run simplified full simulation.
                    prob_simpl = ODEProblem(consumer_dynamics!, u0, tspan, p_simpl)
                    sol_simpl = solve(prob_simpl, Tsit5(), reltol=1e-8, abstol=1e-8)
                    rt_simpl = simulate_press_perturbation(u0, p_simpl, tspan, t_perturb, delta; solver=Tsit5(), plot=false)
                    mean_rt_simpl = mean(skipmissing(rt_simpl))
                    
                    # Step 7: Abundance distribution metrics.
                    fit_ln = fit(LogNormal, fixed_B)
                    mean_ln = exp(fit_ln.μ + (fit_ln.σ^2) / 2)
                    var_ln = exp(2 * fit_ln.μ + fit_ln.σ^2) * (exp(fit_ln.σ^2) - 1)
                    relVar = var_ln / mean_ln^2
                    
                    # Step 8: Emergent network parameters.
                    μ_eff = n * mean(offdiag_vals)
                    σ_eff = sqrt(n) * std(offdiag_vals)
                    pairs = [(A[i,j], A[j,i]) for i in 1:n for j in 1:n if i < j]
                    gamma_eff = length(pairs) > 1 ? cor(getindex.(pairs, 1), getindex.(pairs, 2)) : 0.0
                    
                    # Step 9: Predicted return time.
                    T_pred = predicted_return_time(phi, gamma_eff, relVar)
                    
                    # Step 10: Additional stability metrics.
                    resilience = compute_resilience(B_eq_full, p_model)
                    reactivity = compute_reactivity(B_eq_full, p_model)
                    
                    # Assemble results for this iteration.
                    record = (species_count = n,
                              delta = delta,
                              iteration = iter,
                              connectance = conn,
                              persistence = phi,
                              return_full = mean_rt_full,
                              return_simpl = mean_rt_simpl,
                              return_pred = T_pred,
                              relvar = relVar,
                              resilience = resilience,
                              reactivity = reactivity)
                    
                    lock(results_lock) do
                        push!(sim_results, record)
                    end
                catch e
                    @warn "Iteration $iter skipped for n=$n, conn=$conn, delta=$delta" exception = e
                end
            end
        end
    end
end

#############################
# Post-Processing: Create DataFrame
#############################

df = DataFrame(sim_results)
println(first(df, 10))
