################ INFORMATION #######################
# a) Ladder1.jl is for all the functions
# b) Ladder2.jl is for running the simulations
# c) Ladder3.jl is for post-processing and plotting
include("Ladder1.jl")
using Logging
logger = SimpleLogger(stdout, Logging.Info)
# Master pipeline function.
"""
    run_simulation(R_vals, C_vals, conn_vals, delta_vals; Niter=10, tspan=(0.0,50.0), t_perturb=25.0, plot=false)

For each combination of resource count, consumer count, connectance, and press perturbation intensity,
this function:
  - Generates a target equilibrium state from a LogNormal.
  - Constructs a global interaction matrix A (with resources not preying).
  - Defines model parameters and calibrates them (using calibrate_params, assumed defined).
  - Simulates the full model to obtain equilibrium, computes persistence, resilience, and reactivity.
  - Simulates a press perturbation for both the full and a simplified model.
  - Returns a DataFrame of results.
"""
function run_simulation(
    R_vals::Vector{Int}, C_vals::Vector{Int},
    conn_vals::AbstractRange{Float64}, delta_vals::AbstractRange{Float64};
    Niter::Int=10,
    tspan::Tuple{Float64,Float64}=(0.0,50.0),
    t_perturb::Float64=25.0,
    plot::Bool=false, 
    show_warnings::Bool=true,
    max_calib::Int=10
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
                            m_cons = fill(0.1, C)                   # rand(Uniform(0.1, 1.0), C)    # Mortality rates for consumers.
                            
                            d_res = fill(2.0, R)                       # Scaling factors for resources.
                                
                            epsilon = max.(0, min.(1, rand(Normal(0.3, 0.1), total_species, total_species)))

                            p_tuple = (R, C, m_cons, d_res, epsilon, A)
                            
                            # Step 4: Calibrate thresholds so that fixed_state is an equilibrium.
                            R_eq = fixed_state[1:R]
                            C_eq = fixed_state[R+1:total_species]
                            new_xi_cons, new_r_res, A = calibrate_params(R_eq, C_eq, p_tuple)
                            count_calib = 0
                            
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
                                                        
                            # Step 5: Simulate full model.
                            prob_full = ODEProblem(trophic_ode!, fixed_state, tspan, p_tuple)
                            if show_warnings
                                sol_full = solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)
                            else
                                sol_full = with_logger(logger) do
                                    solve(prob_full, Tsit5(), reltol=1e-8, abstol=1e-8)
                                end
                            end

                            B_eq_full = sol_full.u[end]
                            
                            if sol_full.t[end] < t_perturb || any(isnan, sol_full.u[end]) || any(isinf, sol_full.u[end])
                                error("Error: solution did not finish properly")
                            end

                            println("New consumer thresholds (xi_cons): ", new_xi_cons)
                            println("New resource growth rates (r_res): ", new_r_res)

                            persistence = sum(x -> x > EXTINCTION_THRESHOLD, B_eq_full) / total_species
                            resilience = compute_resilience(B_eq_full, p_tuple)
                            reactivity = compute_reactivity(B_eq_full, p_tuple)
                            
                            # Step 6: Full model press perturbation.
                            rt_full, before_persistence_full, new_equil_full = simulate_press_perturbation(
                                fixed_state, p_tuple, tspan, t_perturb, delta;
                                solver=Tsit5(), plot=plot, show_warnings=show_warnings,
                                full_or_simple=true
                            )

                            mean_rt_full = mean(filter(!isnan, rt_full))
                            after_persitence_full = sum(x -> x > EXTINCTION_THRESHOLD, new_equil_full) / total_species
                            
                            # Step 7: Simplified model version.
                            offdiag_values = [A[i,j] for i in 1:total_species, j in 1:total_species if i != j && A[i,j] > 0.0]
                            avg_off = isempty(offdiag_values) ? 0.0 : mean(offdiag_values)
                            A_simpl = copy(A)
                            # scalar fill for positive entries
                            A_simpl[A .> 0.0] .= avg_off

                            # element‐wise weight for negative entries
                            mask = A .< 0.0
                            A_simpl[mask] .= -avg_off .* epsilon[mask]

                            p_simpl = (R, C, m_cons, new_xi_cons, new_r_res, d_res, epsilon, A_simpl)

                            rt_simpl, before_persistence_simpl, new_equil_simpl = simulate_press_perturbation(
                                fixed_state, p_simpl, tspan, t_perturb, delta;
                                solver=Tsit5(), plot=plot, show_warnings=show_warnings,
                                full_or_simple=false
                            )
                            
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
                                before_persistence_full = before_persistence_full,
                                after_persistence_full = after_persitence_full,
                                before_persistence_simpl = before_persistence_simpl,
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

#############################################################################
#############################################################################
################## Run Simulations and Gather Results #######################
#############################################################################
#############################################################################
# Specify parameter ranges:
R_vals = [20]            # e.g., try 4 and 6 resources
C_vals = [5]            # e.g., try 5 and 8 consumers
connectance_list = 0.1:0.1:0.5
delta_list = 0.1:0.1:0.1

# Run the simulation pipeline.
df_results2 = run_simulation(
    R_vals, C_vals, connectance_list, delta_list;
    Niter=10,
    tspan=(0.0,500.0), t_perturb=250.0, 
    plot=true, show_warnings=true,
    max_calib=10
)
println("First 10 rows of simulation results:")
println(first(df_results, 10))
