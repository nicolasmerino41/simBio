using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")

# ----------------------------------------------------------------
# 5) Stability & survival checks
# -------------------------------------------------------------------------------
function is_locally_stable(J::AbstractMatrix)
    if any(!isfinite, J)
        return false
    end
    lambda = eigvals(J)
    maximum(real.(lambda)) < 0 
end

function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false,sol.u[end]) : (all(sol.u[end] .> EXTINCTION_THRESHOLD),sol.u[end])
end

# --------------------------------------------------------------------------------
# Main sweep: fix S & C, vary structure, collect full & ladder pers.
# --------------------------------------------------------------------------------
function short_ComputingLadder(
    S::Int=50, C::Int=20;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[0.1,1.0],
    IS_vals_B_term=[0.1, 1.0],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[0.1],
    delta_vals=[1.0,3.0],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100,
    B_term = false,
    iterations=1,
    Rmed_iterations=200,
)
    R = S - C
    # A = zeros(S,S)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S,EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        conn_vals,IS_vals,scenarios,delta_vals,eps_scales,mortality_vals,growth_vals,1:iterations
    ))
    println("Number of combinations: ", length(combos))
    
    @threads for (conn,IS,scen,delta,epsi,m_val,g_val,ite) in sample(combos, min(length(combos), number_of_combinations); replace=false)
        
        # 1) build A & epsilon
        local A = make_A(zeros(S,S),R,conn,scen; IS = IS, B_term = B_term)
        # for i in R+1:S, j in R+1:S
        #     if A[i, j] < 0.0
        #         A[i, j] = A[i, j] * 0.01
        #     end
        # end

        local epsilon = clamp.(randn(S,S).*epsi,0,1)

        # --- NEW: compute collectivity phi ---
        local psi = compute_collectivity(copy(A), copy(epsilon))
        # backup_A = deepcopy(A)
        # backup_epsilon = deepcopy(epsilon)
        
        ########### #TODO ACTIVATE THIS PART IF YOU WANT TO GO BACK TO ENSURING FULL FEASIBLE #############
        # 2) calibrate Xi,K ? get (R_eq,C_eq,xi_cons,K_res)
        # 2) generate *all* feasible (Xi,K) sets for this A,epsilon
        old_epsilon = copy(epsilon)
        thr_sets = generate_feasible_thresholds(A, epsilon, R) #; margins = 1.0)
        if isempty(thr_sets)
            # @warn "No feasible thresholds for $conn, $IS, $scen, $epsi"
            continue
        end
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside generate_feasible_thresholds!"

        isempty(thr_sets) && ( @warn "No feasible thresholds for $conn, $IS, $scen, $epsi)" ; continue )

        # 2b) now loop over each feasible threshold-set
        local t = thr_sets[1]
        if any(t.C_eq .<= 0) || any(t.R_eq .<= 0)
            @error "Non-positive equilibrium abundances squized in the loop: $t"
            continue
        end
        xi_cons   = t.xi_cons
        K_res     = t.K_res
        R_eq        = t.R_eq
        C_eq        = t.C_eq
        #########################################################################################
        #########################################################################################
        # set up parameters & run the rest of your pipeline exactly as before
        # m_cons = fill(m_val, C)
        m_cons = abs.(rand(Normal(m_val, 0.2), C))
        # d_res  = ones(R)
        # r_res  = fill(g_val, R)
        r_res  = abs.(rand(Normal(g_val, 0.2), R))

        p     = (R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)
        fixed = vcat(R_eq, C_eq)

        # 3) stability & survival
        ok, B0 = survives!(fixed, p; cb=cb)
        # !ok && continue
        if !all(isapprox.(B0, fixed, atol=1e-3))
            @warn "B0 is not close to fixed point: B0 =  $(B0), and fixed = $(fixed)"
            continue
        end
        D, M = compute_jacobian(fixed, p)
        J_full = D * M
        !is_locally_stable(J_full) && continue
        
        # resilience_full is -max Re(eig(J)), but compute_resilience returns max Re
        old_epsilon = copy(epsilon)
        resilience_full = compute_resilience(fixed, p)
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside compute_resilience!"
        reactivity_full =  compute_reactivity(fixed, p)
        
        original_k_xi = vcat(K_res, xi_cons)
        # but you can also sweep t=0.1, 1, 10, etc.
        t0 = 0.01

        # full-model median return rate
        # Rmed_full = median_return_rate(J_full, fixed; t=t0, n=Rmed_iterations)
        # ssp_rmed_full = species_return_rates(J_full, fixed; t=t0, n=Rmed_iterations)
        analytical_rmed_full = analytical_species_return_rates(J_full; t=0.001)
        ############### TRYING SOMETHING NEW ################
        prob1 = ODEProblem(trophic_ode!, B0, (tspan[1], tpert), p)
        sol1 = solve(prob1, Tsit5(); callback = cb, reltol=1e-8, abstol=1e-8)
        pre_state = sol1.u[end]
        local manual_before_persistence = count(x -> x > EXTINCTION_THRESHOLD, pre_state) / length(pre_state)
        ###############################################################
        # 4) full-model persistence
        local rt_press, os, ire, before_full, after_persistence_full, new_equil, _ = simulate_press_perturbation(
            B0, p, tspan, tpert, delta;
            solver=Tsit5(),
            plot=false,
            show_warnings=true,
            full_or_simple=true,
            cb=cb,
            species_specific_perturbation=false
        )

        if manual_before_persistence != before_full || !isone(manual_before_persistence)
            error("manual_before_persistence = $(manual_before_persistence), before_full = $(before_full)")
            continue
        end

        # 4a) full-model pulse
        rt_pulse, _, after_pulse_full, _, _ = simulate_pulse_perturbation(
            B0, p, tspan, tpert, delta;
            solver=Tsit5(),
            plot=false,
            cb=cb,
            species_specific_perturbation=false
        )
        rt_press_full = mean(filter(!isnan, rt_press))
        rt_pulse_full = mean(filter(!isnan, rt_pulse))
        rt_pulse_full_vector = rt_pulse
        rt_press_full_vector = rt_press
        collectivity_full = psi 
        tau_full = -1 ./ diag(J_full)
        mean_tau_full = mean(tau_full)
        
        J_diff_full = norm(J_full - J_full)
        J_full_norm = norm(J_full)
        
        # min_delta_K_full = fill(0.0, R)
        
        # for i in 1:R
        #     min_delta_K_full[i] = find_min_K_reduction(
        #         i, B0, p, tspan, tpert; cb
        #     )
        # end
        
        # min_delta_xi_full = fill(0.0, C)
        # for i in 1:C
        #     min_delta_xi_full[i] = find_min_xi_increase(
        #         i, B0, p, tspan, tpert; cb
        #     )
        # end
        # mean_min_delta_K_full = mean(min_delta_K_full)
        # mean_min_delta_xi_full = mean(min_delta_xi_full)

        # 5) ladder persistence
        before_persistence_S = Dict(i => NaN for i in 1:11)
        after_persistence_S  = Dict(i => NaN for i in 1:11)
        after_pulse_S = Dict(i => NaN for i in 1:11)
        rt_press_S   = Dict(i => NaN for i in 1:11)
        rt_pulse_S   = Dict(i => NaN for i in 1:11)
        collectivity_S = Dict(i => NaN for i in 1:11)
        resilience_S  = Dict(i=>NaN for i in 1:11)
        reactivity_S  = Dict(i=>NaN for i in 1:11)
        stable_S    = Dict(i=>NaN for i in 1:11)
        
        # Rmed_S    = Dict(i=>NaN for i in 1:11)
        # ssp_rmed_S = Dict(i => Float64[] for i in 1:11)
        analytical_rmed_S = Dict(i => Float64[] for i in 1:11)
        
        tau_S = Dict(i => Float64[] for i in 1:11)
        mean_tau_S = Dict(i => NaN for i in 1:11)
        K_Xi_S = Dict(i => Float64[] for i in 1:11)
        J_diff_S = Dict(i => NaN for i in 1:11)
        min_delta_K_S = Dict(i => Float64[] for i in 1:11)
        min_delta_xi_S = Dict(i => Float64[] for i in 1:11)
        mean_min_delta_K_S = Dict(i => NaN for i in 1:11)
        mean_min_delta_xi_S = Dict(i => NaN for i in 1:11)
        rt_press_vector_S = Dict(i => Float64[] for i in 1:11)
        rt_pulse_vector_S = Dict(i => Float64[] for i in 1:11)
        @info "Running ladder"

        # original equilibrium abundances
        # R_eq_full, C_eq_full = B0[1:R], B0[R+1:S] # B0 is the simulated equilibrium
        R_eq_full, C_eq_full = B0[1:R], B0[R+1:S] # B0 is the calibrated equilibrium

        ###########################################################################################
        ###########################################################################################
        ###########################################################################################
        ########################## SIMPLIFIED MODEL STEPS #########################################
        ###########################################################################################
        ###########################################################################################
        for step in 1:3
            # build simplified A and epsilon
            local A_s, epsilon_s = short_transform_for_ladder_step(step, copy(A), copy(epsilon))
            psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))
            
            if step == 1
                # 1) check that A, epsilon are unchanged
                @assert all(A_s .== A)             "A was mutated on step $step!"
                @assert all(epsilon_s .== epsilon) "epsilon was mutated on step $step!"
                @assert isapprox(psi_s, psi; atol=1e-12)  "collectivity mismatch at step $step: psi=$psi, psi_s=$(psi_s)"
            end
            
            # 5a) Recompute xi_hat
            xi_hat = zeros(C)
            for k in 1:C
                i = R + k
                # 1) feeding gains (A>0)
                gain = sum(epsilon_s[i,j]*A_s[i,j]*R_eq_full[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
                # 2) predation losses: A_s[i,j]<0, but we need "+(-A)B"
                loss = sum(A_s[i,j]*C_eq_full[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                # consumer eq: xi = B_i - gain - loss
                xi_hat[k] = -C_eq_full[k] + gain + loss
            end

            # 5b) Recompute K_hat
            K_hat = zeros(R)
            for i in 1:R
                # resource eq uses A[i,j] (j consumer) directly:
                drain = sum(A_s[i,j]*C_eq_full[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0)
                # K_i = B_i + ? A[i,j] B_j
                K_hat[i] = abs(-R_eq_full[i] + drain)
            end

            # 5c) Solve for new equilibrium
            eq = try
                    calibrate_from_K_xi(xi_hat, K_hat, epsilon_s, A_s)
                catch err
                    @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                    continue
                end
    
            R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
            # also guard against non-finite or non-positive solution
            # if any(!isfinite, eq) || any(x->x<=0, eq)
            #     @warn "Step $step: infeasible equilibrium (non-finite or =0)"
            #     continue
            # end
            
            K_Xi_S[step] = vcat(K_hat, xi_hat)
            
            # d_res_hat = r_res_full ./ K_hat

            # 5d) simulate simplified model
            p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon_s, A_s)

            rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(),
                plot=false,
                cb=cb,
                full_or_simple=false
            )
            rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(),
                plot=false,
                cb=cb,
                species_specific_perturbation=false
            )
            
            before_persistence_S[step] = before_s
            after_persistence_S[step]  = after_s
            after_pulse_S[step] = after_pulse3
            rt_press_S[step]   = mean(filter(!isnan, rt_press2))
            rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
            rt_press_vector_S[step] = rt_press2
            rt_pulse_vector_S[step] = rt_pulse3
            collectivity_S[step] = psi_s
            resilience_S[step] = compute_resilience(B0, p_s)
            reactivity_S[step] = compute_reactivity(B0, p_s)
            # println("Step $step: ")
            # @show maximum(abs, K_hat - K_res)
            # @show maximum(abs, xi_hat - xi_cons)
            # @show maximum(abs, B0 - fixed)
            # @info "tau full is $(tau_full) and tau short is $(mean(filter(!isnan, 1.0 ./ (vcat(r_res,m_cons).*B0./vcat(K_hat,xi_hat)))))"

            D_s, M_s = compute_jacobian(B0, p_s)
            J_s = D_s * M_s
            J_diff_S[step] = norm(J_s - J_full)
            stable_S[step] = is_locally_stable(J_s)
            
            # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
            # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
            analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

            tau_S[step] = -1 ./ diag(J_s)
            mean_tau_S[step] = mean(-1 ./ diag(J_s))

            # min_delta_K_S[step] = fill(0.0, R)
            # for i in 1:R
            #     min_delta_K_S[step][i] = find_min_K_reduction(
            #         i, B0, p_s, tspan, tpert; cb
            #     )
            # end
            # min_delta_xi_S[step] = fill(0.0, C)
            # for i in 1:C
            #     min_delta_xi_S[step][i] = find_min_xi_increase(
            #         i, B0, p_s, tspan, tpert; cb
            #     )
            # end
            # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
            # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
        end
        
        for step in 4:11
            # -------------------------------------------------------------------------
            # step 17: re-randomise m_cons
            # -------------------------------------------------------------------------
            if step == 4
                A_s, epsilon_s = short_transform_for_ladder_step(1, copy(A), copy(epsilon))
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))

                # 1) check that A, epsilon are unchanged
                @assert all(A_s .== A)             "A was mutated on step $step!"
                @assert all(epsilon_s .== epsilon) "epsilon was mutated on step $step!"
                @assert isapprox(psi_s, psi; atol=1e-12)  "collectivity mismatch at step $step: psi=$psi, psi_s=$(psi_s)"

                # new m_cons vector ~ Normal(mean(m_cons), 0.2)
                mum = mean(m_cons)
                mum = min(mum*1.5, 1.0)
                m_cons_s = rand.(Normal(mum, 0.2), C)

                # parameters tuple with updated m_cons
                p_s = (R, C, m_cons_s, xi_cons, r_res, K_res, epsilon_s, A_s)

                # now exactly the same simulate & record as before:
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                species_specific_perturbation=false
                )

                K_Xi_S[step] = vcat(K_res, xi_cons)

                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(B0, p_s)
                reactivity_S[step] = compute_reactivity(B0, p_s)

                D_s, M_s = compute_jacobian(B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)

                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s)) 

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
            # -------------------------------------------------------------------------
            # step 18: re-randomise xi_cons
            # -------------------------------------------------------------------------
            elseif step == 5
                A_s, epsilon_s = short_transform_for_ladder_step(1, copy(A), copy(epsilon))
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))

                # 1) check that A, epsilon are unchanged
                @assert all(A_s .== A)             "A was mutated on step $step!"
                @assert all(epsilon_s .== epsilon) "epsilon was mutated on step $step!"
                @assert isapprox(psi_s, psi; atol=1e-12)  "collectivity mismatch at step $step: psi=$psi, psi_s=$(psi_s)"

                # new xi_cons vector ~ Normal(mean(xi_cons), 0.1*mean(xi_cons))
                muxi = mean(xi_cons)
                sigmaxi = 0.1 * muxi
                xi_cons_s = rand.(Normal(muxi, sigmaxi), C)

                # # recalculate equilibrium with new xi_cons but same K_hat
                # eq_s = calibrate_from_K_xi(xi_cons_s, K_res, epsilon_s, A_s)
                # R_eq_s, C_eq_s = eq_s[1:R], eq_s[R+1:end]

                # parameters tuple with updated xi_cons
                p_s = (R, C, m_cons, xi_cons_s, r_res, K_res, epsilon_s, A_s)

                # same simulate & record
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                species_specific_perturbation=false
                )

                K_Xi_S[step] = vcat(K_res, xi_cons_s)

                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(B0, p_s)
                reactivity_S[step] = compute_reactivity(B0, p_s)

                D_s, M_s = compute_jacobian(B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)

                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s))

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
            # -------------------------------------------------------------------------
            # step 19: re-randomise K_res
            # -------------------------------------------------------------------------
            elseif step == 6
                A_s, epsilon_s = short_transform_for_ladder_step(1, copy(A), copy(epsilon))
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))

                # 1) check that A, epsilon are unchanged
                @assert all(A_s .== A)             "A was mutated on step $step!"
                @assert all(epsilon_s .== epsilon) "epsilon was mutated on step $step!"
                @assert isapprox(psi_s, psi; atol=1e-12)  "collectivity mismatch at step $step: psi=$psi, psi_s=$(psi_s)"

                # new K_res vector ~ Normal(mean(K_res), 0.1*mean(K_res))
                muK = mean(K_res)
                sigmaK = 0.1 * muK
                K_res_s = rand.(Normal(muK, sigmaK), R)

                # # recalculate equilibrium with new K_res but same xi_hat
                # eq_s = calibrate_from_K_xi(xi_hat, K_res_s, epsilon_s, A_s)
                # R_eq_s, C_eq_s = eq_s[1:R], eq_s[R+1:end]

                # parameters tuple with updated K_res
                p_s = (R, C, m_cons, xi_cons, r_res, K_res_s, epsilon_s, A_s)

                # same simulate & record
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(), plot=false, cb=cb,
                species_specific_perturbation=false
                )

                K_Xi_S[step] = vcat(K_res_s, xi_cons)

                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(B0, p_s)
                reactivity_S[step] = compute_reactivity(B0, p_s)

                D_s, M_s = compute_jacobian(B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)
                
                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s))

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
            elseif step == 7
                
                A_s, epsilon_s = short_transform_for_ladder_step(2, copy(A), copy(epsilon))
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))

                meanR = mean(B0[1:R])
                meanC = mean(B0[R+1:S])

                C_eq_avg = fill(meanC, C)
                R_eq_avg = fill(meanR, R)
                new_B0 = vcat(R_eq_avg, C_eq_avg)
                
                # 5a) Recompute xi_hat
                xi_hat = zeros(C)
                for k in 1:C
                    i = R + k
                    # 1) feeding gains (A>0)
                    gain = sum(epsilon_s[i,j]*A_s[i,j]*R_eq_avg[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
                    # 2) predation losses: A_s[i,j]<0, but we need "+(-A)B"
                    loss = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                    # consumer eq: xi = B_i - gain - loss
                    xi_hat[k] = -C_eq_avg[k] + gain + loss
                end

                # 5b) Recompute K_hat
                K_hat = zeros(R)
                for i in 1:R
                    # resource eq uses A[i,j] (j consumer) directly:
                    drain = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0)
                    # K_i = B_i + ? A[i,j] B_j
                    K_hat[i] = abs(-R_eq_avg[i] + drain)
                end

                # 5c) Solve for new equilibrium
                eq = try
                        calibrate_from_K_xi(xi_hat, K_hat, epsilon_s, A_s)
                    catch err
                        @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                        continue
                    end
        
                R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
                # also guard against non-finite or non-positive solution
                # if any(!isfinite, eq) || any(x->x<=0, eq)
                #     @warn "Step $step: infeasible equilibrium (non-finite or =0)"
                #     continue
                # end
                
                # tau_S[step] = 1.0 ./ eq
                K_Xi_S[step] = vcat(K_hat, xi_hat)
                
                # d_res_hat = r_res_full ./ K_hat

                # 5d) simulate simplified model
                p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon_s, A_s)
                
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    species_specific_perturbation=false
                )
                
                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(new_B0, p_s)
                reactivity_S[step] = compute_reactivity(new_B0, p_s)

                D_s, M_s = compute_jacobian(new_B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)
                
                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s))

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])

            elseif step == 8
                
                A_s, epsilon_s = short_transform_for_ladder_step(3, copy(A), copy(epsilon))
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))

                meanR = mean(B0[1:R])
                meanC = mean(B0[R+1:S])

                C_eq_avg = fill(meanC, C)
                R_eq_avg = fill(meanR, R)
                new_B0 = vcat(R_eq_avg, C_eq_avg)
                
                # 5a) Recompute xi_hat
                xi_hat = zeros(C)
                for k in 1:C
                    i = R + k
                    # 1) feeding gains (A>0)
                    gain = sum(epsilon_s[i,j]*A_s[i,j]*R_eq_avg[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
                    # 2) predation losses: A_s[i,j]<0, but we need "+(-A)B"
                    loss = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                    # consumer eq: xi = B_i - gain - loss
                    xi_hat[k] = -C_eq_avg[k] + gain + loss
                end

                # 5b) Recompute K_hat
                K_hat = zeros(R)
                for i in 1:R
                    # resource eq uses A[i,j] (j consumer) directly:
                    drain = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0)
                    # K_i = B_i + ? A[i,j] B_j
                    K_hat[i] = abs(-R_eq_avg[i] + drain)
                end

                # 5c) Solve for new equilibrium
                eq = try
                        calibrate_from_K_xi(xi_hat, K_hat, epsilon_s, A_s)
                    catch err
                        @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                        continue
                    end
        
                R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
                # also guard against non-finite or non-positive solution
                # if any(!isfinite, eq) || any(x->x<=0, eq)
                #     @warn "Step $step: infeasible equilibrium (non-finite or =0)"
                #     continue
                # end
                
                K_Xi_S[step] = vcat(K_hat, xi_hat)
                
                # d_res_hat = r_res_full ./ K_hat

                # 5d) simulate simplified model
                p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon_s, A_s)
                
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    species_specific_perturbation=false
                )
                
                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(new_B0, p_s)
                reactivity_S[step] = compute_reactivity(new_B0, p_s)

                D_s, M_s = compute_jacobian(new_B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)
                
                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                J_diff_S[step] = norm(J_s - J_full)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s))

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
            elseif step in 9:11
                if step == 9
                    A_s = make_A(copy(A), R, conn, scen; IS = IS, B_term = B_term)
                elseif step == 10
                    A_s = make_A(copy(A), R, conn, :ER; IS = IS, B_term = B_term)
                elseif step == 11
                    A_s = make_A(copy(A), R, rand(Uniform(0.1, 0.9)), :ER; IS = IS, B_term = B_term)
                end
                
                epsilon_s = copy(epsilon)
                
                psi_s = compute_collectivity(copy(A_s), copy(epsilon_s))
                
                meanR = mean(B0[1:R])
                meanC = mean(B0[R+1:S])

                C_eq_avg = fill(meanC, C)
                R_eq_avg = fill(meanR, R)
                new_B0 = vcat(R_eq_avg, C_eq_avg)
                
                # 5a) Recompute xi_hat
                xi_hat = zeros(C)
                for k in 1:C
                    i = R + k
                    # 1) feeding gains (A>0)
                    gain = sum(epsilon_s[i,j]*A_s[i,j]*R_eq_avg[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
                    # 2) predation losses: A_s[i,j]<0, but we need "+(-A)B"
                    loss = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                    # consumer eq: xi = B_i - gain - loss
                    xi_hat[k] = -C_eq_avg[k] + gain + loss
                end

                # 5b) Recompute K_hat
                K_hat = zeros(R)
                for i in 1:R
                    # resource eq uses A[i,j] (j consumer) directly:
                    drain = sum(A_s[i,j]*C_eq_avg[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0)
                    # K_i = B_i + ? A[i,j] B_j
                    K_hat[i] = abs(-R_eq_avg[i] + drain)
                end

                # 5c) Solve for new equilibrium
                eq = try
                        calibrate_from_K_xi(xi_hat, K_hat, epsilon_s, A_s)
                    catch err
                        @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                        continue
                    end
        
                R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
                # also guard against non-finite or non-positive solution
                # if any(!isfinite, eq) || any(x->x<=0, eq)
                #     @warn "Step $step: infeasible equilibrium (non-finite or =0)"
                #     continue
                # end
                
                K_Xi_S[step] = vcat(K_hat, xi_hat)
                
                # d_res_hat = r_res_full ./ K_hat

                # 5d) simulate simplified model
                p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon_s, A_s)
                
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    full_or_simple=false
                )
                rt_pulse3, _, after_pulse3, _, _ = simulate_pulse_perturbation(
                    new_B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    species_specific_perturbation=false
                )
                
                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                after_pulse_S[step] = after_pulse3
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                rt_press_vector_S[step] = rt_press2
                rt_pulse_vector_S[step] = rt_pulse3
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(new_B0, p_s)
                reactivity_S[step] = compute_reactivity(new_B0, p_s)

                D_s, M_s = compute_jacobian(new_B0, p_s)
                J_s = D_s * M_s
                J_diff_S[step] = norm(J_s - J_full)
                stable_S[step] = is_locally_stable(J_s)
                
                # Rmed_S[step] = median_return_rate(J_s, B0; t=t0, n=Rmed_iterations)
                # ssp_rmed_S[step] = species_return_rates(J_s, B0; t=t0, n=Rmed_iterations)
                analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)

                J_diff_S[step] = norm(J_s - J_full)

                tau_S[step] = -1 ./ diag(J_s)
                mean_tau_S[step] = mean(-1 ./ diag(J_s))

                # min_delta_K_S[step] = fill(0.0, R)
                # for i in 1:R
                #     min_delta_K_S[step][i] = find_min_K_reduction(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # min_delta_xi_S[step] = fill(0.0, C)
                # for i in 1:C
                #     min_delta_xi_S[step][i] = find_min_xi_increase(
                #         i, new_B0, p_s, tspan, tpert; cb
                #     )
                # end
                # mean_min_delta_K_S[step] = mean(min_delta_K_S[step])
                # mean_min_delta_xi_S[step] = mean(min_delta_xi_S[step])
            end
        end

        step_pairs = collect(Iterators.flatten(
            ([
                Symbol("before_persistence_S$i") => before_persistence_S[i],
                Symbol("after_persistence_S$i") => after_persistence_S[i],
                Symbol("after_pulse_S$i") => after_pulse_S[i],
                Symbol("rt_press_S$i") => rt_press_S[i],
                Symbol("rt_pulse_S$i") => rt_pulse_S[i],
                Symbol("collectivity_S$i") => collectivity_S[i],
                Symbol("resilience_S$i") => resilience_S[i],
                Symbol("reactivity_S$i") => reactivity_S[i],

                Symbol("stable_S$i") => stable_S[i],
                
                # Symbol("Rmed_S$i") => Rmed_S[i],
                # Symbol("ssp_rmed_S$i") => ssp_rmed_S[i],
                Symbol("analytical_rmed_S$i") => analytical_rmed_S[i],
                
                Symbol("tau_S$i") => tau_S[i],
                Symbol("mean_tau_S$i") => mean_tau_S[i],

                Symbol("K_Xi_S$i") => K_Xi_S[i],
                Symbol("J_diff_S$i") => J_diff_S[i],
                
                # Symbol("min_delta_K_S$i") => min_delta_K_S[i],
                # Symbol("min_delta_xi_S$i") => min_delta_xi_S[i],
                # Symbol("mean_min_delta_K_S$i") => mean_min_delta_K_S[i],
                # Symbol("mean_min_delta_xi_S$i") => mean_min_delta_xi_S[i],

                Symbol("rt_press_vector_S$i") => rt_press_vector_S[i],
                Symbol("rt_pulse_vector_S$i") => rt_pulse_vector_S[i]
            ] for i in 1:11)
        ))

        rec = (
            conn=conn, IS=IS, scen=scen, delta =delta, epsi=epsi, m_val=m_val, g_val=g_val, ite =ite,
            before_persistence_full=before_full, after_persistence_full=after_persistence_full, after_pulse_full=after_pulse_full,
            rt_press_full=rt_press_full, rt_pulse_full=rt_pulse_full,
            collectivity_full=collectivity_full, resilience_full=resilience_full, reactivity_full=reactivity_full,
            # Rmed_full=Rmed_full, ssp_rmed_full=ssp_rmed_full
            analytical_rmed_full=analytical_rmed_full,
            mean_tau_full=mean_tau_full, tau_full=tau_full, J_full = J_full,
            J_diff_full=J_diff_full, J_full_norm=J_full_norm,
            rt_pulse_vector_full=rt_pulse_full_vector, rt_press_vector_full=rt_press_full_vector,
            # mean_min_delta_K_full = mean_min_delta_K_full, mean_min_delta_xi_full = mean_min_delta_xi_full,
            # min_delta_K_full = min_delta_K_full, min_delta_xi_full = min_delta_xi_full,
            step_pairs...,  # Properly flattened pairs
            p_final = p,
            R_eq = R_eq,
            C_eq = C_eq,
            B0 = B0,
            m_cons = m_cons, r_res = r_res,
            K_Xi_full=original_k_xi,
        )

        lock(locki) do
            push!(results, rec)
        end
    end
    @info "Finished computing $number_of_combinations combinations"
    return DataFrame(results)
end

# --------------------------------------------------------------------------------
# 9) Run it
# --------------------------------------------------------------------------------
R = short_ComputingLadder(
    50, 20;
    conn_vals=0.01:0.04:0.9,
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    IS_vals_B_term=[0.1, 1.0],
    scenarios=[:ER],
    delta_vals=[0.1, 0.3, 0.5, 0.75, 0.01, 0.9],
    eps_scales=[1.0, 0.5, 0.1],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    tspan=(0.,500.0), tpert=250.0,
    number_of_combinations = 5000,
    B_term = false,
    iterations=1,
    Rmed_iterations=1
)

@info "we reached here"
# serialize("Ladder/Outputs/T.jls", T)

T = deserialize("ThePaper/Ladder/Outputs/T.jls")
F = deserialize("ThePaper/Ladder/Outputs/exploring_min_extinction_5000_withsssp_rmed.jls")
