using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")

# Compute SL = B ./ K with B = (I - A) \ K
function compute_SL(A::Matrix{Float64}, K::Vector{Float64})
    B = (I - A) \ K
    return B ./ K
end
# Extract sigma/min(d)
function sigma_over_min_d(A, J)
    d = -diag(J)
    min_d = minimum(d)
    offs = [A[i,j] for i in 1:size(A,1), j in 1:size(A,1) if i!=j]
    sigma = std(offs)
    return sigma/min_d
end

# rewire in place
function rewire_A!(A, mask, sigma)
    A[mask] .= randn(sum(mask)) .* sigma
end

# ----------------------------------------------------------------
# 5) Stability & survival checks
# -------------------------------------------------------------------------------
function is_locally_stable(J::AbstractMatrix)
    if any(!isfinite, J)
        return false
    end
    lambda = eigvals(J)
    maximum(real.(lambda)) <= 0 
end

function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false,sol.u[end]) : (all(sol.u[end] .> EXTINCTION_THRESHOLD),sol.u[end])
end

# --------------------------------------------------------------------------------
# Main sweep: fix S & C, vary structure, collect full & ladder pers.
# --------------------------------------------------------------------------------
function NF_checking(
    S::Int=50, C::Int=20;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    IS_vals_B_term=[0.1,1.0],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[0.1],
    delta_vals=[1.0,3.0],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100,
    B_term = false,
    B_term_IS = [0.1],
    iterations=1,
    Rmed_iterations=200,
    pareto_exponents=[1.25,1.75,2.0,3.0,4.0,5.0],
    pareto_minimum_degrees=[1.0],
    mod_gammas = [1.0,2.0,3.0,5.0,10.0],
    minimum_n_of_consumers=5
)
    R = S - C
    # A = zeros(S,S)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S,EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        conn_vals,IS_vals,scenarios,delta_vals,eps_scales,mortality_vals,growth_vals,pareto_exponents,pareto_minimum_degrees,mod_gammas,1:iterations
    ))
    println("Number of combinations: ", length(combos))
    
    @threads for (conn,IS,scen,delta,epsi,m_val,g_val,pex,p_min_deg,mod_gamma,ite) in sample(combos, min(length(combos), number_of_combinations); replace=false)
        
        # 1) build A & epsilon
        local A = make_A(zeros(S,S),R,conn,scen; IS=IS,B_term=B_term,B_term_IS=B_term_IS[1],pareto_exponent=pex,pareto_minimum_degree=p_min_deg,mod_gamma=mod_gamma)
        # if scen == :PL
        #     local A = degree_scaled_A(A, R, scen)
        # end

        local epsilon = clamp.(randn(S,S).*epsi,0,1)
        # @info "BP1"
        # --- NEW: compute collectivity phi ---
        local psi = compute_collectivity(copy(A), copy(epsilon))
        # backup_A = deepcopy(A)
        # backup_epsilon = deepcopy(epsilon)
        
        #### #TODO DEACTIVATE THIS PART IF YOU WANT TO GO BACK TO ENSURING FULL FEASIBILITY #####
        K_res, xi_cons = sample_Xi_K_from_A(A, epsilon, R;
                            scale_factor = 1.0)
        
        eq = try
            calibrate_from_K_xi(xi_cons, K_res, epsilon, A)
            catch err
                @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                continue
        end
        # @info "BP2"
        R_eq, C_eq = eq[1:R], eq[R+1:S]

        # set up parameters & run the rest of your pipeline exactly as before
        # m_cons = fill(m_val, C)
        m_cons = abs.(rand(Normal(m_val, 0.2), C))
        # d_res  = ones(R)
        # r_res  = fill(g_val, R)
        r_res  = abs.(rand(Normal(g_val, 0.2), R))

        p     = (R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)
        fixed = vcat(R_eq, C_eq)
        fixed[fixed .< 0.0] .= 0.0
        # @info "BP3"
        # 3) stability & survival
        ok, B0 = survives!(fixed, p; cb=cb)
        R_eq_full, C_eq_full = B0[1:R], B0[R+1:S]
        realised_S_full = sum(B0 .> 0.0)
        
        if count(x->x>EXTINCTION_THRESHOLD, B0[31:50]) < minimum_n_of_consumers
            continue
        end
        # println("The initial abundances after cleaning are: \n R_B0=$(B0[1:R]) \n C_B0=$(B0[R+1:S])")
        # println("In step0: \n xi_cons=$xi_cons \n K_res=$K_res")

        # println("ok: $ok, B0: $B0")
        # !ok && continue
        # if !all(isapprox.(B0, fixed, atol=1e-3))
        #     @warn "B0 is not close to fixed point: B0 =  $(B0), and fixed = $(fixed)"
        #     # continue
        # end
        # error("The error is line 125")
        D, M = compute_jacobian(B0, p)
        J_full = D * M
        extant = findall(bi -> bi > EXTINCTION_THRESHOLD, B0)
        J_sub = J_full[extant, extant]
        if !is_locally_stable(J_sub) 
            @info "J_full is not locally stable"
            continue
        end
        # @info "BP5"
        # resilience_full is -max Re(eig(J)), but compute_resilience returns max Re
        old_epsilon = copy(epsilon)
        resilience_full = compute_resilience(B0, p; extinct_species = true)
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside compute_resilience!"
        reactivity_full =  compute_reactivity(B0, p; extinct_species = true)
        
        original_k_xi = vcat(K_res, xi_cons)
        # but you can also sweep t=0.1, 1, 10, etc.
        t0 = 0.01

        # full-model median return rate
        # Rmed_full = median_return_rate(J_full, fixed; t=t0, n=Rmed_iterations)
        # ssp_rmed_full = species_return_rates(J_full, fixed; t=t0, n=Rmed_iterations)
        ssp_analytical_rmed_full = analytical_species_return_rates(J_full; t=0.001)
        analytical_rmed_full = analytical_median_return_rate(J_full; t=0.001)
        # @info "BP6"
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
        # @info "BP7"
        # if manual_before_persistence != before_full || !isone(manual_before_persistence)
        #     error("manual_before_persistence = $(manual_before_persistence), before_full = $(before_full)")
        #     continue
        # end

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
        # rt_pulse_full_vector = rt_pulse
        # rt_press_full_vector = rt_press
        collectivity_full = psi

        tau_full = diag(J_full) .|> x -> x == 0.0 ? 0.0 : -1 / x
        mean_tau_full = mean(filter(!iszero, tau_full))

        inverse_tau_full = -diag(J_full)
        mean_inverse_tau_full = mean(inverse_tau_full)

        SL_full = compute_SL(A, vcat(K_res,xi_cons))
        mean_SL_full = mean(compute_SL(A, vcat(K_res,xi_cons)))
        # tau_full = compute_SL(A, vcat(K_res,xi_cons))
        # mean_tau_full = mean(compute_SL(A, vcat(K_res,xi_cons)))
        sigma_over_min_d_full = sigma_over_min_d(A, J_sub)

        # J_diff_full = norm(J_full - J_full)
        # J_full_norm = norm(J_full)
        
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
        # before_persistence_S = Dict(i => NaN for i in 1:4)
        # after_persistence_S  = Dict(i => NaN for i in 1:4)
        # after_pulse_S = Dict(i => NaN for i in 1:4)
        # rt_press_S   = Dict(i => NaN for i in 1:4)
        # rt_pulse_S   = Dict(i => NaN for i in 1:4)
        realised_S_S = Dict(i => NaN for i in 1:4)

        collectivity_S = Dict(i => NaN for i in 1:4)
        resilience_S  = Dict(i=>NaN for i in 1:4)
        reactivity_S  = Dict(i=>NaN for i in 1:4)
        # stable_S    = Dict(i=>NaN for i in 1:4)
        
        # Rmed_S    = Dict(i=>NaN for i in 1:4)
        # ssp_rmed_S = Dict(i => Float64[] for i in 1:4)
        analytical_rmed_S = Dict(i => NaN for i in 1:4)
        ssp_analytical_rmed_S = Dict(i => Float64[] for i in 1:4)
        
        tau_S = Dict(i => Float64[] for i in 1:4)
        mean_tau_S = Dict(i => NaN for i in 1:4)

        inverse_tau_S = Dict(i => Float64[] for i in 1:4)
        mean_inverse_tau_S = Dict(i => NaN for i in 1:4)

        SL_S = Dict(i => Float64[] for i in 1:4)
        mean_SL_S = Dict(i => NaN for i in 1:4)
        # K_Xi_S = Dict(i => Float64[] for i in 1:4)
        # J_diff_S = Dict(i => NaN for i in 1:4)
        # min_delta_K_S = Dict(i => Float64[] for i in 1:4)
        # min_delta_xi_S = Dict(i => Float64[] for i in 1:4)
        # mean_min_delta_K_S = Dict(i => NaN for i in 1:4)
        # mean_min_delta_xi_S = Dict(i => NaN for i in 1:4)
        # rt_press_vector_S = Dict(i => Float64[] for i in 1:4)
        # rt_pulse_vector_S = Dict(i => Float64[] for i in 1:4)
        sigma_over_min_d_S = Dict(i => NaN for i in 1:4)
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
        for step in 1:4
            
            # K_Xi_S[step] = vcat(K_res, xi_cons)
            A_s = copy(A)
            if step==1
                A_s = make_A(
                    A_s,R,conn,scen;
                    IS=IS,pareto_exponent=pex,pareto_minimum_degree=p_min_deg,
                    mod_gamma=mod_gamma,
                    B_term=B_term,B_term_IS=B_term_IS[1]
                )
                # if scen == :PL
                #     A_s = degree_scaled_A(A_s, R, scen)
                # end
            elseif step==2
                A_s = make_A(
                    A_s,R,min(conn*2, 1.0),scen;
                    IS=IS,pareto_exponent=pex,pareto_minimum_degree=p_min_deg,
                    mod_gamma=mod_gamma,
                    B_term=B_term,B_term_IS=B_term_IS[1]
                )
                # if scen == :PL
                #     A_s = degree_scaled_A(A_s, R, scen)
                # end
            elseif step==3
                A_s = make_A(
                    A_s,R,conn,scen;
                    IS=IS*2,pareto_exponent=pex,pareto_minimum_degree=p_min_deg,
                    mod_gamma=mod_gamma,
                    B_term=B_term,B_term_IS=B_term_IS[1]
                )
                # if scen == :PL
                #     A_s = degree_scaled_A(A_s, R, scen)
                # end
            elseif step==4
                A_s = make_A(
                    A_s,R,min(conn*2, 1.0),scen;
                    IS=IS*2,pareto_exponent=pex,pareto_minimum_degree=p_min_deg,
                    mod_gamma=mod_gamma,
                    B_term=B_term,B_term_IS=B_term_IS[1]
                )
                # if scen == :PL
                #     A_s = degree_scaled_A(A_s, R, scen)
                # end
            end
            
            # 5a) Recompute xi_hat
            xi_hat = zeros(C)
            for k in 1:C
                i = R + k
                # 1) feeding gains (A>0)
                gain = sum(epsilon[i,j]*A_s[i,j]*R_eq_full[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
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
                    calibrate_from_K_xi(xi_hat, K_hat, epsilon, A_s)
                catch err
                    @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                    continue
                end
    
            R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
            eq[eq .< 0.0] .= 0.0
            p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon, A_s)

            ok, new_B0 = survives!(eq, p_s; cb=cb)
            realised_S_S[step] = sum(new_B0 .> 0.0)
              
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
            
            # before_persistence_S[step] = before_s
            # after_persistence_S[step]  = after_s
            # after_pulse_S[step] = after_pulse3
            # rt_press_S[step]   = mean(filter(!isnan, rt_press2))
            # rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
            # rt_press_vector_S[step] = rt_press2
            # rt_pulse_vector_S[step] = rt_pulse3
            collectivity_S[step] = compute_collectivity(A_s, epsilon)
            try 
                resilience_S[step] = compute_resilience(new_B0, p_s; extinct_species = true)
            catch err
                @warn "Step $step: resilience computation failed (singular or NaNs)"
                resilience_S[step] = NaN
            end
            try
               reactivity_S[step] = compute_reactivity(new_B0, p_s; extinct_species = true) 
            catch
                @warn "Step $step: reactivity computation failed (singular or NaNs)"
                reactivity_S[step] = NaN
            end

            # tau_S[step] = compute_SL(A_s, vcat(K_res,xi_cons))
            # mean_tau_S[step] = mean(compute_SL(A_s, vcat(K_res,xi_cons)))
            
            # println("Step $step: ")
            # @show maximum(abs, K_hat - K_res)
            # @show maximum(abs, xi_hat - xi_cons)
            # @show maximum(abs, B0 - fixed)
            # @info "tau full is $(tau_full) and tau short is $(mean(filter(!isnan, 1.0 ./ (vcat(r_res,m_cons).*new_B0./vcat(K_hat,xi_hat)))))"

            D_s, M_s = compute_jacobian(new_B0, p_s)
            J_s = D_s * M_s
            extant = findall(bi -> bi > EXTINCTION_THRESHOLD, new_B0)
            J_s_sub = J_s[extant, extant]

            # J_diff_S[step] = norm(J_s - J_full)
            # stable_S[step] = is_locally_stable(J_s)
            
            sigma_over_min_d_S[step] = sigma_over_min_d(A_s, J_s_sub)

            # Rmed_S[step] = median_return_rate(J_s, new_B0; t=t0, n=Rmed_iterations)
            # ssp_rmed_S[step] = species_return_rates(J_s, new_B0; t=t0, n=Rmed_iterations)
            ssp_analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.001)
            analytical_rmed_S[step] = analytical_median_return_rate(J_s; t=0.001)

            tau_S[step] = diag(J_s) .|> x -> x == 0.0 ? 0.0 : -1 / x
            mean_tau_S[step] = mean(filter(!iszero, tau_S[step]))

            inverse_tau_S[step] = -diag(J_s)
            mean_inverse_tau_S[step] = mean(inverse_tau_S[step])

            SL_S[step] = compute_SL(A_s, vcat(K_hat,xi_hat))
            mean_SL_S[step] = mean(compute_SL(A_s, vcat(K_hat,xi_hat)))

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

        step_pairs = collect(Iterators.flatten(
            ([
                # Symbol("before_persistence_S$i") => before_persistence_S[i],
                # Symbol("after_persistence_S$i") => after_persistence_S[i],
                # Symbol("after_pulse_S$i") => after_pulse_S[i],
                # Symbol("rt_press_S$i") => rt_press_S[i],
                # Symbol("rt_pulse_S$i") => rt_pulse_S[i],
                Symbol("realised_S_S$i") => realised_S_S[i],
                Symbol("collectivity_S$i") => collectivity_S[i],
                Symbol("resilience_S$i") => resilience_S[i],
                Symbol("reactivity_S$i") => reactivity_S[i],

                # Symbol("stable_S$i") => stable_S[i],
                
                # Symbol("Rmed_S$i") => Rmed_S[i],
                # Symbol("ssp_rmed_S$i") => ssp_rmed_S[i],
                Symbol("analytical_rmed_S$i") => analytical_rmed_S[i],
                Symbol("ssp_analytical_rmed_S$i") => ssp_analytical_rmed_S[i],
                
                Symbol("tau_S$i") => tau_S[i],
                Symbol("mean_tau_S$i") => mean_tau_S[i],

                Symbol("inverse_tau_S$i") => inverse_tau_S[i],
                Symbol("mean_inverse_tau_S$i") => mean_inverse_tau_S[i],

                Symbol("SL_S$i") => SL_S[i],
                Symbol("mean_SL_S$i") => mean_SL_S[i],

                # Symbol("K_Xi_S$i") => K_Xi_S[i],
                # Symbol("J_diff_S$i") => J_diff_S[i],
                
                # Symbol("min_delta_K_S$i") => min_delta_K_S[i],
                # Symbol("min_delta_xi_S$i") => min_delta_xi_S[i],
                # Symbol("mean_min_delta_K_S$i") => mean_min_delta_K_S[i],
                # Symbol("mean_min_delta_xi_S$i") => mean_min_delta_xi_S[i],

                # Symbol("rt_press_vector_S$i") => rt_press_vector_S[i],
                # Symbol("rt_pulse_vector_S$i") => rt_pulse_vector_S[i],

                Symbol("sigma_over_min_d_S$i") => sigma_over_min_d_S[i]
            ] for i in 1:4)
        ))

        rec = (
            conn=conn, IS=IS, scen=scen, delta =delta, epsi=epsi, m_val=m_val, g_val=g_val, ite =ite,
            pex=pex, p_min_deg=p_min_deg, mod_gamma=mod_gamma,
            # before_persistence_full=before_full, after_persistence_full=after_persistence_full, after_pulse_full=after_pulse_full,
            # rt_press_full=rt_press_full, rt_pulse_full=rt_pulse_full,
            realised_S_full=realised_S_full, collectivity_full=collectivity_full, 
            resilience_full=resilience_full, reactivity_full=reactivity_full,
            # Rmed_full=Rmed_full, ssp_rmed_full=ssp_rmed_full, 
            analytical_rmed_full=analytical_rmed_full, ssp_analytical_rmed_full=ssp_analytical_rmed_full,
            sigma_over_min_d_full,
            tau_full=tau_full, inverse_tau_full=inverse_tau_full, 
            mean_tau_full=mean_tau_full, mean_inverse_tau_full=mean_inverse_tau_full,
            SL_full=SL_full, mean_SL_full=mean_SL_full,
            # J_full = J_full,
            # J_diff_full=J_diff_full, J_full_norm=J_full_norm,
            # rt_pulse_vector_full=rt_pulse_full_vector, rt_press_vector_full=rt_press_full_vector,
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
O = NF_checking(
    50, 20;
    conn_vals=0.01:0.02:0.5,
    IS_vals=[0.001, 0.01, 0.1, 1.0, 2.0],
    IS_vals_B_term=[0.1, 1.0],
    scenarios=[:ER],
    delta_vals=[MAX_DELTA], #[0.1, 0.3, 0.5, 0.75, 0.01, 0.9],
    eps_scales=[1.0, 0.5, 0.1],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    tspan=(0.,500.0), tpert=250.0,
    number_of_combinations = 10000,
    B_term = false,
    iterations=10,
    Rmed_iterations=10,
    pareto_exponents=[1.25], #[1.25,1.75,2.0,3.0,4.0,5.0],
    pareto_minimum_degrees=[1.0], #[1.0,2.0,3.0,4.0,5.0,6.0],
    mod_gammas=[1.0], #[1.0,2.0,3.0,5.0,10.0]
    minimum_n_of_consumers=5
)

desired = [
  :conn, :IS, :scen, :delta, :epsi, :m_val, :g_val, :ite, :pex, :p_min_deg, :mod_gamma,
  :realised_S_full, :resilience_full, :reactivity_full, :collectivity_full, :tau_full, :mean_tau_full, :sigma_over_min_d_full, :SL_full, :mean_SL_full, :inverse_tau_full, :mean_inverse_tau_full, :analytical_rmed_full, :ssp_analytical_rmed_full,
  :realised_S_S1, :resilience_S1, :reactivity_S1, :collectivity_S1, :tau_S1, :mean_tau_S1, :sigma_over_min_d_S1, :SL_S1, :mean_SL_S1, :inverse_tau_S1, :mean_inverse_tau_S1, :analytical_rmed_S1, :ssp_analytical_rmed_S1,
  :realised_S_S2, :resilience_S2, :reactivity_S2, :collectivity_S2, :tau_S2, :mean_tau_S2, :sigma_over_min_d_S2, :SL_S2, :mean_SL_S2, :inverse_tau_S2, :mean_inverse_tau_S2, :analytical_rmed_S2, :ssp_analytical_rmed_S2,
  :realised_S_S3, :resilience_S3, :reactivity_S3, :collectivity_S3, :tau_S3, :mean_tau_S3, :sigma_over_min_d_S3, :SL_S3, :mean_SL_S3, :inverse_tau_S3, :mean_inverse_tau_S3, :analytical_rmed_S3, :ssp_analytical_rmed_S3,
  :realised_S_S4, :resilience_S4, :reactivity_S4, :collectivity_S4, :tau_S4, :mean_tau_S4, :sigma_over_min_d_S4, :SL_S4, :mean_SL_S4, :inverse_tau_S4, :mean_inverse_tau_S4, :analytical_rmed_S4, :ssp_analytical_rmed_S4
]

G = O[!, desired]

serialize("NF_checking10000ER.jls", O)