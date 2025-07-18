
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
    if isempty(d)
        return NaN
    end
    min_d = minimum(d)
    offs = [A[i,j] for i in 1:size(A,1), j in 1:size(A,1) if i != j]
    if isempty(offs)
        return NaN
    end
    sigma = std(offs)
    return sigma / min_d
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
function checking_for_jacobians(
    S::Int=50, C::Int=20;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    IS_vals_B_term=[0.1,1.0],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[1.0],
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
        !ok && continue
        # if !all(isapprox.(B0, fixed, atol=1e-3))
        #     @warn "B0 is not close to fixed point: B0 =  $(B0), and fixed = $(fixed)"
        #     continue
        # end
        D, M = compute_jacobian(B0, p)
        J_full = D * M
        !is_locally_stable(J_full) && continue
        S_full = sum(B0 .> 0.0; init=0.0)
        # if S_full < 50 && is_locally_stable(J_full)
        #     @error("Hey, we got here!!")
        # end
        # resilience_full is -max Re(eig(J)), but compute_resilience returns max Re
        old_epsilon = copy(epsilon)
        resilience_full = compute_resilience(B0, p; extinct_species=false)
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside compute_resilience!"
        reactivity_full =  compute_reactivity(B0, p; extinct_species=false)
        
        original_k_xi = vcat(K_res, xi_cons)
        # but you can also sweep t=0.1, 1, 10, etc.
        t0 = 0.01

        # full-model median return rate
        # Rmed_full = median_return_rate(J_full, fixed; t=t0, n=Rmed_iterations)
        # ssp_rmed_full = species_return_rates(J_full, fixed; t=t0, n=Rmed_iterations)
        ssp_analytical_rmed_full = analytical_species_return_rates(J_full; t=0.001)
        analytical_rmed_full = analytical_median_return_rate(J_full; t=0.001)

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
        # rt_pulse_full_vector = rt_pulse
        # rt_press_full_vector = rt_press
        collectivity_full = psi 
        tau_full = -1 ./ diag(J_full)
        mean_tau_full = mean(tau_full)

        inverse_tau_full = -diag(J_full)
        mean_inverse_tau_full = mean(inverse_tau_full)

        SL_full = compute_SL(A, vcat(K_res,xi_cons))
        mean_SL_full = mean(compute_SL(A, vcat(K_res,xi_cons)))
        # tau_full = compute_SL(A, vcat(K_res,xi_cons))
        # mean_tau_full = mean(compute_SL(A, vcat(K_res,xi_cons)))
        sigma_over_min_d_full = sigma_over_min_d(A, J_full)

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
        S_S = Dict(i => NaN for i in 1:4)
        collectivity_S = Dict(i => NaN for i in 1:4)
        resilience_S  = Dict(i=>NaN for i in 1:4)
        reactivity_S  = Dict(i=>NaN for i in 1:4)
        # stable_S    = Dict(i=>NaN for i in 1:4)
        
        # Rmed_S    = Dict(i=>NaN for i in 1:4)
        # ssp_rmed_S = Dict(i => Float66[] for i in 1:4)
        analytical_rmed_S = Dict(i => NaN for i in 1:4)
        ssp_analytical_rmed_S = Dict(i => Float64[] for i in 1:4)
        
        tau_S = Dict(i => Float64[] for i in 1:4)
        mean_tau_S = Dict(i => NaN for i in 1:4)

        inverse_tau_S = Dict(i => Float64[] for i in 1:4)
        mean_inverse_tau_S = Dict(i => NaN for i in 1:4)

        SL_S = Dict(i => Float64[] for i in 1:4)
        mean_SL_S = Dict(i => NaN for i in 1:4)
        # K_Xi_S = Dict(i => Float64[] for i in 1:6)
        # J_diff_S = Dict(i => NaN for i in 1:6)
        # min_delta_K_S = Dict(i => Float64[] for i in 1:6)
        # min_delta_xi_S = Dict(i => Float64[] for i in 1:6)
        # mean_min_delta_K_S = Dict(i => NaN for i in 1:6)
        # mean_min_delta_xi_S = Dict(i => NaN for i in 1:6)
        # rt_press_vector_S = Dict(i => Float64[] for i in 1:6)
        # rt_pulse_vector_S = Dict(i => Float64[] for i in 1:6)
        sigma_over_min_d_S = Dict(i => NaN for i in 1:4)
        @info "Running ladder"

        # original equilibrium abundances
        # R_eq_full, C_eq_full = B0[1:R], B0[R+1:S] # B0 is the simulated equilibrium
        local R_eq_full, C_eq_full = B0[1:R], B0[R+1:S] # B0 is the calibrated equilibrium

        ###########################################################################################
        ###########################################################################################
        ###########################################################################################
        ########################## SIMPLIFIED MODEL STEPS #########################################
        ###########################################################################################
        ###########################################################################################
        for step in 1:3
            
            A_s = copy(A)
            
            # 5a) Recompute xi_hat
            xi_hat = xi_cons
            
            K_hat = K_res

            # 5c) Solve for new equilibrium
            eq = try
                    calibrate_from_K_xi(xi_hat, K_hat, epsilon, A_s)
                catch err
                    @warn "Step $step: equilibrium solve failed (singular or NaNs)" 
                    continue
                end
    
            R_eq_s, C_eq_s = eq[1:R], eq[R+1:S]
            
            p_s = (R, C, m_cons, xi_hat, r_res, K_hat, epsilon, A_s)
            # if step == 5
            #     println(
            #         "r_res: ", length(r_res), " m_cons: ", length(m_cons), " K_hat: ", length(K_hat), " xi_hat: ", length(xi_hat), " C: ", C, " R: ", R, "R_eq_full: ", length(R_eq_full), " C_eq_full: ", length(C_eq_full)
            #     )
            # end
            ok, new_B0 = survives!(eq, p_s; cb=cb)
            S_S[step] = sum(new_B0 .> 0.0; init=0.0)
            
            # after_pulse_S[step] = after_pulse3
            # rt_press_S[step]   = mean(filter(!isnan, rt_press2))
            # rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
            # rt_press_vector_S[step] = rt_press2
            # rt_pulse_vector_S[step] = rt_pulse3
            collectivity_S[step] = compute_collectivity(A_s, epsilon)
            
            # ——— reshuffle J_s according to `step` ———
            # J_s :: AbstractMatrix{<:Real}
            # R      = number of resources (first R rows/cols)
            # step   ∈ {1,2,3,4}
            J_s = copy(J_full)    
            S = size(J_s,1)
            resources = 1:R
            consumers = (R+1):S

            if step == 1
                # 1) shuffle only the diagonal, separately for R and C
                vals = [J_s[i,i] for i in resources]
                shuffle!(vals)
                for (i,v) in zip(resources, vals)
                    J_s[i,i] = v
                end

                vals = [J_s[i,i] for i in consumers]
                shuffle!(vals)
                for (i,v) in zip(consumers, vals)
                    J_s[i,i] = v
                end

            elseif step == 2
                # 2) shuffle only off-diagonals, separately for R and C rows
                # resources
                pos = [(i,j) for i in resources, j in 1:S if j != i]
                vals = [J_s[i,j] for (i,j) in pos]
                shuffle!(vals)
                for ((i,j), v) in zip(pos, vals)
                    J_s[i,j] = v
                end

                # consumers
                pos = [(i,j) for i in consumers, j in 1:S if j != i]
                vals = [J_s[i,j] for (i,j) in pos]
                shuffle!(vals)
                for ((i,j), v) in zip(pos, vals)
                    J_s[i,j] = v
                end

            elseif step == 3
                # 3) shuffle both diag & off-diag, but still separately by R vs C
                # resources
                pos = [(i,j) for i in resources, j in 1:S]
                vals = [J_s[i,j] for (i,j) in pos]
                shuffle!(vals)
                for ((i,j), v) in zip(pos, vals)
                    J_s[i,j] = v
                end

                # consumers
                pos = [(i,j) for i in consumers, j in 1:S]
                vals = [J_s[i,j] for (i,j) in pos]
                shuffle!(vals)
                for ((i,j), v) in zip(pos, vals)
                    J_s[i,j] = v
                end

            elseif step == 4
                # 4) shuffle every entry in J_s together
                pos = [(i,j) for i in 1:S, j in 1:S]
                vals = [J_s[i,j] for (i,j) in pos]
                shuffle!(vals)
                for ((i,j), v) in zip(pos, vals)
                    J_s[i,j] = v
                end

            else
                error("Invalid step = $step; must be 1–4")
            end
    

            # tau_S[step] = compute_SL(A_s, vcat(K_res,xi_cons))
            # mean_tau_S[step] = mean(compute_SL(A_s, vcat(K_res,xi_cons)))
            
            # println("Step $step: ")
            # @show maximum(abs, K_hat - K_res)
            # @show maximum(abs, xi_hat - xi_cons)
            # @show maximum(abs, B0 - fixed)
            # @info "tau full is $(tau_full) and tau short is $(mean(filter(!isnan, 1.0 ./ (vcat(r_res,m_cons).*new_B0./vcat(K_hat,xi_hat)))))"

            # D_s, M_s = compute_jacobian(new_B0, p_s)
            # J_s = D_s * M_s
            extant = findall(bi -> bi > EXTINCTION_THRESHOLD, new_B0)
            J_s_sub = J_s[extant, extant]
            # J_diff_S[step] = norm(J_s - J_full)
            # stable_S[step] = is_locally_stable(J_s)
            
            sigma_over_min_d_S[step] = sigma_over_min_d(A_s, J_s_sub)

            # Rmed_S[step] = median_return_rate(J_s, new_B0; t=t0, n=Rmed_iterations)
            # ssp_rmed_S[step] = species_return_rates(J_s, new_B0; t=t0, n=Rmed_iterations)
            ssp_analytical_rmed_S[step] = analytical_species_return_rates(J_s; t=0.1)
            analytical_rmed_S[step] = analytical_median_return_rate(J_s; t=0.1)

            # 3) RESILIENCE
            ev = eigvals(J_s)
            resilience_S[step] = maximum(real.(ev))
            # 4) REACTIVITY
            J_sym = (J_s + J_s') / 2
            ev_sym = eigvals(J_sym)
            if isempty(ev_sym) && return 
                reactivity_S[step] = NaN
            else
                reactivity_S[step] = maximum(real.(ev_sym))
            end

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
                Symbol("S_S$i") => S_S[i],
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
            after_persistence_full=after_persistence_full,
            # before_persistence_full=before_full, after_persistence_full=after_persistence_full, after_pulse_full=after_pulse_full,
            # rt_press_full=rt_press_full, rt_pulse_full=rt_pulse_full,
            S_full=S_full, collectivity_full=collectivity_full, 
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
R = checking_for_jacobians(
    50, 20;
    conn_vals=0.01:0.01:1.0,
    IS_vals=[0.0001, 0.001, 0.01, 0.1, 1.0, 2.0],
    IS_vals_B_term=[0.1, 1.0],
    scenarios=[:ER],
    delta_vals=[0.5, 0.3, 0.9, 0.99], #[0.1, 0.3, 0.5, 0.75, 0.01, 0.9],
    eps_scales=[1.0, 0.5, 0.1],
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    tspan=(0.,500.0), tpert=250.0,
    number_of_combinations = 10000,
    B_term = false,
    iterations=5,
    Rmed_iterations=5,
    pareto_exponents=[1.0],#[1.0,1.5,2.0,3.0,4.0,5.0],
    pareto_minimum_degrees=[1.0],#[1.0,2.0,3.0,4.0,5.0,6.0],
    mod_gammas=[1.0],#[1.0,2.0,3.0,5.0,10.0]
)

desired = [
  :conn, :IS, :scen, :delta, :epsi, :m_val, :g_val, :ite, :pex, :p_min_deg, :mod_gamma,
  :S_full, :resilience_full, :reactivity_full, :collectivity_full, :tau_full, :mean_tau_full, :sigma_over_min_d_full, :SL_full, :mean_SL_full, :inverse_tau_full, :mean_inverse_tau_full, :analytical_rmed_full, :ssp_analytical_rmed_full,
  :S_S1,:resilience_S1, :reactivity_S1, :collectivity_S1, :tau_S1, :mean_tau_S1, :sigma_over_min_d_S1, :SL_S1, :mean_SL_S1, :inverse_tau_S1, :mean_inverse_tau_S1, :analytical_rmed_S1, :ssp_analytical_rmed_S1,
  :S_S2,:resilience_S2, :reactivity_S2, :collectivity_S2, :tau_S2, :mean_tau_S2, :sigma_over_min_d_S2, :SL_S2, :mean_SL_S2, :inverse_tau_S2, :mean_inverse_tau_S2, :analytical_rmed_S2, :ssp_analytical_rmed_S2,
  :S_S3,:resilience_S3, :reactivity_S3, :collectivity_S3, :tau_S3, :mean_tau_S3, :sigma_over_min_d_S3, :SL_S3, :mean_SL_S3, :inverse_tau_S3, :mean_inverse_tau_S3, :analytical_rmed_S3, :ssp_analytical_rmed_S3,
  :S_S4,:resilience_S4, :reactivity_S4, :collectivity_S4, :tau_S4, :mean_tau_S4, :sigma_over_min_d_S4, :SL_S4, :mean_SL_S4, :inverse_tau_S4, :mean_inverse_tau_S4, :analytical_rmed_S4, :ssp_analytical_rmed_S4,
  :S_S5,:resilience_S5, :reactivity_S5, :collectivity_S5, :tau_S5, :mean_tau_S5, :sigma_over_min_d_S5, :SL_S5, :mean_SL_S5, :inverse_tau_S5, :mean_inverse_tau_S5, :analytical_rmed_S5, :ssp_analytical_rmed_S5
]

G = R[!, desired]

serialize("checking_changing_groups_ALL1000000_plussubgroups.jls", R)