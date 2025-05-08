function feasibility_search(
    S_vals, conn_vals, C_ratios, IS_vals,
    d_vals, mortality_vals, epsilon_vals, delta_vals;
    degree_distribution_types = [:ER, :PL, :MOD],
    pareto_exponents        = [1.1, 2.0, 5.0],   # only used for :PL
    mod_gammas              = [1.0, 5.0, 10.0],  # only used for :MOD
    pyramid_skewness = [1.0, 0.1, 0.01],
    abundance_distribution = [:Log, :Normal],
    B_term = [false],
    tspan=(0.0,50.0),
    t_perturb=25.0,
    max_calib=10,
    abundance_mean=1.0,
    atol=1.0,
    xi_threshold=0.7,
    number_of_combos=1000,
    calibrate_params_constraints = true,
    species_specific_perturbation = false,
)
    # build all combinations, but only pair pexs for :PL and mod_gamma for :MOD
    combos = Tuple[]
    for S in S_vals, conn in conn_vals, C_ratio in C_ratios,
        IS in IS_vals, d in d_vals, m in mortality_vals,
        eps in epsilon_vals, delta in delta_vals,
        scenario in degree_distribution_types,
        skew in pyramid_skewness,
        abund in abundance_distribution,
        B_term in B_term

        if scenario == :PL
            for pexs in pareto_exponents
                push!(combos, (S,conn,C_ratio,IS,d,m,eps,delta,scenario,pexs,0.0,skew,abund,B_term))
            end
        elseif scenario == :MOD
            for modg in mod_gammas
                push!(combos, (S,conn,C_ratio,IS,d,m,eps,delta,scenario,0.0,modg,skew,abund,B_term))
            end
        else  # :ER
            push!(combos, (S,conn,C_ratio,IS,d,m,eps,delta,scenario,0.0,0.0,skew,abund,B_term))
        end
    end

    T = nthreads()
    local_recs = [Vector{NamedTuple}() for _ in 1:T]

    @threads for idx in sample(1:length(combos), min(number_of_combos, length(combos)), replace=false)
        tid = threadid()
        S, conn, C_ratio, IS, d_val, m_val, eps_mean, delta, scenario, pexs, modg, skew, abund, B_term = combos[idx]

        # resource / consumer split
        C = clamp(round(Int, S*C_ratio), 1, S-1)
        R = S - C

        # pick correct callback
        cb = eval(Symbol("cb_no_trigger$(S)"))

        # 1) interaction matrix
        A = zeros(S,S)
        A = make_A(A, R, conn, scenario;
                    pareto_exponent=pexs,
                    mod_gamma=modg,
                    B_term=B_term)
        A = A .* IS

        # 2) equilibrium draws: pyramid skew via eps_mean*modg? no, keep independent
        if abund == :Log
            R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
        elseif abund == :Normal
            R_eq = abs.(rand(Normal(abundance_mean, abundance_mean), R))
            C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
        end
        fixed = vcat(R_eq, C_eq)

        # 3) rates
        m_cons      = fill(m_val, C)
        d_res       = fill(d_val, R)
        ε_mat       = clamp.(rand(Normal(eps_mean, eps_mean), S, S), 0, 1)

        # 4) calibrate
        xi_cons, r_res = calibrate_params(R_eq, C_eq,
            (R,C,m_cons,d_res,ε_mat,A);
            xi_threshold=xi_threshold,
            constraints=calibrate_params_constraints
        )
        tries = 1
        while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
            A .= 0
            A = make_A(A, R, conn, scenario;
                        pareto_exponent=pexs,
                        mod_gamma=modg,
                        B_term=B_term)
            A = A .* IS
            
            if abund == :Log
                R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
                C_eq = abs.(rand(LogNormal(log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
            elseif abund == :Normal
                R_eq = abs.(rand(Normal(abundance_mean, abundance_mean), R))
                C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
            end
            fixed = vcat(R_eq, C_eq)
            xi_cons, r_res = calibrate_params(R_eq, C_eq,
                (R,C,m_cons,d_res,ε_mat,A);
                xi_threshold=xi_threshold,
                constraints=calibrate_params_constraints
            )
            tries += 1
        end

        # if calibration fails
        if calibrate_params_constraints && any(isnan, xi_cons) || any(isnan, r_res)
            push!(local_recs[tid], (
                S=S, conn=conn, C_ratio=C_ratio,
                IS=IS, d=d_val, m=m_val, epsilon=eps_mean,
                delta=delta, scenario=scenario,
                R=R, C=C,
                feasible=false,
                constraints = calibrate_params_constraints,
                ssp = species_specific_perturbation,
                before_p=0.0, after_p=0.0,
                pexs=pexs, mod_gamma=modg,
                scorr=0.0, skew = skew,
                abundance_distribution = abund,
                B_term = B_term
            ))
            continue
        elseif !calibrate_params_constraints && any(isnan, r_res)
            push!(local_recs[tid], (
                S=S, conn=conn, C_ratio=C_ratio,
                IS=IS, d=d_val, m=m_val, epsilon=eps_mean,
                delta=delta, scenario=scenario,
                R=R, C=C,
                feasible=false,
                constraints = calibrate_params_constraints,
                ssp = species_specific_perturbation,
                before_p=0.0, after_p=0.0,
                pexs=pexs, mod_gamma=modg,
                scorr=0.0, skew = skew,
                abundance_distribution = abund,
                B_term = B_term
            ))
            continue
        end

        # 5) simulate unperturbed
        p_full = (R,C,m_cons,xi_cons,r_res,d_res,ε_mat,A)
        sol    = solve(ODEProblem(trophic_ode!, fixed, tspan, p_full),
                        Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
        B_eq   = sol.u[end]
        
        if sol.t[end] < t_perturb || any(isnan, sol.u[end]) || any(isinf, sol.u[end]) ||
            any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:S]) || any(sol.u[end] .< EXTINCTION_THRESHOLD)
            push!(local_recs[tid], (
                S=S, conn=conn, C_ratio=C_ratio,
                IS=IS, d=d_val, m=m_val, epsilon=eps_mean,
                delta=delta, scenario=scenario,
                R=R, C=C,
                feasible=false,
                constraints = calibrate_params_constraints,
                ssp = species_specific_perturbation,
                before_p=0.0, after_p=0.0,
                pexs=pexs, mod_gamma=modg,
                scorr=0.0, skew = skew,
                abundance_distribution = abund,
                B_term = B_term
            ))
            continue
        else
            
            before_p = mean(B_eq .> EXTINCTION_THRESHOLD)

            # 6) press perturbation
            rt, os, ire, _, B2, _ = simulate_press_perturbation(
                fixed, p_full, tspan, t_perturb, delta;
                solver=Tsit5(), plot=false,
                show_warnings=true,
                cb=cb,
                species_specific_perturbation=species_specific_perturbation
            )

            # 7) analytic vs sim scorr
            # 1) get the Jacobian at the “true” equilibrium B_post0
            D, Mstar   = compute_jacobian(B_eq, p_full)
            J_simp          = D * Mstar

            I_mat    = I(R+C)  # identity matrix (R+C)×(R+C)

            # 2) build the press‐vector properly: ∂f_C/∂ξ_i = –B_C* at equilibrium
            press_vec      = zeros(R+C)
            press_vec[R+1:R+C] .= 1.0

            # analytic per-unit sensitivity
            V            = -inv(I_mat .- A)    # (R+C)×(R+C)
            deltaB_ana  = V * press_vec          # length R+C

            # simulated per-unit response (length R+C)
            deltaB_sim = (B2 .- B_eq) ./ delta

            zero_idx = findall(iszero.(deltaB_sim))
                
            deltaB_sim[zero_idx] .= NaN
            deltaB_ana[zero_idx] .= NaN
            filter!(x -> !isnan(x), deltaB_sim)
            filter!(x -> !isnan(x), deltaB_ana)

            scorr = round(cor(deltaB_ana, deltaB_sim), digits=3)

            after_p = mean(B2 .> EXTINCTION_THRESHOLD)

            push!(local_recs[tid], (
                S=S, conn=conn, C_ratio=C_ratio,
                IS=IS, d=d_val, m=m_val, epsilon=eps_mean,
                delta=delta, scenario=scenario,
                R=R, C=C,
                feasible=true,
                constraints = calibrate_params_constraints,
                ssp = species_specific_perturbation,
                before_p=before_p, after_p=after_p,
                pexs=pexs, mod_gamma=modg,
                scorr=scorr, skew = skew,
                abundance_distribution = abund,
                B_term = B_term
            ))
        end
    end

    # gather and return
    return DataFrame(reduce(vcat, local_recs))
end
