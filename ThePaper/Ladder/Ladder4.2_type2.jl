# -----------------------------------------------------------
# Main simulation pipeline with Holling‐Type II response
# -----------------------------------------------------------
function rp2(
    S_vals, conn_vals, delta_vals;
    max_RC_steps    = 10,
    ladder_steps    = 1:12,
    Niter           = 10,
    tspan           = (0.0,50.0),
    t_perturb       = 25.0,
    max_calib       = 10,
    plot_full       = false,
    plot_simple     = false,
    plot_steps      = false,
    atol            = 10.0,
    d_value         = 1.0,
    IS_reduction    = 1.0,
    abundance_mean  = 1.0,
    pred_mortality  = 0.1,
    epsilon_mean    = 0.2,
    handling_time   = 0.1
)
    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()

    suffix(step) = step == 1 ? "Full" : "S$(step-1)"
    S_conn_pairs = [(S,conn) for S in S_vals for conn in conn_vals]

    Threads.@threads for (S, conn) in S_conn_pairs
        for delta in delta_vals
            for R in pick_steps(1, S-5, max_RC_steps)
                C = S - R
                total = S
                for iter in 1:Niter
                    try
                        # 1) build random adjacency
                        A = zeros(total,total)
                        for i in (R+1):total, j in 1:total
                            if i!=j && rand()<conn
                                A[i,j]= rand()*IS_reduction
                                A[j,i]= -rand()*IS_reduction
                            end
                        end

                        # 2) draw equilibrium biomasses
                        R_eq = abs.(rand(LogNormal(log(abundance_mean),abundance_mean*2), R))
                        C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1),abundance_mean*0.2), C))
                        fixed = vcat(R_eq, C_eq)

                        # 3) parameters
                        m_cons = fill(pred_mortality, C)
                        d_res  = fill(d_value, R)
                        ε_full = clamp.(rand(LogNormal(log(epsilon_mean),epsilon_mean), total, total), 0,1)
                        h      = fill(handling_time, total, total)

                        # 4) calibrate once
                        p_cal = (R, C, m_cons, d_res, ε_full, A, h)
                        xi_cons, r_res = calibrate_two(R_eq, C_eq, p_cal)
                        tries = 1
                        while (any(isnan, xi_cons) || any(isnan, r_res)) && tries<max_calib
                            # rebuild A, then recalibrate
                            A .= 0
                            for i in (R+1):total, j in 1:total
                                if i!=j && rand()<conn
                                    A[i,j]= rand()*IS_reduction
                                    A[j,i]= -rand()*IS_reduction
                                end
                            end
                            xi_cons, r_res = calibrate_two(R_eq, C_eq, (R,C,m_cons,d_res,ε_full,A,h))
                            tries+=1
                        end
                        if any(isnan, xi_cons) || any(isnan, r_res)
                            continue
                        end

                        # 5) structural metric
                        g = SimpleGraph(A .!= 0)
                        degs = degree(g)
                        degree_cv = std(degs)/mean(degs)

                        # 6) full model solve & metrics
                        p_full = (R, C, m_cons, xi_cons, r_res, d_res, ε_full, A, h)
                        sol = solve(ODEProblem(trophic_two!, fixed, tspan, p_full),
                                    Tsit5(); reltol=1e-8, abstol=1e-8)
                        if sol.t[end]<t_perturb || any(!isfinite, sol.u[end])
                            continue
                        end

                        Beq = sol.u[end]
                        rt_full, os_full, ire_full, _, B2 = simulate_press_perturbation(
                            fixed, p_full, tspan, t_perturb, delta;
                            solver=Tsit5(), plot=plot_full,
                            show_warnings=true, full_or_simple=true)

                        # analytical sensitivity
                        J_full       = build_jacobian(fixed, p_full)
                        V_analytical = compute_analytical_V(J_full, R, C, m_cons, xi_cons)
                        sens_corr_full = cor(vec(sum(V_analytical,dims=2)),
                                            vec((B2 .- fixed)./delta))

                        # derive compound error
                        resi = compute_resilience(Beq, p_full)
                        t_pred = 1/resi
                        t_obs  = mean(filter(!isnan, rt_full))
                        err_T  = abs(t_obs-t_pred)/((t_obs+t_pred)/2)
                        err_OS = mean(filter(!isnan, os_full))
                        err_IRE= mean(filter(!isnan, ire_full))
                        comp_err_full = (err_T + err_OS + err_IRE)/3

                        metrics = Dict{String,NamedTuple}()
                        metrics["Full"] = (
                            compound_error  = comp_err_full,
                            return_time     = t_obs,
                            overshoot       = err_OS,
                            ire             = err_IRE,
                            before_p        = mean(Beq .> EXTINCTION_THRESHOLD),
                            after_p         = mean(B2 .> EXTINCTION_THRESHOLD),
                            resilience      = resi,
                            reactivity      = compute_reactivity(Beq,p_full),
                            sensitivity_corr= sens_corr_full
                        )

                        # 7) simplified steps
                        for step in ladder_steps[2:end]
                            suf = suffix(step)
                            A_s, ε_s = transform_for_ladder_step(step, A, ε_full)
                            p_s = (R, C, m_cons, xi_cons, r_res, d_res, ε_s, A_s, h)
                            sol2 = solve(ODEProblem(trophic_two!, fixed, tspan, p_s),
                                        Tsit5(); reltol=1e-8, abstol=1e-8)
                            if sol2.t[end]<t_perturb || any(!isfinite, sol2.u[end])
                                metrics[suf] = NamedTuple()
                                continue
                            end
                            Beq2 = sol2.u[end]
                            rt2, os2, ire2, _, B2_s = simulate_press_perturbation(
                            fixed, p_s, tspan, t_perturb, delta;
                            solver=Tsit5(), plot=plot_simple,
                            show_warnings=true, full_or_simple=false)
                            J_s = build_jacobian(fixed, p_s)
                            V_s = compute_analytical_V(J_s, R, C, m_cons, xi_cons)

                            corr_s  = cor(vec(sum(V_s,dims=2)), vec((B2_s.-fixed)./delta))
                            res2    = compute_resilience(Beq2, p_s)
                            t_pred2 = 1/res2
                            t_obs2  = mean(filter(!isnan, rt2))
                            err_T2  = abs(t_obs2-t_pred2)/((t_obs2+t_pred2)/2)
                            err_OS2 = mean(filter(!isnan, os2))
                            err_IRE2= mean(filter(!isnan, ire2))
                            comp_err2 = (err_T2+err_OS2+err_IRE2)/3

                            metrics[suf] = (
                            compound_error   = comp_err2,
                            return_time      = t_obs2,
                            overshoot        = err_OS2,
                            ire              = err_IRE2,
                            before_p         = mean(Beq2 .> EXTINCTION_THRESHOLD),
                            after_p          = mean(B2_s .> EXTINCTION_THRESHOLD),
                            resilience       = res2,
                            reactivity       = compute_reactivity(Beq2,p_s),
                            sensitivity_corr = corr_s
                            )
                        end

                        # 8) flatten
                        base = (
                            species_count  = total,
                            resource_count = R,
                            consumer_count = C,
                            connectance    = conn,
                            perturb_delta  = delta,
                            iteration      = iter,
                            degree_cv      = degree_cv
                        )
                        rec = base
                        for step in ladder_steps
                            suf = suffix(step)
                            m   = metrics[suf]
                            rec = merge(rec,
                            (; Symbol("compound_error_$suf")   => m.compound_error,
                                Symbol("return_time_$suf")      => m.return_time,
                                Symbol("overshoot_$suf")        => m.overshoot,
                                Symbol("ire_$suf")              => m.ire,
                                Symbol("before_persistence_$suf")=> m.before_p,
                                Symbol("after_persistence_$suf") => m.after_p,
                                Symbol("resilience_$suf")        => m.resilience,
                                Symbol("reactivity_$suf")        => m.reactivity,
                                Symbol("sensitivity_corr_$suf")  => m.sensitivity_corr)
                            )
                        end

                        lock(results_lock) do
                            push!(results, rec)
                        end

                    catch e
                        @warn "Iteration error" S conn delta R iter e
                    end
                end
            end
        end
    end

    return DataFrame(results)
end

begin
    # Specify parameter ranges:
    S_vals = [30]
    conn_vals = 0.5:0.1:1.0
    delta_vals = 0.1:0.1:0.1

    new_results13 = rp2(
        S_vals, conn_vals, delta_vals;
        max_RC_steps    = 5,
        ladder_steps    = 1:12,
        Niter           = 10,
        tspan           = (0.0,500.0),
        t_perturb       = 250.0,
        max_calib       = 1,
        plot_full       = false,
        plot_simple     = false,
        plot_steps      = false,
        atol            = 10.0,
        d_value         = 1.0,
        IS_reduction    = 1.0,
        abundance_mean  = 1.0,
        pred_mortality  = 0.1,
        epsilon_mean    = 0.2,
        handling_time   = 0.1
    )
end
