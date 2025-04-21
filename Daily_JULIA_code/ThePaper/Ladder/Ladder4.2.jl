include("Ladder4.1_drago1.jl")

# -----------------------------------------------------------
# Main simulation pipeline with compound‐error and degree‐CV
# -----------------------------------------------------------
function rppp(
    S_vals, conn_vals, delta_vals;
    max_RC_steps = 10,
    ladder_steps = 1:12,
    Niter        = 10,
    tspan        = (0.0,50.0),
    t_perturb    = 25.0,
    max_calib    = 10,
    plot_full    = false,
    plot_simple  = false,
    plot_steps   = false,
    atol         = 10.0,
    d_value      = 1.0,
    IS_reduction = 1.0,
    abundance_mean = 1.0,
    pred_mortality = 0.1,
    epsilon_mean   = 0.2
)
    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()

    suffix(step) = step == 1 ? "Full" : "S$(step-1)"
    S_conn_pairs = [(S,conn) for S in S_vals for conn in conn_vals]

    Threads.@threads for (S, conn) in S_conn_pairs
        for delta in delta_vals
            for R in pick_steps(5, S-1, max_RC_steps)
                C = S - R
                total = S

                for iter in 1:Niter
                    println("Running iteration $iter for R=$R, C=$C, conn=$conn, delta=$delta")
                    try
                        A_adj = zeros(total,total)
                        for i in (R+1):total, j in 1:total
                            if i !=j && rand() < conn && iszero(A_adj[i,j])
                                A_adj[i,j]= rand() * IS_reduction
                                A_adj[j,i]= -rand() * IS_reduction
                            end
                        end
                        
                        # 1) target equilibrium
                        fixed = vcat(abs.(rand(LogNormal(abundance_mean, abundance_mean*2), R)), abs.(rand(LogNormal(abundance_mean*0.1, abundance_mean*0.2), C)))
                        R_eq, C_eq = fixed[1:R], fixed[R+1:end]

                        # 2) fixed params
                        m_cons = fill(pred_mortality, C)
                        d_res  = fill(d_value, R)

                        # 3) full ε
                        # ε_full = fill(epsilon_mean, total, total)
                        ε_full = clamp.(rand(LogNormal(epsilon_mean, epsilon_mean), total, total), 0, 1)

                        # 4) calibrate ONCE with full A & ε
                        p_cal = (R, C, m_cons, d_res, ε_full, A_adj)
                        xi_cons, r_res = calibrate_params(R_eq, C_eq, p_cal)

                        # maybe retry if NaN
                        tries=1
                        while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
                            # println("i tried")
                            # build binary adjacency once
                            A_adj = zeros(total,total)
                            for i in (R+1):total, j in 1:total
                                if i !=j && rand() < conn && iszero(A_adj[i,j])
                                    A_adj[i,j]= rand() * IS_reduction
                                    A_adj[j,i]= -rand() * IS_reduction # -A_adj[i,j]
                                end
                            end

                            fixed = vcat(abs.(rand(LogNormal(abundance_mean, abundance_mean*2), R)), abs.(rand(LogNormal(abundance_mean*0.1, abundance_mean*0.2), C)))
                            p_cal = (R, C, m_cons, d_res, ε_full, A_adj)
                            
                            R_eq, C_eq = fixed[1:R], fixed[R+1:end]
                            xi_cons, r_res = calibrate_params(R_eq, C_eq, p_cal)

                            tries+=1
                        end

                        if any(isnan,xi_cons)||any(isnan,r_res)
                            # @warn "Calibration failed; skipping iteration"
                            continue
                        end
                        RY_vector = C_eq./xi_cons
                        # dict to hold all step‐metrics
                        metrics = Dict{String,NamedTuple}()

                        # run FULL first (step=1)
                        p_full = (R, C, m_cons, xi_cons, r_res, d_res, ε_full, A_adj)
                        
                        g = SimpleGraph(A_adj .!= 0)
                        degs = degree(g)
                        degree_cv = std(degs) / mean(degs)
                        ###################################################################
                        ################### EXPERIMENTAL PART #############################
                        #  A) Build Jacobian & analytical V for the full model
                        J_full       = build_jacobian(fixed, p_full)
                        V_analytical = compute_analytical_V(J_full, R, C, m_cons, xi_cons)

                        prob = ODEProblem(trophic_ode!, fixed, tspan, p_full)
                        sol = solve(prob, reltol=1e-8,abstol=1e-8)
                        # if any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=1e-3) for i in 1:total])
                        #     println("The final abundances do not match the target equilibrium at conn = $conn")
                        # else
                        #     println("The final abundances match the target equilibrium at conn = $conn")
                        # end
                        if sol.t[end] < t_perturb || any(isnan, sol.u[end]) || any(isinf, sol.u[end]) || any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:total])
                            # @warn "Error: solution did not finish properly"
                            continue
                        end
                        # begin
                        #     fig = Figure(; size=(800, 600), title="Trophic Dynamics")
                        #     ax = Axis(fig[1, 1], xlabel="Time", ylabel="Population Size", title="Trophic Dynamics")
                        #     for i in 1:(R)
                        #         MK.lines!(ax, sol.t, sol[i, :], label="Resource $i", color=:blue, linewidth=1)
                        #     end
                        #     for i in 1:(C)
                        #         MK.lines!(ax, sol.t, sol[R+i, :], label="Consumer $i", color=:red, linewidth=1)
                        #     end
                        #     display(fig)
                        # end

                        Beq = sol.u[end]
                        pers = mean(Beq .> EXTINCTION_THRESHOLD)
                        resi = compute_resilience(Beq, p_full)
                        reac = compute_reactivity(Beq, p_full)
                        rt_full, os_full, ire_full, _, B2 = simulate_press_perturbation(
                            fixed, p_full, tspan, t_perturb, delta;
                            solver=Tsit5(), plot=plot_full,
                            show_warnings=true,
                            full_or_simple = true
                        )

                        pred_resp = sum(V_analytical, dims=2)
                        obs_resp  = (B2 .- fixed) ./ delta

                        # flatten both to plain vectors and correlate
                        sens_corr_full = cor(vec(pred_resp), vec(obs_resp))

                        # derive predicted return time from resilience
                        t_pred   = 1/resi
                        t_obs    = mean(filter(!isnan, rt_full))
                        err_T    = abs(t_obs - t_pred) / ((t_obs + t_pred)/2)
                        err_OS   = mean(filter(!isnan, os_full))
                        err_IRE  = mean(filter(!isnan, ire_full))
                        compound_error_full = 0.0
                        
                        ### 6) Store full‐model metrics ###
                        metrics = Dict{String,NamedTuple}()
                        metrics["Full"] = (
                            compound_error  = compound_error_full,
                            return_time     = mean(filter(!isnan, rt_full)),
                            overshoot       = mean(filter(!isnan, os_full)),
                            ire             = mean(filter(!isnan, ire_full)),
                            before_p        = mean(sol.u[end] .> EXTINCTION_THRESHOLD),
                            after_p         = mean(B2 .> EXTINCTION_THRESHOLD),
                            resilience      = resi,
                            reactivity      = reac,
                            sensitivity_corr       = sens_corr_full,
                        )

                        # now each simplified step
                        for step in ladder_steps[2:end]
                            suf = suffix(step)
                            A_used, ε_used = transform_for_ladder_step(step, A_adj, ε_full)
                            p_simp = (R,C,m_cons,xi_cons,r_res,d_res,ε_used,A_used)
                            J_simp         = build_jacobian(fixed, p_simp)
                            V_ana_simp     = compute_analytical_V(J_simp, R, C, m_cons, xi_cons)
                            
                            sol2 = solve(
                                ODEProblem(trophic_ode!, fixed, tspan, p_simp), Tsit5(),
                                reltol=1e-8,abstol=1e-8
                            )

                            if sol2.t[end] < t_perturb || any(isnan, sol2.u[end]) || any(isinf, sol2.u[end])
                                metrics[suf] = (
                                compound_error=NaN,
                                return_time=NaN,
                                overshoot=NaN,
                                ire=NaN,
                                before_p=NaN,
                                after_p=NaN,
                                resilience=NaN,
                                reactivity=NaN,
                                sensitivity_corr=NaN
                            )
                                continue
                            end
                            if plot_steps
                                fig = Figure(; size=(800, 600), title="Trophic Dynamics Step $suf")
                                ax = Axis(fig[1, 1], xlabel="Time", ylabel="Population Size", title="Trophic Dynamics Step $suf")
                                for i in 1:(R)
                                    MK.lines!(ax, sol2.t, sol2[i, :], label="Resource $i", color=:blue, linewidth=1)
                                end
                                for i in 1:(C)
                                    MK.lines!(ax, sol2.t, sol2[R+i, :], label="Consumer $i", color=:red, linewidth=1)
                                end
                                display(fig)
                            end

                            Beq2 = sol2.u[end]
                            before_p = mean(sol2.u[end] .> EXTINCTION_THRESHOLD)
                            res2  = compute_resilience(Beq2,p_simp)
                            rea2  = compute_reactivity(Beq2,p_simp)
                            # get the full per-species press matrix just like in the full model
                            rt2, os2, ire2, _, B2_simp = simulate_press_perturbation(
                                    fixed, p_simp, tspan, t_perturb, delta;
                                    solver=Tsit5(), plot=plot_simple,
                                    show_warnings=true,
                                    full_or_simple=true
                                )
            
                            # in each step, after computing V_ana_simp and B2_simp
                            pred_resp_simp = sum(V_ana_simp, dims=2)
                            obs_resp_simp  = (B2_simp .- fixed) ./ delta
                            corr_simp      = cor(vec(pred_resp_simp), vec(obs_resp_simp))

                            # compound error, step
                            t_pred2 = mean(filter(!isnan, rt2))
                            t_obs2  = mean(filter(!isnan, rt2))
                            err_T2   = abs(t_obs2 - t_pred2) / ((t_obs2 + t_pred2)/2)
                            os_pred2 = mean(filter(!isnan, os_full))
                            os_obs = mean(filter(!isnan, os2))
                            err_OS2  = abs(os_obs - os_pred2) / ((os_obs + os_pred2)/2)
                            ire_pred2 = mean(filter(!isnan, ire_full))
                            ire_obs = mean(filter(!isnan, ire2))
                            err_IRE2  = abs(ire_obs - ire_pred2) / ((ire_obs + ire_pred2)/2)
                            comp_err2 = (err_T2 + err_OS2 + err_IRE2)/3

                            metrics[suf] = (
                                compound_error = comp_err2,
                                return_time    = mean(filter(!isnan, rt2)),
                                overshoot      = mean(filter(!isnan, os2)),
                                ire            = mean(filter(!isnan, ire2)),
                                before_p       = before_p,
                                after_p        = mean(B2_simp .> EXTINCTION_THRESHOLD),
                                resilience     = res2,
                                reactivity     = rea2,
                                sensitivity_corr      = corr_simp,
                            )
                        end

                        ### 8) Flatten and save ###
                        base = (
                            species_count   = total,
                            resource_count  = R,
                            consumer_count  = C,
                            connectance     = conn,
                            perturb_delta   = delta,
                            iteration       = iter,
                            RY_vector       = RY_vector,
                            degree_cv       = degree_cv
                        )
                        rec = base
                        for step in ladder_steps
                            suf = suffix(step)
                            m = metrics[suf]
                            rec = merge(rec,
                                (; Symbol("compound_error_$suf")  => m.compound_error,
                                   Symbol("return_time_$suf")     => m.return_time,
                                   Symbol("overshoot_$suf")       => m.overshoot,
                                   Symbol("ire_$suf")             => m.ire,
                                   Symbol("before_persistence_$suf") => m.before_p,
                                   Symbol("after_persistence_$suf")  => m.after_p,
                                   Symbol("resilience_$suf")       => m.resilience,
                                   Symbol("reactivity_$suf")       => m.reactivity),
                                   (; Symbol("sensitivity_corr_$suf") => m.sensitivity_corr)
                            )
                        end

                        lock(results_lock) do
                            push!(results, rec)
                        end

                    catch e
                        @warn "Iteration error" S conn iter e
                    end
                end
            end
        end
    end

    return DataFrame(results)
end

# helper: pick up to max_steps integers between minR and maxR
function pick_steps(minR::Int, maxR::Int, max_steps::Int)
    allR = collect(minR:maxR)
    if length(allR) <= max_steps
        return allR
    else
        # LinRange gives you max_steps values from minR to maxR, but as Float64
        # round to Int and unique to make sure you still span the ends
        picked = unique(round.(Int, LinRange(minR, maxR, max_steps)))
        # in pathological cases rounding could collapse two endpoints; 
        # force inclusion of both ends:
        push!(picked, minR, maxR)
        return sort(unique(picked))
    end
end

begin
    # Specify parameter ranges:
    S_val = [50]
    connectance_list = 0.1:0.2:1.0
    delta_list = 0.1:0.1:0.1

    new_results10 = rppp(
        S_val, connectance_list, delta_list;
        max_RC_steps = 5,
        ladder_steps=1:16,
        Niter=1,
        tspan=(0.0, 500.0),
        t_perturb=250.0,
        max_calib=5,
        plot_full=false,
        plot_simple=false,
        atol = 10.0,
        d_value = 1.0,
        IS_reduction = 1.0,
        abundance_mean = 1.0,
        pred_mortality = 0.1,
        epsilon_mean = 0.1
    )
end

serialize("Results/new_results10.jls", new_results10)
CSV.write("Results/new_results10.csv", new_results10)