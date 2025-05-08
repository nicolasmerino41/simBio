include("Ladder4.1_drago1.jl")

# -----------------------------------------------------------
# Main simulation pipeline
# -----------------------------------------------------------
function rpp(
    S_vals, conn_vals, delta_vals;
    max_RC_steps = 10,
    ladder_steps=1:12,
    Niter=10,
    tspan=(0.0,50.0),
    t_perturb=25.0,
    max_calib=10,
    plot_full = false,
    plot_simple = false,
    atol = 10.0,
    d_value = 1.0,
    IS_reduction = 1.0,
    homogenise = nothing,
    plot_steps = false,
    abundance_mean = 1.0,
    pred_mortality = 0.1,
    epsilon_mean = 0.2
)
    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()

    # suffix helper
    suffix(step) = step==1 ? "Full" : "S$(step-1)"
    
    S_conn_pairs = [(S, conn) for S in S_vals for conn in conn_vals]
    # Generate all combinations of S and conn
    for (S, conn) in S_conn_pairs
        for delta in delta_vals
            # Generate all valid (R, C) combinations such that R + C = S_val
            for R in pick_steps(5, S-1, max_RC_steps)
                C = S - R
                total = S
                Threads.@threads for iter in 1:Niter
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
                            if !isnothing(homogenise)
                                A_adj = zeros(total,total)
                                for i in (R+1):total, j in 1:total
                                    if i !=j && rand() < conn && iszero(A_adj[i,j])
                                        A_adj[i,j]= homogenise
                                        A_adj[j,i]= -homogenise
                                    end
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
                        rt, overshoot, ire, _, B2 = simulate_press_perturbation(
                            fixed, p_full, tspan, t_perturb, delta;
                            solver=Tsit5(), plot=plot_full,
                            show_warnings=true,
                            full_or_simple = true
                        )

                        pred_resp = sum(V_analytical, dims=2)
                        obs_resp  = (B2 .- fixed) ./ delta

                        # flatten both to plain vectors and correlate
                        sens_corr_full = cor(vec(pred_resp), vec(obs_resp))

                        metrics["Full"] = (
                            persistence=pers,
                            return_time=mean(filter(!isnan, rt)),
                            overshoot=mean(filter(!isnan, overshoot)),
                            ire=mean(filter(!isnan, ire)),
                            before_p=pers,
                            after_p=mean(B2 .> EXTINCTION_THRESHOLD),
                            resilience=resi,
                            reactivity=reac
                        )
                        # E) Attach to your Full‐step metrics
                        metrics["Full"] = merge(
                            metrics["Full"],
                            (sensitivity_corr = sens_corr_full,)
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
                                persistence=NaN,
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
                            pers2 = mean(Beq2 .> EXTINCTION_THRESHOLD)
                            res2  = compute_resilience(Beq2,p_simp)
                            rea2  = compute_reactivity(Beq2,p_simp)
                            # get the full per-species press matrix just like in the full model
                            rt2, overshoot, ire, _, B2_simp = simulate_press_perturbation(
                                    fixed, p_simp, tspan, t_perturb, delta;
                                    solver=Tsit5(), plot=plot_simple,
                                    show_warnings=true,
                                    full_or_simple=true
                                )
            
                            # in each step, after computing V_ana_simp and B2_simp
                            pred_resp_simp = sum(V_ana_simp, dims=2)
                            obs_resp_simp  = (B2_simp .- fixed) ./ delta
                            corr_simp      = cor(vec(pred_resp_simp), vec(obs_resp_simp))

                            metrics[suf] = (
                                persistence=pers2,
                                return_time=mean(filter(!isnan, rt2)),
                                overshoot=mean(filter(!isnan, overshoot)),
                                ire=mean(filter(!isnan, ire)),
                                before_p=pers2,
                                after_p=mean(B2_simp .> EXTINCTION_THRESHOLD),
                                resilience=res2,
                                reactivity=rea2
                            )
                            # update metrics for this step
                            metrics[suf] = merge(
                                metrics[suf],
                                (sensitivity_corr = corr_simp,)
                            )
                        end

                        # flatten into one NamedTuple
                        base = (
                            true_stable = !any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=1e-3) for i in 1:total]),
                            species_count=total,
                            resource_count=R,
                            consumer_count=C,
                            connectance=conn,
                            perturb_delta=delta,
                            iteration=iter,
                            RY_vector=RY_vector
                        )
                        rec = base
                        for step in ladder_steps
                            suf = suffix(step)
                            m = metrics[suf]
                            rec = merge(
                                rec,
                                (; Symbol("persistence_$suf")      => m.persistence,
                                Symbol("return_time_$suf")      => m.return_time,
                                Symbol("overshoot_$suf")        => m.overshoot,
                                Symbol("ire_$suf")              => m.ire,
                                Symbol("before_persistence_$suf")=> m.before_p,
                                Symbol("after_persistence_$suf") => m.after_p,
                                Symbol("resilience_$suf")        => m.resilience,
                                Symbol("reactivity_$suf")        => m.reactivity),
                                (; Symbol("sensitivity_corr_$suf") => m.sensitivity_corr)
                            )
                        end

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
    println("$(nrow(results)) / $(Niter*length(connectance_list)*length(delta_list)), simulations succesful")
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
    S_val = [75]
    connectance_list = 0.1:0.1:1.0
    delta_list = 0.1:0.1:0.1

    new_results9 = rpp(
        S_val, connectance_list, delta_list;
        ladder_steps=1:16,
        Niter=50,
        tspan=(0.0, 500.0),
        t_perturb=250.0,
        max_calib=10,
        plot_full=false,
        plot_simple=false,
        atol = 10.0,
        d_value = 2.0,
        IS_reduction = 1.0,
        homogenise = nothing,
        plot_steps = false,
        abundance_mean = 10.0,
        pred_mortality = 0.05,
        epsilon_mean = 0.1
    )
end
serialize("Results/new_results9.jls", new_results9)
CSV.write("Results/new_results9.csv", new_results9)
