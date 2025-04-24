include("Ladder1.jl")

################################################################################
################################################################################
######################### RUN PIPELINE #########################################
################################################################################
################################################################################
function rp(
    R_vals, C_vals, conn_vals, delta_vals;
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

    Threads.@threads for conn in conn_vals
        for R in R_vals, C in C_vals, delta in delta_vals
            total = R+C

            for iter in 1:Niter
                # println("run")
                try
                    A_adj = zeros(total,total)
                    for i in (R+1):total, j in 1:total
                        if i !=j && rand() < conn && iszero(A_adj[i,j])
                            A_adj[i,j]= rand() * IS_reduction
                            A_adj[j,i]= -rand() * IS_reduction
                        end
                    end
                    
                    # 1) target equilibrium
                    fixed = vcat(abs.(rand(Normal(abundance_mean, abundance_mean*2), R)), abs.(rand(Normal(abundance_mean*0.1, abundance_mean*0.2), C)))
                    R_eq, C_eq = fixed[1:R], fixed[R+1:end]

                    # 2) fixed params
                    m_cons = fill(pred_mortality, C)
                    d_res  = fill(d_value, R)

                    # 3) full ε
                    # ε_full = fill(epsilon_mean, total, total)
                    ε_full = clamp.(rand(Normal(epsilon_mean, epsilon_mean), total, total), 0, 1)

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

                        fixed = vcat(abs.(rand(Normal(abundance_mean, abundance_mean*2), R)), abs.(rand(Normal(abundance_mean*0.1, abundance_mean*0.2), C)))
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
                    rt,_,B2 = simulate_press_perturbation(
                        fixed, p_full, tspan, t_perturb, delta;
                        solver=Tsit5(), plot=plot_full,
                        show_warnings=true,
                        full_or_simple = true
                    )

                    # sum across columns to get the predicted deltaB for a uniform δξ
                    pred_resp = sum(V_analytical, dims=2)           # S×1 matrix
                    obs_resp  = (B2 .- fixed) ./ delta              # S‑vector

                    # flatten both to plain vectors and correlate
                    sens_corr_full = cor(vec(pred_resp), vec(obs_resp))

                    metrics["Full"] = (
                        persistence=pers,
                        return_time=mean(filter(!isnan, rt)),
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
                        rt2, _, B2_simp = simulate_press_perturbation(
                                fixed, p_simp, tspan, t_perturb, delta;
                                solver=Tsit5(), plot=plot_simple,
                                show_warnings=true,
                                full_or_simple=true
                            )
                        # now B2_simp is total×total, matching V_ana_simp
                        V_sim_simp = (B2_simp .- fixed) ./ delta
        
                        # in each step, after computing V_ana_simp and B2_simp
                        pred_resp_simp = sum(V_ana_simp, dims=2)
                        obs_resp_simp  = (B2_simp .- fixed) ./ delta
                        corr_simp      = cor(vec(pred_resp_simp), vec(obs_resp_simp))


                        metrics[suf] = (
                            persistence=pers2,
                            return_time=mean(filter(!isnan, rt2)),
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
    println("$(nrow(results)) / $(Niter*length(connectance_list)*length(delta_list)), simulations succesful")
    return DataFrame(results)
end

##########################################################################
##########################################################################
################### Ladder‐of‐simplification transforms ##################
##########################################################################
##########################################################################
function transform_for_ladder_step(step, A_adj, ε_full)
    total = size(A_adj,1)

    if step == 1 # Full
        return A_adj, ε_full
    elseif step == 2 # S2 epsilon mean per row
        ε2 = similar(ε_full)
        for i in 1:total
            ε2[i,:] .= mean(ε_full[i,:])
        end
        return A_adj, ε2
    elseif step == 3 # S3 epsilon global mean
        return A_adj, fill(mean(ε_full), total, total)
    elseif step == 4 # S4 epsilon re-randomised
        ε4 = similar(ε_full)
        for i in 1:total, j in 1:total
            ε4[i, j] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
        end
        return A_adj, ε4
    elseif 5 ≤ step ≤ 8 # A is averaged separately for positive and negative values
        
        pos, neg = mean(A_adj[A_adj.>0]), mean(A_adj[A_adj.<0])
        
        A5 = map(x-> x>0 ? pos : x<0 ? neg : 0, A_adj)
        if step==5 # ϵ is full
            return A5, ε_full
        elseif step==6 # ϵ is mean per row
            ε6 = similar(ε_full)
            for i in 1:total 
                ε6[i,:] .= mean(ε_full[i,:])
            end
            return A5, ε6
        elseif step==7 # ϵ is global mean
            ε7 = fill(mean(ε_full), total, total)
            return A5, ε7
        elseif step==8 # ϵ is re-randomised
            ε8 = similar(ε_full)
            for i in 1:total, j in 1:total
                ε8[i, j] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
            end
            return A5, ε8
        end
    elseif 9 ≤ step ≤ 12 # A is averaged across all non-zero values
        
        m = mean(abs.(A_adj[A_adj .!= 0]))
        A6 = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)
        
        if step==9 # ϵ is full
            return A6, ε_full
        elseif step==10 # ϵ is mean per row
            ε10 = similar(ε_full)
            for i in 1:total
                ε10[i,:] .= mean(ε_full[i,:])
            end
            return A6, ε10
        elseif step==11 # ϵ is global mean
            return A6, fill(mean(ε_full), total, total)
        elseif step==12 # ϵ is re-randomised
            ε12 = similar(ε_full)
            for i in 1:total*total
                ε12[i] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
            end
            return A6, ε12
        end
    elseif 13 ≤ step ≤ 16 # A is fully randomised with same sparsity pattern as A_adj
        # Fully randomized A matrix using Normal(mean, std) of abs non-zero A entries
        A_vals = abs.(A_adj[A_adj .!= 0])
        A_mean, A_std = mean(A_vals), std(A_vals)
    
        A_rand = similar(A_adj)
        for i in 1:total, j in 1:total
            if A_adj[i, j] != 0
                A_rand[i, j] = rand(Normal(A_mean, A_std)) * sign(A_adj[i, j])
            else
                A_rand[i, j] = 0.0
            end
        end
    
        if step == 13 # ϵ is full
            return A_rand, ε_full
        elseif step == 14 # ϵ is mean per row
            ε14 = similar(ε_full)
            for i in 1:total
                ε14[i, :] .= mean(ε_full[i, :])
            end
            return A_rand, ε14
        elseif step == 15 # ϵ is global mean
            return A_rand, fill(mean(ε_full), total, total)
        elseif step == 16 # ϵ is re-randomised
            ε16 = similar(ε_full)
            for i in 1:total, j in 1:total
                ε16[i, j] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
            end
            return A_rand, ε16
        end    
    else
        error("Unknown ladder step $step")
    end
end

if false
    # Specify parameter ranges:
    R_vals = [10, 20, 30, 40]            # e.g., try 4 and 6 resources
    C_vals = [5, 10, 15, 20]            # e.g., try 5 and 8 consumers
    connectance_list = 0.1:0.1:1.0
    delta_list = 0.1:0.1:0.1

    new_results3 = rp(
        R_vals, C_vals, connectance_list, delta_list;
        ladder_steps=1:16,
        Niter=100,
        tspan=(0.0, 500.0),
        t_perturb=250.0,
        max_calib=10,
        plot_full=false,
        plot_simple=false,
        atol = 10.0,
        d_value = 1.0,
        IS_reduction = 1.0,
        homogenise = nothing,
        plot_steps = false,
        abundance_mean = 1.0,
        pred_mortality = 0.1,
        epsilon_mean = 0.2
    )
end
# serialize("Daily_JULIA_Code/ThePaper/Ladder/new_results3.jls", new_results3)
new_results3 = deserialize("Daily_JULIA_Code/ThePaper/Ladder/new_results3.jls")

println("First 10 rows of simulation results:")
println(first(df_results, 10))