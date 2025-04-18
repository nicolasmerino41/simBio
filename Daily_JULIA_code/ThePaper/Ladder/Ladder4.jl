"""
    cp(R_eq, C_eq, p_calib)

Given:
  - R_eq::Vector{Float64}: target resource biomasses (length R)
  - C_eq::Vector{Float64}: target consumer biomasses (length C)
  - p_calib = (R, C, m_cons, d_res, epsilon, A_used)

where
  * R, C        are numbers of resources/consumers,
  * m_cons      is a length‑C vector of mortalities,
  * d_res       is a length‑R vector of resource self‑damping,
  * epsilon     is a (R+C)x(R+C) efficiency matrix,
  * A_used      is the already‑weighted interaction matrix (positive for “i eats j”, negative for “j eats i”).

Returns `(new_xi_cons, new_r_res)` of lengths C and R respectively.
"""
function cp(R_eq::Vector{Float64},
            C_eq::Vector{Float64},
            p_calib)
    R, C, m_cons, d_res, epsilon, A_used = p_calib
    total = R + C
    B_eq = vcat(R_eq, C_eq)

    new_xi_cons = zeros(C)
    new_r_res   = zeros(R)

    # 1) Consumer thresholds
    for k in 1:C
        i = R + k

        # Gains from prey (only A_used[i,j]>0)
        prey_sum = sum(
            epsilon[i,j] * A_used[i,j] * B_eq[j]
            for j in 1:total if A_used[i,j] > 0.0;
            init = 0.0
        )

        # Losses to predators (only A_used[i,j]<0)
        pred_sum = sum(
            A_used[i,j] * B_eq[j]
            for j in 1:total if A_used[i,j] < 0.0;
            init = 0.0
        )

        ξ = -C_eq[k] + prey_sum - pred_sum
        new_xi_cons[k] = ξ > 0 ? ξ : NaN

        if prey_sum - pred_sum > new_xi_cons[k]
            # println("Hey, worked and RY is $(C_eq[k]/ξ) because C is $(C_eq[k]) and ξ is $(ξ)")
            if C_eq[k]/ξ > 0.3
                new_xi_cons[k] = NaN
            end
        else
            # println("Nope and RY is $(C_eq[k]/ξ) because C is $(C_eq[k]) and ξ is $(ξ)")
            new_xi_cons[k] = NaN
        end
    end

    # 2) Resource growth rates
    for i in 1:R
        # Consumer rows are R+1:total; losses where those entries <0
        pred_sum = sum(
            A_used[i,j] * B_eq[j]
            for j in 1:R+C if A_used[i,j] < 0.0;
            init = 0.0
        )
        r = d_res[i] * (R_eq[i] + pred_sum)
        new_r_res[i] = r > 0 ? r : NaN

        # if new_r_res[i]  >=  d_res[i]*(B_eq[i] + pred_sum)
        #     println("Resources are fine and r/d is $(new_r_res[i]/d_res[i])")
        # else
        #     println("Resources are not fine and r/d is $(new_r_res[i]/d_res[i])")
        # end
    end

    return new_xi_cons, new_r_res
end

# -----------------------------------------------------------
# Main simulation pipeline
# -----------------------------------------------------------
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
                    ε_full = clamp.(rand(Normal(0.2, 0.1), total, total), 0, 1)

                    # 4) calibrate ONCE with full A & ε
                    p_cal = (R, C, m_cons, d_res, ε_full, A_adj)
                    xi_cons, r_res = cp(R_eq, C_eq, p_cal)

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
                        xi_cons, r_res = cp(R_eq, C_eq, p_cal)

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

# -----------------------------------------------------------
# Ladder‐of‐simplification transforms
# -----------------------------------------------------------
function transform_for_ladder_step(step, A_adj, ε_full)
    total = size(A_adj,1)

    if step == 1
        return A_adj, ε_full
    elseif step == 2
        ε2 = similar(ε_full)
        for i in 1:total
            ε2[i,:] .= mean(ε_full[i,:])
        end
        return A_adj, ε2
    elseif step == 3
        return A_adj, fill(mean(ε_full), total, total)
    elseif step == 4
        ε4 = similar(ε_full)
        for i in 1:total, j in 1:total
            ε4[i, j] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
        end
        return A_adj, ε4
    elseif 5 ≤ step ≤ 8
        # 5: ε_full; 6: random per row; 7: global .5
        pos, neg = mean(A_adj[A_adj.>0]), mean(A_adj[A_adj.<0])
        # println("pos: $pos, neg: $neg for step $step")
        # if isnan(neg)
        #     println("Negative values of A are $(A_adj[A_adj.<0])")
        # end
        A5 = map(x-> x>0 ? pos : x<0 ? neg : 0, A_adj)
        if step==5
            return A5, ε_full
        elseif step==6
            ε6 = similar(ε_full)
            for i in 1:total
                ε6[i,:] .= mean(ε_full[i,:])
            end
            return A5, ε6
        elseif step==7
            ε7 = fill(mean(ε_full), total, total)
            return A5, ε7
        elseif step==8
            ε8 = similar(ε_full)
            for i in 1:total, j in 1:total
                ε8[i, j] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
            end
            return A5, ε8
        end
    elseif 9 ≤ step ≤ 12
        # 9: ε_full; 10: Average per row; 11: Global .5; 12: re-randomised
        m = mean(abs.(A_adj[A_adj .!= 0]))
        A6 = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)
        # println("m: $m for step $step")
        if step==9
            return A6, ε_full
        elseif step==10
            ε10 = similar(ε_full)
            for i in 1:total
                ε10[i,:] .= mean(ε_full[i,:])
            end
            return A6, ε10
        elseif step==11
            return A6, fill(mean(ε_full), total, total)
        else
            ε12 = similar(ε_full)
            for i in 1:total*total
                ε12[i] = clamp(rand(Normal(mean(ε_full), std(ε_full))), 0, 1)
            end
            return A6, ε12
        end
    elseif 13 ≤ step ≤ 16
        # Fully randomized A matrix using Normal(mean, std) of abs non-zero A entries
        A_vals = abs.(A_adj[A_adj .!= 0])
        A_mean, A_std = mean(A_vals), std(A_vals)
        # println("Randomised A from N($A_mean, $A_std) for step $step")
    
        # Randomise A with same sparsity pattern as A_adj
        A_rand = similar(A_adj)
        for i in 1:total, j in 1:total
            if A_adj[i, j] != 0
                A_rand[i, j] = rand(Normal(A_mean, A_std)) * sign(A_adj[i, j])
            else
                A_rand[i, j] = 0.0
            end
        end
    
        if step == 13
            return A_rand, ε_full
        elseif step == 14
            ε14 = similar(ε_full)
            for i in 1:total
                ε14[i, :] .= mean(ε_full[i, :])
            end
            return A_rand, ε14
        elseif step == 15
            return A_rand, fill(mean(ε_full), total, total)
        else
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

begin
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
serialize("Daily_JULIA_Code/ThePaper/Ladder/new_results3.jls", new_results3)

println("First 10 rows of simulation results:")
println(first(df_results, 10))

"""
    build_jacobian(B_eq, p)

Analytically constructs the (R+C)×(R+C) Jacobian at equilibrium B_eq,
given p = (R,C,m_cons,xi_cons,r_res,d_res,ε,A).
"""
function build_jacobian(B_eq, p)
    R, C, m_cons, xi_cons, r_res, d_res, ε, A = p
    total = R + C
    J = zeros(total, total)

    # Resource block
    for i in 1:R
        # ∂(du_i)/∂B_i = - d_res[i] * B_eq[i]
        J[i,i] = -d_res[i]*B_eq[i]
        # ∂(du_i)/∂B_j (j a consumer)
        for j in (R+1):total
            J[i,j] =  d_res[i]*B_eq[i]*A[i,j]
        end
    end

    # Consumer block
    for k in 1:C
        i = R + k
        ψ = B_eq[i] / xi_cons[k]
        prefactor = m_cons[k] * ψ
        for j in 1:total
            # -δ_{ij} + ε[i,j]*A[i,j]  - A[j,i]
            J[i,j] = prefactor * ( (i==j ? -1 : 0) + ε[i,j]*A[i,j] - A[j,i] )
        end
    end

    return J
end

"""
    compute_analytical_V(J, R, C, m_cons, xi_cons)

Return V = J^{-1} * D, with D_{ii}=m_i/xi_i for consumers (i>R), D_{ii}=0 for resources.
"""
function compute_analytical_V(J, R, C, m_cons, xi_cons)
    total = size(J,1)
    D = zeros(total, total)
    for k in 1:C
        D[R+k, R+k] = m_cons[k] / xi_cons[k]
    end
    # Solve J * V = D  ⇒  V = J \ D
    return J \ D
end
