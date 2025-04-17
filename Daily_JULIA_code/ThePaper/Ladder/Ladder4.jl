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
            abs(A_used[i,j]) * B_eq[j]
            for j in 1:total if A_used[i,j] < 0.0;
            init = 0.0
        )

        ξ = -C_eq[k] + prey_sum - pred_sum
        new_xi_cons[k] = ξ > 0 ? ξ : NaN
    end

    # 2) Resource growth rates
    for i in 1:R
        # Consumer rows are R+1:total; losses where those entries <0
        pred_sum = sum(
            abs(A_used[i,j]) * B_eq[j]
            for j in 1:R+C if A_used[i,j] < 0.0;
            init = 0.0
        )
        r = d_res[i] * (R_eq[i] + pred_sum)
        new_r_res[i] = r > 0 ? r : NaN
    end

    return new_xi_cons, new_r_res
end

# -----------------------------------------------------------
# Main simulation pipeline
# -----------------------------------------------------------
function rp(
    R_vals, C_vals, conn_vals, delta_vals;
    ladder_steps=1:10,
    Niter=10,
    tspan=(0.0,50.0),
    t_perturb=25.0,
    max_calib=10,
    plot_full = false,
    plot_simple = false,
    atol = 10.0,
    d_value = 1.0,
    IS_reduction = 1.0,
    homogenise = nothing
)
    results = Vector{NamedTuple}()
    # lock = ReentrantLock()

    # suffix helper
    suffix(step) = step==1 ? "Full" : "S$(step-1)"

    for R in R_vals, C in C_vals, conn in conn_vals, delta in delta_vals
        total = R+C

        for iter in 1:Niter
            # println("run")

            A_adj = zeros(total,total)
            for i in (R+1):total, j in 1:total
                if i !=j && rand() < conn && iszero(A_adj[i,j])
                    A_adj[i,j]= rand() * IS_reduction
                    A_adj[j,i]= -rand() * IS_reduction
                end
            end
            
            # 1) target equilibrium
            fixed = abs.(rand(LogNormal(0.5, 1.0), total))
            R_eq, C_eq = fixed[1:R], fixed[R+1:end]

            # 2) fixed params
            m_cons = fill(0.1, C)
            d_res  = fill(d_value, R)

            # 3) full ε
            ε_full = clamp.(rand(Normal(0.3, 0.1), total, total), 0, 1)

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
                        A_adj[j,i]= -rand() * IS_reduction
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

                fixed = abs.(rand(LogNormal(0.5,1.0), total))
                p_cal = (R, C, m_cons, d_res, ε_full, A_adj)
                
                R_eq, C_eq = fixed[1:R], fixed[R+1:end]
                xi_cons, r_res = cp(R_eq, C_eq, p_cal)

                tries+=1
            end

            if any(isnan,xi_cons)||any(isnan,r_res)
                @warn "Calibration failed; skipping iteration"
                continue
            end

            # dict to hold all step‐metrics
            metrics = Dict{String,NamedTuple}()

            # run FULL first (step=1)
            p_full = (R, C, m_cons, xi_cons, r_res, d_res, ε_full, A_adj)

            prob = ODEProblem(trophic_ode!, fixed, tspan, p_full)
            sol = solve(prob, reltol=1e-8,abstol=1e-8)
            println("\n", any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=1e-3) for i in 1:total]))
            if sol.t[end] < t_perturb || any(isnan, sol.u[end]) || any(isinf, sol.u[end]) || any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:total])
                @warn "Error: solution did not finish properly"
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
            metrics["Full"] = (
                persistence=pers,
                return_time=mean(skipmissing(rt)),
                before_p=pers,
                after_p=mean(B2 .> EXTINCTION_THRESHOLD),
                resilience=resi,
                reactivity=reac
            )

            # now each simplified step
            for step in ladder_steps[2:end]
                suf = suffix(step)
                A_used, ε_used = transform_for_ladder_step(step, A_adj, ε_full)
                p_simp = (R,C,m_cons,xi_cons,r_res,d_res,ε_used,A_used)
                sol2 = solve(
                    ODEProblem(trophic_ode!, fixed, tspan, p_simp), Tsit5(),
                    reltol=1e-8,abstol=1e-8
                )
                Beq2 = sol2.u[end]
                pers2 = mean(Beq2 .> EXTINCTION_THRESHOLD)
                res2  = compute_resilience(Beq2,p_simp)
                rea2  = compute_reactivity(Beq2,p_simp)
                rt2,_,B3 = simulate_press_perturbation(
                    fixed, p_simp, tspan, t_perturb, delta;
                    solver=Tsit5(), plot=plot_simple,
                    show_warnings=true,
                    full_or_simple = false
                )
                metrics[suf] = (
                    persistence=pers2,
                    return_time=mean(skipmissing(rt2)),
                    before_p=pers2,
                    after_p=mean(B3 .> EXTINCTION_THRESHOLD),
                    resilience=res2,
                    reactivity=rea2
                )
            end

            # flatten into one NamedTuple
            base = (
                resource_count=R,
                consumer_count=C,
                connectance=conn,
                perturb_delta=delta,
                iteration=iter
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
                    Symbol("reactivity_$suf")        => m.reactivity)
                )
            end

            # lock(lock) do
            push!(results, rec)
            # end
        end
    end

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
        ε3 = zeros(total,total)
        for i in 1:total
            ε3[i,:] .= clamp.(rand(Normal(0.3, 0.1), total), 0, 1)
        end
        return A_adj, ε3
    elseif step == 4
        return A_adj, fill(mean(ε_full), total, total)
    elseif 5 ≤ step ≤ 7
        # 5: ε_full; 6: random per row; 7: global .5
        pos, neg = mean(A_adj[A_adj.>0]), mean(A_adj[A_adj.<0])
        A5 = map(x-> x>0 ? pos : x<0 ? neg : 0, A_adj)
        if step==5
            return A5, ε_full
        elseif step==6
            ε6 = rand(total,total)
            return A5, ε6
        else
            return A5, fill(0.5, total, total)
        end
    elseif 8 ≤ step ≤ 10
        # 8: ε_full; 9: random; 10: global .5
        m = mean(abs.(A_adj[A_adj .!= 0]))
        A6 = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)
        if step==8
            return A6, ε_full
        elseif step==9
            return A6, rand(total,total)
        else
            return A6, fill(0.5, total, total)
        end
    else
        error("Unknown ladder step $step")
    end
end

begin
# for i in 1:20 
    # Specify parameter ranges:
    R_vals = [20]            # e.g., try 4 and 6 resources
    C_vals = [5]            # e.g., try 5 and 8 consumers
    connectance_list = 0.1:0.1:0.5
    delta_list = 0.1:0.1:0.2

    new_results = rp(
        R_vals, C_vals, connectance_list, delta_list;
        ladder_steps=1:9,
        Niter=3,
        tspan=(0.0, 500.0),
        t_perturb=250.0,
        max_calib=50,
        plot_full=true,
        plot_simple=true,
        atol = 10.0,
        d_value = 1.0,
        IS_reduction = 1.0,
        homogenise = 0.3
    )
end

println("First 10 rows of simulation results:")
println(first(df_results, 10))