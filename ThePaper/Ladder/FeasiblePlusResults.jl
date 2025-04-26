include("Ladder4.1_drago1.jl")
using DifferentialEquations, Random, Statistics, DataFrames, CSV
import Base.Threads: @threads, nthreads, threadid
include("ExploringFeasibleSpace.jl")

function make_A(
    A::AbstractMatrix{<:Real},
    R::Int,
    conn::Float64,
    scenario::Symbol
)
    S = size(A,1)
    @assert 1 <= R < S "Must have at least one resource and one consumer"
    C = S - R

    # zero-out
    fill!(A, 0.0)

    if scenario == :ER
        # each possible consumer?any link with prob=conn
        for i in (R+1):S, j in 1:S
            if i != j && rand() < conn && iszero(A[i,j])
                A[i,j] = abs(rand(Normal()))
                A[j,i] = -abs(rand(Normal()))
            end
        end
    elseif scenario == :PL
        # power-law consumer out-degree: only consumers get ks
        a     = 1.75
        k_max = S-1
        # sample desired out-degrees
        raw   = rand(Pareto(1, a), C)
        ks    = clamp.(round.(Int, raw), 1, k_max)

        for (idx, k) in enumerate(ks)
            ci = R + idx  # consumer index
            # pick k distinct prey-nodes from 1:S except ci
            js = sample(filter(x->x!=ci, 1:S), k; replace=false)
            for j in js
                A[ci, j] = abs(rand(Normal()))
                A[j,  ci] = -abs(rand(Normal()))
            end
        end
    else
        error("Unknown scenario $scenario")
    end

    return A
end

# -----------------------------------------------------------
# 3) rppp_guided: only over pre-screened param tuples
# -----------------------------------------------------------
function rppp_guided(df_params::DataFrame, delta_vals;
                     ladder_steps=1:16,
                     Niter=10, tspan=(0.0,50.0), t_perturb=25.0,
                     max_calib=10, plot_full=false, plot_simple=false,
                     atol=10.0,
                     abundance_mean=1.0,
                     combinations=1000
)

    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()
    suffix(step) = step==1 ? "Full" : "S$(step-1)"

    @threads for row in eachrow(df_params)[sample(1:nrow(df_params), min(combinations, nrow(df_params)), replace=false)]
        S        = row.S
        conn     = round(row.conn, digits=3)
        C_ratio  = row.C_ratio
        IS_red   = row.IS
        d_val    = row.d
        m_val    = row.m
        epsilon_mean   = row.epsilon
        C        = row.C
        R        = row.R
        scenario = row.scenario
        constraints = row.constraints
        ssp = row.ssp

        if S == 25
            cb = cb_no_trigger25
        elseif S == 50
            cb = cb_no_trigger50
        elseif S == 75
            cb = cb_no_trigger75
        elseif S == 100
            cb = cb_no_trigger100
        end
        for delta in delta_vals
            for iter in 1:Niter
                try
                    # build A & epsilon_full
                    A = zeros(S,S)
                    # for i in (R+1):S, j in 1:S
                    #     if i!=j && rand()<conn
                    #         A[i,j]=abs(rand(Normal(0,IS_red)))
                    #         A[j,i]=-abs(rand(Normal(0,IS_red)))
                    #     end
                    # end

                    A = make_A(A, R, conn, scenario)

                    epsilon_full = clamp.(rand(LogNormal(epsilon_mean, epsilon_mean), S, S), 0, 1)

                    # target eq
                    R_eq = abs.(rand(LogNormal(log(abundance_mean)-abundance_mean^2/2, abundance_mean), R))
                    C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1)-(abundance_mean*0.2)^2/2, abundance_mean*0.2), C))
                    fixed = vcat(R_eq, C_eq)

                    m_cons  = fill(m_val, C)
                    d_res   = fill(d_val, R)

                    # calibrate
                    xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,epsilon_full,A); xi_threshold=0.7, constraints=constraints)
                    tries=1
                    while (any(isnan, xi_cons)||any(isnan, r_res)) && tries<max_calib
                        A .= 0
                        # for i in (R+1):S, j in 1:S
                        #     if i!=j && rand()<conn
                        #         A[i,j]=abs(rand(Normal(0,IS_red)))
                        #         A[j,i]=-abs(rand(Normal(0, IS_red)))
                        #     end
                        # end
                        A = make_A(A, R, conn, scenario)

                        R_eq = abs.(rand(LogNormal(log(abundance_mean)-abundance_mean^2/2, abundance_mean), R))
                        C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1)-(abundance_mean*0.2)^2/2, abundance_mean*0.2), C))
                        fixed = vcat(R_eq, C_eq)
                        xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,epsilon_full,A); xi_threshold=0.7, constraints=constraints)
                        tries+=1
                    end
                    if any(isnan, xi_cons)||any(isnan, r_res)
                        continue
                    end

                    # fullâ€model solve & metrics
                    p_full = (R,C,m_cons,xi_cons,r_res,d_res,epsilon_full,A)

                    prob   = ODEProblem(trophic_ode!, fixed, tspan, p_full)
                    sol    = solve(prob,Tsit5(); callback=cb, reltol=1e-8, abstol=1e-8)
                    if sol.t[end] < t_perturb || any(isnan, sol.u[end]) || any(isinf, sol.u[end]) || any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:S])
                        @warn "Error: solution did not finish properly"
                        continue
                    end
                    B_eq = sol.u[end]

                    g = SimpleGraph(A .!= 0)
                    degs = degree(g)
                    degree_cv = std(degs) / mean(degs)

                    meanB = mean(B_eq)
                    varB = var(B_eq)
                    RelVar = varB / (meanB^2 + 1e-8)
                    
                    resi = compute_resilience(B_eq, p_full)
                    reac = compute_reactivity(B_eq, p_full)
                    # press perturb
                    rt_full, os_full, ire_full, _, B2, sp_rt_full = simulate_press_perturbation(
                        fixed, p_full, tspan, t_perturb, delta;
                        solver=Tsit5(), plot=plot_full,
                        show_warnings=true, full_or_simple=true,
                        cb = cb,
                        species_specific_perturbation=ssp
                    )
                    # sensitivity corr
                    J_full   = build_jacobian(fixed, p_full)
                    V_ana    = compute_analytical_V(J_full, R, C, m_cons, xi_cons)
                    pred     = sum(V_ana,dims=2)
                    obs      = (B2 .- fixed)/delta
                    sens_corr_full = cor(vec(pred), vec(obs))

                    # record metrics for each step
                    metrics = Dict("Full" => (
                        return_time=mean(filter(!isnan, rt_full)),
                        overshoot  =mean(filter(!isnan, os_full)),
                        ire        =mean(filter(!isnan, ire_full)),
                        resilience =resi,
                        reactivity =reac,
                        sens_corr  =sens_corr_full,
                        before_p   =mean(sol.u[end] .> 1e-6),
                        after_p    =mean(B2 .> 1e-6),
                        sp_rt_mean     =mean(filter(!isnan, sp_rt_full)),
                        sp_rt_cv       =std(filter(!isnan, sp_rt_full)) / mean(filter(!isnan, sp_rt_full))
                        ))

                    # ladder steps 
                    for step in ladder_steps[2:end]
                        suf = suffix(step)
                        A_s, epsilon_s = transform_for_ladder_step(step, A, epsilon_full)
                        p_simp   = (R,C,m_cons,xi_cons,r_res,d_res,epsilon_s,A_s)

                        J_simp         = build_jacobian(fixed, p_simp)
                        V_ana_simp     = compute_analytical_V(J_simp, R, C, m_cons, xi_cons)

                        sol2  = solve(
                            ODEProblem(
                                trophic_ode!, fixed, tspan, p_simp
                            ),
                            Tsit5(); 
                            callback=cb, reltol=1e-8, abstol=1e-8
                        )
                        if sol2.t[end] < t_perturb || any(isnan, sol2.u[end]) || any(isinf, sol2.u[end])
                            metrics[suf] = (
                                return_time=NaN, overshoot=NaN, ire=NaN,
                                resilience=NaN, reactivity=NaN,
                                sens_corr=NaN, before_p=NaN, after_p=NaN
                            )
                            continue
                        end

                        B_eq2 = sol2.u[end]
                        before_p = mean(sol2.u[end] .> EXTINCTION_THRESHOLD)
                        res2  = compute_resilience(B_eq2,p_simp)
                        rea2  = compute_reactivity(B_eq2,p_simp)
                        # get the full per-species press matrix just like in the full model
                        rt2, os2, ire2, _, B2_simp, sp_rt_simp = simulate_press_perturbation(
                                fixed, p_simp, tspan, t_perturb, delta;
                                solver=Tsit5(), plot=plot_simple,
                                show_warnings=true,
                                full_or_simple=true,
                                cb = cb,
                                species_specific_perturbation=ssp
                            )
        
                        # in each step, after computing V_ana_simp and B2_simp
                        pred_resp_simp = sum(V_ana_simp, dims=2)
                        obs_resp_simp  = (B2_simp .- fixed) ./ delta
                        corr_simp      = cor(vec(pred_resp_simp), vec(obs_resp_simp))

                        # compound error, step
                        # t_pred2 = mean(filter(!isnan, rt2))
                        # t_obs2  = mean(filter(!isnan, rt2))
                        # err_T2   = abs(t_obs2 - t_pred2) / ((t_obs2 + t_pred2)/2)
                        # os_pred2 = mean(filter(!isnan, os_full))
                        # os_obs = mean(filter(!isnan, os2))
                        # err_OS2  = abs(os_obs - os_pred2) / ((os_obs + os_pred2)/2)
                        # ire_pred2 = mean(filter(!isnan, ire_full))
                        # ire_obs = mean(filter(!isnan, ire2))
                        # err_IRE2  = abs(ire_obs - ire_pred2) / ((ire_obs + ire_pred2)/2)
                        # comp_err2 = (err_T2 + err_OS2 + err_IRE2)/3

                        metrics[suf] = (
                            # compound_error = comp_err2,
                            return_time    = mean(filter(!isnan, rt2)),
                            overshoot      = mean(filter(!isnan, os2)),
                            ire            = mean(filter(!isnan, ire2)),
                            resilience     = res2,
                            reactivity     = rea2,
                            sens_corr      = corr_simp,
                            before_p       = before_p,
                            after_p        = mean(B2_simp .> EXTINCTION_THRESHOLD),
                            sp_rt_mean         = mean(filter(!isnan, sp_rt_simp)),
                            sp_rt_cv           = std(filter(!isnan, sp_rt_simp)) / mean(filter(!isnan, sp_rt_simp))ºº
                        )
                    end

                    # assemble record
                    base = (
                        S=S, conn=conn, C_ratio=C_ratio,
                        IS=IS_red, d=d_val, m=m_val, epsilon=epsilon_mean,
                        delta=delta, iteration=iter, R=R, C=C, degree_cv=degree_cv, RelVar=RelVar, 
                        scenario=scenario, constraints=constraints, ssp=ssp
                    )
                    rec = base
                    for step in ladder_steps
                        suf = suffix(step)
                        m = metrics[suf]
                        rec = merge(rec,
                            (; Symbol("rt_$suf")      => m.return_time,
                            Symbol("os_$suf")      => m.overshoot,
                            Symbol("ire_$suf")     => m.ire,
                            Symbol("res_$suf")     => m.resilience,
                            Symbol("rea_$suf")     => m.reactivity,
                            Symbol("scorr_$suf")   => m.sens_corr,
                            Symbol("bef_$suf")     => m.before_p,
                            Symbol("aft_$suf")     => m.after_p,
                            Symbol("sp_rt_$suf")   => m.sp_rt_mean,
                            Symbol("sp_rt_cv_$suf")   => m.sp_rt_cv)
                        )
                    end

                    lock(results_lock) do
                        push!(results, rec)
                    end

                catch e
                    @warn "Error in rppp_guided" row delta iter e
                end
            end  # iter
        end  # delta
    end  # threads

    return DataFrame(results)
end
# -----------------------------------------------------------
# 4) Example invocation
# -----------------------------------------------------------
# cb_no_trigger25 = build_callbacks(25, EXTINCTION_THRESHOLD)
# cb_no_trigger50 = build_callbacks(50, EXTINCTION_THRESHOLD)
cb_no_trigger75 = build_callbacks(75, EXTINCTION_THRESHOLD)
# cb_no_trigger100 = build_callbacks(100, EXTINCTION_THRESHOLD)

begin
    S_vals      = [75]
    conn_vals   = range(0.01, 0.4, length=10)
    C_ratios    = 0.1:0.1:0.5
    IS_vals     = [0.001, 0.01, 0.1, 1.0, 2.0]
    d_vals      = [0.01, 0.1, 1.0, 2.0]
    mort_vals   = [0.05, 0.1, 0.2, 0.5]
    epsilon_vals= [0.01, 0.1, 0.2, 0.5]
    delta_vals  = [0.1]

    df_feas = feasibility_search(
        S_vals, conn_vals, C_ratios, IS_vals,
        d_vals, mort_vals, epsilon_vals, delta_vals;
	    degree_distribution_types=[:ER, :PL],
        Niter=1,
        tspan=(0.0,100.0),
        t_perturb=50.0,
        max_calib=50,
        abundance_mean=1.0,
        atol=10.0,
        xi_threshold=2.0,
        number_of_combos=100000, # 10000 is good 
        calibrate_params_constraints=false,
        species_specific_perturbation=true
    )
end

# (b) Extract good combos
df_good = unique(
    df_feas[df_feas.feasible .== true,
            [:S, :conn, :C_ratio, :IS, :d, :m, :epsilon, :R, :C, :scenario, :constraints]]
)

# (c) Guided rppp
begin
    delta_vals = [0.1]
    df_results = rppp_guided(
        df_good, delta_vals;
        ladder_steps=1:16,
        Niter=1, max_calib=15, # 10 is good
        tspan=(0.0,500.0), t_perturb=250.0,
        atol=10.0,
        abundance_mean=1.0,
        combinations=40000 # You want 10000
    )
end

# (d) Save or inspect
CSV.write("guided_results_ERplusPL1_75_noconstraint_speciesSpecificRTs.csv", df_results)
