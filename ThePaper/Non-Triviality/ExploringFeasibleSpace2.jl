# using DifferentialEquations, Random, Statistics, DataFrames, CSV
# import Base.Threads: @threads, nthreads, threadid
# include("../Ladder/Ladder4.1_drago1.jl")
# --- helper to build A according to scenario ---
function make_A(
    A::AbstractMatrix{<:Real},
    R::Int,
    conn::Float64,
    scenario::Symbol,
    param::Float64
)
    S = size(A,1)
    @assert 1 ≤ R < S "Must have at least one resource and one consumer"
    C = S - R

    # zero‐out
    fill!(A, 0.0)

    if scenario == :ER
        # each possible consumer→any link with prob=conn
        for i in (R+1):S, j in 1:S
            if i != j && rand() < conn
                A[i,j] = abs(rand(Normal()))
                A[j,i] = -abs(rand(Normal()))
            end
        end
    elseif scenario == :PL
        # power‐law consumer out‐degree: only consumers get ks
        α     = 3.0
        k_max = S-1
        # sample desired out‐degrees
        raw   = rand(Pareto(1, α), C)
        ks    = clamp.(round.(Int, raw), 1, k_max)

        for (idx, k) in enumerate(ks)
            ci = R + idx  # consumer index
            # pick k distinct prey‐nodes from 1:S except ci
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

# --- modified feasibility_search2 ---
function feasibility_search2(
    S_vals, conn_vals, C_ratios, IS_vals,
    d_vals, m_vals, eps_vals, delta_vals;
    # now pass scenarios as (Symbol,Float) tuples:
    network_scenarios = [(:ER, 0.0), (:PL, 2.5), (:MOD, 2.0)],
    pyramid_skews     = [0.5, 1.0, 2.0],
    Niter=10, tspan=(0.0,500.0), t_perturb=25.0,
    max_calib=10, abundance_mean=1.0,
    atol=1e-3, xi_threshold=0.7,
    number_of_combos=1000
)
    combos = collect(Iterators.product(
      S_vals, conn_vals, C_ratios, IS_vals,
      d_vals, m_vals, eps_vals, delta_vals,
      network_scenarios, pyramid_skews
    ))
    T = nthreads()
    local_recs = [Vector{NamedTuple}() for _ in 1:T]

    @threads for idx in sample(1:length(combos), min(number_of_combos,length(combos)), replace=false)
        tid = threadid()
        S, conn, C_ratio, IS, d_val, m_val, epsilon_mean, delta, (scenario,param), skew = combos[idx]
        C = clamp(round(Int,S*C_ratio),1,S-1)
        R = S - C
        if S == 25
            cb = cb_no_trigger25
        elseif S == 50
            cb = cb_no_trigger50
        elseif S == 75
            cb = cb_no_trigger75
        elseif S == 100
            cb = cb_no_trigger100
        end
        for iter in 1:Niter
            # 1) build A
            A = zeros(S,S)
            A = make_A(A, R,conn, scenario, param)

            # 2) draw eq with pyramid skew
            cons_CV = 0.2
            res_CV  = skew * cons_CV
            R_eq = abs.(rand(LogNormal(log(abundance_mean) - res_CV^2/2, res_CV), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1) - cons_CV^2/2, cons_CV), C))
            fixed = vcat(R_eq, C_eq)

            # --- rates ---
            m_cons = fill(m_val, C)
            d_res   = fill(d_val, R)
            epsilon_mat   = clamp.(rand(LogNormal(epsilon_mean, epsilon_mean), S, S), 0, 1)

            # --- calibrate with up to max_calib retries ---
            pcal = (R, C, m_cons, d_res, epsilon_mat, A)
            xi_cons, r_res = calibrate_params(R_eq, C_eq, pcal; xi_threshold=xi_threshold)
            
            tries = 1
            while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
                A .= 0
                A = make_A(A, R, conn, scenario, param)
                R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
                C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1) - ((abundance_mean*0.2)^2)/2, abundance_mean*0.2), C))
                fixed = vcat(R_eq, C_eq)
                pcal = (R, C, m_cons, d_res, epsilon_mat, A)
                xi_cons, r_res = calibrate_params(R_eq, C_eq, pcal; xi_threshold=xi_threshold)
                tries += 1
            end

            # if calibration failed, record infeasible
            if any(isnan, xi_cons) || any(isnan, r_res)
                push!(local_recs[tid], (
                    S=S, conn=conn, C_ratio=C_ratio,
                    IS=IS, d=d_val, m=m_val,
                    epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                    feasible=false, before_p = 0.0, after_p = 0.0, scenario=scenario, param=param, skew=skew
                ))
                continue
            end

            # --- simulate the unperturbed dynamics ---
            pfull = (R, C, m_cons, xi_cons, r_res, d_res, epsilon_mat, A)
            prob   = ODEProblem(trophic_ode!, fixed, tspan, pfull)
            sol    = solve(prob, Tsit5(); callback=cb, reltol=1e-8, abstol=1e-8)

            # --- feasibility check ---
            if sol.t[end] < t_perturb || any(isnan, sol.u[end]) || any(isinf, sol.u[end]) || any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:S])
                ok = false
                push!(local_recs[tid], (
                S=S, conn=conn, C_ratio=C_ratio,
                IS=IS, d=d_val, m=m_val,
                epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                feasible=ok, before_p = 0.0, after_p = 0.0, scenario=scenario, param=param, skew=skew
            ))
            else
                ok = true

                Beq = sol.u[end]
                before_pers = mean(Beq .> EXTINCTION_THRESHOLD)

                rt_full, os_full, ire_full, _, B2 = simulate_press_perturbation(
                            fixed, pfull, (0.0, 500.0), 250.0, delta;
                            solver=Tsit5(), plot=false,
                            show_warnings=true,
                            full_or_simple = true,
                            cb = cb
                        )

                push!(local_recs[tid], (
                    S=S, conn=conn, C_ratio=C_ratio,
                    IS=IS, d=d_val, m=m_val,
                    epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                    feasible=ok, before_p = before_pers, after_p = mean(B2 .> EXTINCTION_THRESHOLD), scenario=scenario, param=param, skew=skew
                ))
            end
        end
    end

    return DataFrame(reduce(vcat, local_recs))
end

if false
# cb_no_trigger25 = build_callbacks(25, EXTINCTION_THRESHOLD)
# serialize("ThePaper/Callbacks_library/cb_no_trigger25.jls", cb_no_trigger25)
cb_no_trigger25 = deserialize("ThePaper/Callbacks_library/cb_no_trigger25.jls")
# cb_no_trigger50 = build_callbacks(50, EXTINCTION_THRESHOLD)
# serialize("ThePaper/Callbacks_library/cb_no_trigger50.jls", cb_no_trigger50)
cb_no_trigger50 = deserialize("ThePaper/Callbacks_library/cb_no_trigger50.jls")
# cb_no_trigger75 = build_callbacks(75, EXTINCTION_THRESHOLD)
# serialize("ThePaper/Callbacks_library/cb_no_trigger75.jls", cb_no_trigger75)
cb_no_trigger75 = deserialize("ThePaper/Callbacks_library/cb_no_trigger75.jls")
# cb_no_trigger100 = build_callbacks(100, EXTINCTION_THRESHOLD)
# serialize("ThePaper/Callbacks_library/cb_no_trigger100.jls", cb_no_trigger100)
cb_no_trigger100 = deserialize("ThePaper/Callbacks_library/cb_no_trigger100.jls")
end

begin
    S_vals      = [25, 50, 75, 100]
    conn_vals   = 0.1:0.2:1.0
    C_ratios    = 0.1:0.2:0.9
    IS_vals     = [0.001, 0.01, 0.1, 1.0, 2.0]
    d_vals      = [0.01, 0.1, 1.0, 2.0]
    mort_vals   = [0.05, 0.1, 0.2, 0.5]
    epsilon_vals= [0.01, 0.1, 0.2, 0.5]
    delta_vals  = [0.1]
    network_scenarios=[(:PL,1.0), (:PL,2.5), (:PL,5.0), (:PL,10.0)]

    A_feas = feasibility_search2(
        S_vals, conn_vals, C_ratios, IS_vals,
        d_vals, mort_vals, epsilon_vals, delta_vals;
        network_scenarios=network_scenarios, 
        pyramid_skews=[0.5,1.0,2.0],
        Niter=1,
        tspan=(0.0,100.0),
        t_perturb=50.0,
        max_calib=10,
        abundance_mean=1.0,
        atol=1.0,
        xi_threshold=0.7,
        number_of_combos=10000 # 10000 is good 
    )
    # (b) Extract good combos
    A_good = unique(
        A_feas[A_feas.feasible .== true,
                [:S, :conn, :C_ratio, :IS, :d, :m, :epsilon, :R, :C, :scenario, :param, :skew]]
    )
end