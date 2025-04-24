"""
    feasibility_search(
        S_vals, conn_vals, C_ratios, IS_vals,
        d_vals, mortality_vals, epsilon_vals;
        Niter=10, tspan=(0.0,50.0), t_perturb=25.0,
        max_calib=10, abundance_mean=1.0,
        atol=1e-3, xi_threshold=0.7
    )

Same as before but fully threaded:
- builds all (S,conn,C_ratio,IS,d,m,epsilon) combos
- splits them across Julia threads
- each thread collects its results in a local Vector{NamedTuple}
- at the end, we vcat them and build a DataFrame
"""
function feasibility_search(
    S_vals, conn_vals, C_ratios, IS_vals,
    d_vals, mortality_vals, epsilon_vals, delta_vals;
    Niter=10,
    tspan=(0.0,50.0),
    t_perturb=25.0,
    max_calib=10,
    abundance_mean=1.0,
    atol=1e-3,
    xi_threshold=0.7,
    number_of_combos=1000
)
    # 1) build the full list of parameter combinations
    combos = collect(Iterators.product(
      S_vals, conn_vals, C_ratios,
      IS_vals, d_vals, mortality_vals,
      epsilon_vals, delta_vals
    ))
    T = nthreads()
    # one local result buffer per thread
    local_recs = [Vector{NamedTuple}() for _ in 1:T]
    
    @threads for idx in sample(1:length(combos), min(number_of_combos, length(combos)), replace=false)
        tid = threadid()
        (S, conn, C_ratio, IS, d_value, pred_mortality, epsilon_mean , delta) = combos[idx]
        # resource / consumer split
        C = clamp(round(Int, S*C_ratio), 1, S-1)
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
            # --- build interaction matrix ---
            A = zeros(S,S)
            for i in R+1:S, j in 1:S
                if i!=j && rand() < conn

                    A[i,j] =  abs(rand(Normal(0, IS)))
                    A[j,i] = -abs(rand(Normal(0, IS)))
                end
            end

            # --- target equilibrium ---
            R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1) - ( (abundance_mean*0.2)^2)/2, abundance_mean*0.2), C))
            fixed = vcat(R_eq, C_eq)

            # --- rates ---
            m_cons = fill(pred_mortality, C)
            d_res   = fill(d_value, R)
            epsilon_mat   = clamp.(rand(LogNormal(epsilon_mean, epsilon_mean*0.1), S, S), 0, 1)

            # --- calibrate with up to max_calib retries ---
            pcal = (R, C, m_cons, d_res, epsilon_mat, A)
            xi_cons, r_res = calibrate_params(R_eq, C_eq, pcal; xi_threshold=xi_threshold)
            
            tries = 1
            while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
                A .= 0
                for i in R+1:S, j in 1:S
                    if i!=j && rand() < conn && iszero(A[i,j])
                        A[i,j] = abs(rand(Normal(0, IS)))
                        A[j,i] = -abs(rand(Normal(0, IS)))
                    end
                end
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
                    IS=IS, d=d_value, m=pred_mortality,
                    epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                    feasible=false, before_p = 0.0, after_p = 0.0
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
                IS=IS, d=d_value, m=pred_mortality,
                epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                feasible=ok, before_p = 0.0, after_p = 0.0
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
                    IS=IS, d=d_value, m=pred_mortality,
                    epsilon=epsilon_mean, delta=delta, iter=iter, R=R, C=C,
                    feasible=ok, before_p = before_pers, after_p = mean(B2 .> EXTINCTION_THRESHOLD)
                ))
            end
        end
    end

    # 3) gather & return
    allrecs = reduce(vcat, local_recs)
    return DataFrame(allrecs)
end

begin
    S_vals      = [50]
    conn_vals   = 0.1:0.2:1.0
    C_ratios    = 0.1:0.2:0.9
    IS_vals     = [0.001, 0.01, 0.1, 1.0, 2.0]
    d_vals      = [0.01, 0.1, 1.0, 2.0]
    mort_vals   = [0.05, 0.1, 0.2, 0.5]
    epsilon_vals= [0.01, 0.1, 0.2, 0.5]
    delta_vals  = [0.1, 0.2, 0.3, 0.4, 0.5]

    df1 = feasibility_search(
        S_vals, conn_vals, C_ratios, IS_vals,
        d_vals, mort_vals, epsilon_vals, delta_vals;
        Niter=20,
        tspan=(0.0,100.0),
        t_perturb=50.0,
        max_calib=5,
        abundance_mean=1.0,
        atol=1.0,
        xi_threshold=0.7,
        number_of_combos=100
    )
end

percent_feasible = nrow(filter(r -> r.feasible, df1))/nrow(df1)
subset_df = filter(r -> r.feasible, df1)
