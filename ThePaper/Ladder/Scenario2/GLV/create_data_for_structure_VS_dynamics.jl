function exploring(
    S::Int=50, C::Int=20;
    conn_vals=[0.05, 0.1, 0.2],
    scenarios=[:ER, :PL, :MOD],
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    delta_vals=[1.0, 3.0],
    margins=[1.0],
    number_of_combinations::Int=100,
    iterations::Int=1,
    tspan=(0.0, 500.0),
    tpert::Float64=250.0,
    pareto_exponents=[1.25,1.75,2.0,3.0,4.0,5.0],
    pareto_minimum_degrees=[1.0],
    mod_gammas = [1.0,2.0,3.0,5.0,10.0]
)
    R = S - C
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    combos = collect(Iterators.product(conn_vals, scenarios, IS_vals, delta_vals, margins, 1:iterations, pareto_exponents, pareto_minimum_degrees, mod_gammas))
    @info "Computing $(length(combos)) combinations"
    global cb = build_callbacks(50, 1e-6)
    @threads for (conn, scen, IS, delta, marg, ite, pex, p_min_deg, mod_gamma) in
            sample(combos, min(length(combos), number_of_combinations); replace=false)

        # 1) build A & callback
        A  = make_A(
            zeros(S,S), R, conn, scen; IS=IS,
            pareto_exponent=pex,pareto_minimum_degree=p_min_deg,mod_gamma=mod_gamma
        )
        
        # 2) collectivity
        phi = compute_collectivity(A)

        # 3) find feasible K & equilibrium
        thr_sets = generate_feasible_thresholds(A, R; margins=[marg])
        isempty(thr_sets) && continue
        tset = thr_sets[1]
        K    = tset.K
        u_eq = tset.u_eq

        # 4) stability & survival
        ok, u0 = survives!(u_eq, (K,A); cb=cb)
        !ok && continue
        J = compute_jacobian_glv(u0, (K,A))
        # J    = D * M
        !is_locally_stable(J) && continue

        # 5) full-model metrics
        S_full           = count(x->x>1e-6, u0)
        resilience_full  = compute_resilience(u0, (K,A))
        reactivity_full  = compute_reactivity(u0, (K,A))
        # SL_full          = compute_SL(A, K)
        SL_full         = diag(J) .|> x -> x == 0.0 ? 0.0 : -1 / x
        mean_SL_full    = mean(SL_full)
        sigma_full       = sigma_over_min_d(A, J)

        # 5a) press perturbation
        rt_press_vec, before_press, after_press, _ =
            simulate_press_perturbation_glv(u0, (K,A), tspan, tpert, delta, R; cb=cb)
        rt_press_full   = mean(skipmissing(rt_press_vec))

        # 5b) pulse perturbation
        rt_pulse_vec, before_pulse, after_pulse, _ =
            simulate_pulse_perturbation_glv(u0, (K,A), tspan, tpert, delta; cb=cb)
        rt_pulse_full   = mean(skipmissing(rt_pulse_vec))

        rmed_full = analytical_median_return_rate(J; t=1.0)

        # assemble record
        rec = (
            conn=conn, scen=scen, IS=IS, delta=delta, marg=marg, ite=ite,
            S_full=S_full,
            resilience_full=resilience_full,
            reactivity_full=reactivity_full,
            collectivity_full=phi,
            SL_full=SL_full,
            mean_SL_full=mean_SL_full,
            sigma_over_min_d_full=sigma_full,
            rmed_full=rmed_full,
            after_press_full=after_press,
            after_pulse_full=after_pulse,
            rt_press_full=rt_press_full,
            rt_pulse_full=rt_pulse_full,
            p_final=(K,A),
            B_eq = u0,
            Beq_cv = std(u0) / mean(u0)
        )

        lock(locki) do
            push!(results, rec)
        end
    end

    @info "Finished computing $number_of_combinations combinations"
    return DataFrame(results)
end

function run_all_exploring()
    R_all = DataFrame()
    for i in [(100, 40), (50, 20), (200, 80), (300, 120)]
        a, b = i[1], i[2]
        R = exploring(
            a, b;
            conn_vals=0.25,
            scenarios=[:ER, :PL, :MOD],
            IS_vals=[0.1],
            delta_vals=[0.1],
            margins=[1.0, 2.0, 3.0, 4.0, 5.0, 0.01],
            number_of_combinations=100,
            iterations=10,
            pareto_exponents=[1.0, 1.25, 1.75, 2.0, 3.0, 4.0, 5.0],
            pareto_minimum_degrees=[5.0, 10.0, 15.0, 20.0],
            mod_gammas=[1.0,2.0,3.0,5.0,10.0]
        )
        R_all = vcat(R_all, R)
    end
    return R_all
end

R_all = run_all_exploring()
serialize("data_for_structureVsDynamics.jls", R_all)
R_all = deserialize("data_for_structureVsDynamics.jls")
# serialize("checking_glv_50000perSpeciesN.jls", R_all)

# include("orderColumns.jl")
# R = reorder_columns(R)
R = deserialize("checking_glv_50000ALL.jls")
R = deserialize("checking_glv_50000perSpeciesN.jls")
# R = deserialize("checking_10000ALL.jls")

# Assuming df1 and df2 are already defined DataFrames with the same column names
# merged_df = vcat(G, T)
