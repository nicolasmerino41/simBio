function sigma_over_min_d(A, J)
    d = -diag(J)             # from J = D(-I + A)
    min_d = minimum(d)
    offs = [A[i,j] for i in 1:size(A,1), j in 1:size(A,1) if i!=j]
    σ = std(offs)
    return σ/min_d
end

function test_rewire_wide(
    ; S=50, C=20, conns=0.1:0.1:0.9, sigmas=[0.001,0.01,0.1,0.5,1.0,10.0], reps=100,
    mortality_vals=[0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.5, 1.0, 3.0, 5.0, 7.0],
    number_of_combinations=3000
)
    R = S - C
    scenarios = [:ER, :PL, :MOD]
    cols = [:scenario, :conn, :sigma, :rep,
            :res0,:rea0,:sl0,:sig0,
            :res1,:rea1,:sl1,:sig1,
            :res2,:rea2,:sl2,:sig2]

    results = DataFrame(scenario=Symbol[], conn=Float64[], sigma=Float64[], rep=Int[],
                        res0=Float64[], rea0=Float64[], sl0=Float64[], sig0=Float64[],
                        res1=Float64[], rea1=Float64[], sl1=Float64[], sig1=Float64[],
                        res2=Float64[], rea2=Float64[], sl2=Float64[], sig2=Float64[])

    Random.seed!(123)
    cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        scenarios, 1:reps, sigmas, conns, mortality_vals, growth_vals
    ))

    locki = ReentrantLock()
    println("Number of combinations: ", length(combos))

    @threads for (scen, rep, sigma, conn, m_val, g_val) in sample(combos, min(length(combos), number_of_combinations); replace=false)
    # 1) sample A using make_A!!
        local A = zeros(S,S)
        make_A!(A, R, conn, scen; IS=sigma)
        A[diagind(A)] .= 0
        # mask of nonzero links
        mask = A .!= 0

        # sample epsilon but unused here, just placeholder
        local epsilon = clamp.(randn(S,S).*0.1,0,1)

        ########### #TODO ACTIVATE THIS PART IF YOU WANT TO GO BACK TO ENSURING FULL FEASIBLE #############
        # 2) calibrate Xi,K ? get (R_eq,C_eq,xi_cons,K_res)
        # 2) generate *all* feasible (Xi,K) sets for this A,epsilon
        old_epsilon = copy(epsilon)
        thr_sets = generate_feasible_thresholds(A, epsilon, R) #; margins = 1.0)
        if isempty(thr_sets)
            @warn "No feasible thresholds for $conn, $sigma, $scen"
            continue
        end
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside generate_feasible_thresholds!"

        # 2b) now loop over each feasible threshold-set
        local t = thr_sets[1]
        if any(t.C_eq .<= 0) || any(t.R_eq .<= 0)
            @error "Non-positive equilibrium abundances squized in the loop: $t"
            continue
        end
        xi_cons   = t.xi_cons
        K_res     = t.K_res
        R_eq        = t.R_eq
        C_eq        = t.C_eq
        #########################################################################################
        #########################################################################################
        # set up parameters & run the rest of your pipeline exactly as before
        # m_cons = fill(m_val, C)
        m_cons = abs.(rand(Normal(m_val, 0.2), C))
        # d_res  = ones(R)
        # r_res  = fill(g_val, R)
        r_res  = abs.(rand(Normal(g_val, 0.2), R))

        p     = (R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)
        fixed = vcat(R_eq, C_eq)

        # 3) stability & survival
        ok, B0 = survives!(fixed, p; cb=cb)
        # !ok && continue
        if !all(isapprox.(B0, fixed, atol=1e-3))
            # @warn "B0 is not close to fixed point: B0 =  $(B0), and fixed = $(fixed)"
            # continue
        end
        D, M = compute_jacobian(fixed, p)
        J_full = D * M
        !is_locally_stable(J_full) && continue

        # now for each step collect metrics
        row = Any[scen, conn, sigma, rep]
        for step in 0:2
            A_step = copy(A)
            if step==1
                rewire_A!(A_step, A .!= 0, sigma)
            elseif step==2
                make_A!(A_step,R,conn*2,scen; IS=sigma*2)
            end
            
            p_step = (R,C,p[3], p[4], p[5], p[6], epsilon, A_step)
            D_s,M_s = compute_jacobian(fixed,p_step)
            J_s = D_s*M_s

            lock(locki) do
                push!(row, compute_resilience(fixed,p_step),
                    compute_reactivity(fixed,p_step),
                    mean(compute_SL(A_step, vcat(K_res,xi_cons))),
                    sigma_over_min_d(A_step, J_s))
            end
        end
    end

    # filter: keep only rows where any res*>0
    # results = filter(r-> r.res0>0 || r.res1>0 || r.res2>0, results)
    # step_keys = ["0","1","2"]
    # res_cols = Symbol.("res" .* step_keys)
    # results = filter(row -> all(row[c] < 0 for c in res_cols), results)
    # println("subset size: ", nrow(results))

    # plotting
    metrics = [(:res, "Resilience"), (:rea,"Reactivity"), (:sl,"Mean SL"), (:sig,"σ/min(d)")]
    for scen in scenarios
        df = results[results.scenario .== scen, :]
        fig = Figure(; size=(1000,450))
        for (i,(sym,label)) in enumerate(metrics)
            for (j,step) in enumerate((1,2))  # comparing step0 to step1 and step2
                ax = Axis(fig[j, i];
                    title = "$(scen) $label: orig vs step$(step)",
                    xlabel="orig", ylabel="step$(step)")
                x = df[!, Symbol(string(sym,"0"))]
                y = df[!, Symbol(string(sym,step))]
                scatter!(ax, x, y; alpha=0.3)
                # low, high = extrema(x)
                # lines!(ax, [low,high],[low,high]; color=:black, linestyle=:dash)
            end
        end
    display(fig)
    end

    return results
end

# call it
res_wide = test_rewire_wide()
