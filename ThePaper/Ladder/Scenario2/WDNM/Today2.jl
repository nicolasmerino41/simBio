using Random, LinearAlgebra, Statistics, DataFrames, CairoMakie, Distributions

# — helper: build Jacobian from A and d
function build_J_from_A(d::Vector{Float64}, A::Matrix{Float64})
    D = diagm(d)
    return D * ( -I(size(A,1)) + A )
end

# — rewire in place
function rewire_A!(A, mask, σ)
    A[mask] .= randn(sum(mask)) .* σ
end

# — make_A from your pipeline
function make_A!(A,R,conn,scenario; IS=1.0, pareto_exponent=1.75, mod_gamma=5.0, B_term=false)
    S = size(A,1); C=S-R
    fill!(A,0)
    if scenario==:ER
        for i in R+1:S, j in 1:R
            if rand()<conn
                A[i,j]=abs(randn()*IS); A[j,i]=-abs(randn()*IS)
            end
        end
    elseif scenario==:PL
        raw=rand(Pareto(1,pareto_exponent),C); ks=clamp.(floor.(Int,raw),1,R)
        for (idx,k) in enumerate(ks)
            i=R+idx
            for j in sample(1:R,min(k,R);replace=false)
                A[i,j]=abs(randn()*IS); A[j,i]=-abs(randn()*IS)
            end
        end
    elseif scenario==:MOD
        halfR,halfC=fld(R,2),fld(C,2)
        res1,res2=1:halfR,(halfR+1):R
        con1,con2=R+1:R+halfC,(R+halfC+1):S
        for i in R+1:S, j in 1:R
            same=(i in con1 && j in res1)||(i in con2 && j in res2)
            p=clamp(conn*(same ? mod_gamma : 1/mod_gamma),0,1)
            if rand()<p
                A[i,j]=abs(randn()*IS); A[j,i]=-abs(randn()*IS)
            end
        end
    else error("Unknown scenario") end
    if B_term
        for i in R+1:S, j in R+1:S
            if i!=j && rand()<conn
                A[i,j]=abs(randn()*IS*0.1); A[j,i]=-abs(randn()*IS*0.1)
            end
        end
    end
    return A
end

# — main pipeline
function test_rewiring_scenarios(; S=50, C=20, conns=0.1:0.1:0.9, sigmas=[0.01, 0.1, 0.5, 1.0, 10.0], reps=100)
    R = S - C
    scenarios = [:ER, :PL, :MOD]
    results = DataFrame(scenario=Symbol[], step=Int[], resilience=Float64[],
                        reactivity=Float64[], sl=Float64[], sigma_ratio=Float64[])
    Random.seed!(123)
    cb = build_callbacks(S,EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        scenarios, 1:reps, sigmas, conns
    ))

    for (scen, rep, sigma, conn) in combos
        # 1) sample A using make_A!!
        A = zeros(S,S)
        make_A!(A, R, conn, scen; IS=sigma)
        A[diagind(A)] .= 0
        # mask of nonzero links
        mask = A .!= 0

        # sample epsilon but unused here, just placeholder
        epsilon = clamp.(randn(S,S)*0.1,0,1)

        ########### #TODO ACTIVATE THIS PART IF YOU WANT TO GO BACK TO ENSURING FULL FEASIBLE #############
        # 2) calibrate Xi,K ? get (R_eq,C_eq,xi_cons,K_res)
        # 2) generate *all* feasible (Xi,K) sets for this A,epsilon
        old_epsilon = copy(epsilon)
        thr_sets = generate_feasible_thresholds(A, epsilon, R) #; margins = 1.0)
        if isempty(thr_sets)
            # @warn "No feasible thresholds for $conn, $IS, $scen, $epsi"
            continue
        end
        @assert all(old_epsilon .== epsilon) "epsilon was mutated inside generate_feasible_thresholds!"

        isempty(thr_sets) && ( @warn "No feasible thresholds for $conn, $IS, $scen, $epsi)" ; continue )

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
        m_cons = abs.(rand(Normal(0.1, 0.2), C))
        # d_res  = ones(R)
        # r_res  = fill(g_val, R)
        r_res  = abs.(rand(Normal(1.0, 0.2), R))

        p     = (R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)
        fixed = vcat(R_eq, C_eq)

        # 3) stability & survival
        ok, B0 = survives!(fixed, p; cb=cb)
        # !ok && continue
        if !all(isapprox.(B0, fixed, atol=1e-3))
            @warn "B0 is not close to fixed point: B0 =  $(B0), and fixed = $(fixed)"
            continue
        end
        D, M = compute_jacobian(fixed, p)
        J_full = D * M
        !is_locally_stable(J_full) && continue

        for step in 0:2
            # prepare A_step exactly as before
            A_step = copy(A)
            if step == 1
                rewire_A!(A_step, mask, sigma)
            elseif step == 2
                new_conn = conn*2
                new_IS   = sigma*2
                mask2    = rand(S,S) .< new_conn
                mask2   .&= .!I(S)
                make_A!(A_step, R, new_conn, scen; IS=new_IS)
                mask = mask2
            end

            # rebuild parameter tuple exactly as in short_ComputingLadder
            p_step = (R, C,
                      m_cons,               # from your outer scope
                      xi_cons,
                      r_res,
                      K_res,
                      epsilon,
                      A_step)

            # compute Jacobian & local equilibrium B0 (already fixed)
            D_step, M_step = compute_jacobian(fixed, p_step)
            J_step = D_step * M_step

            # now measure metrics
            ρ = compute_resilience(fixed, p_step)
            η = compute_reactivity(fixed, p_step)
            sl = mean(compute_SL(A_step, vcat(K_res, xi_cons)))
            σr = sigma_over_min_d(A_step, J_step)

            push!(results, (scen, step, ρ, η, sl, σr))
        end
    end

    

    # — plot with 3 rows (scenarios) × 4 columns (metrics)
    fig = Figure(; size=(1000, 465))
    metrics = [
    (:resilience, "Resilience"),
    (:reactivity, "Reactivity"),
    (:sl,         "Mean SL"),
    (:sigma_ratio,"σ/min(d)")
    ]
    colors = Dict(:ER=>:blue, :PL=>:green, :MOD=>:purple)

    for (r, scen) in enumerate(scenarios)
        sub = results[results.scenario .== scen, :]
        for (c, (col, label)) in enumerate(metrics)
            ax = Axis(fig[r, c];
                title = string(scen, " — ", label),
                xlabel = "Step", ylabel = label,
                xticks = (0:2, ["orig","rewire","rewire+"]),
            )
            for step in 0:2
                vals = sub[sub.step .== step, col]
                scatter!(ax, fill(step, length(vals)), vals;
                        color = colors[scen], alpha = 0.3)
                m = mean(vals)
                lines!(ax, [step-0.2, step+0.2], [m, m];
                    color = colors[scen], linewidth = 2)
            end
        end
    end

    display(fig)

    for scen in scenarios
        # filter once
        df_s = results[results.scenario .== scen, :]

        fig = Figure(resolution=(1000, 300))
        for (c, (col, _)) in enumerate(metrics)
            # orig vs rewire
            ax1 = Axis(
                fig[1, c];
                title   = "$(scen) $(col): orig vs rewire",
                xlabel  = "orig", ylabel="rewire",
                titlesize=10, xticksize=6, yticksize=6,
                xlabelsize=8, ylabelsize=8
            )
            M0 = df_s[df_s.step .== 0, col]
            M1 = df_s[df_s.step .== 1, col]
            scatter!(ax1, M0, M1; alpha=0.3)
            low, high = minimum(M0), maximum(M0)
            lines!(ax1, [low, high], [low, high]; color=:black, linestyle=:dash)

            # orig vs rewire+
            ax2 = Axis(
                fig[2, c];
                title   = "$(scen) $(col): orig vs rewire+",
                xlabel  = "orig", ylabel="rewire+",
                titlesize=10, xticksize=6, yticksize=6,
                xlabelsize=8, ylabelsize=8
            )
            M2 = df_s[df_s.step .== 2, col]
            scatter!(ax2, M0, M2; alpha=0.3)
            lines!(ax2, [low, high], [low, high]; color=:black, linestyle=:dash)
        end

        display(fig)
    end
    
    return results
end

# run
res = test_rewiring_scenarios()
