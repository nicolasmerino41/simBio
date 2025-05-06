include("Ladder4.1_drago1.jl")
using DifferentialEquations, Random, Statistics, DataFrames, CSV
import Base.Threads: @threads, nthreads, threadid
include("ExploringFeasibleSpace.jl")
# include("Collectivity.jl")
cb_no_trigger30 = build_callbacks(30, EXTINCTION_THRESHOLD)
cb_no_trigger40 = build_callbacks(40, EXTINCTION_THRESHOLD)
cb_no_trigger50 = build_callbacks(50, EXTINCTION_THRESHOLD)

# ----------------------------------------------------------------------------
# 1) Check local stability: max real part of Jacobian < 0
# ----------------------------------------------------------------------------
is_locally_stable(J) = maximum(real(eigvals(J))) < 0

# ----------------------------------------------------------------------------
# 2) Simulate unperturbed dynamics to check 100% persistence
# ----------------------------------------------------------------------------
function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    if sol.t[end] < tspan[2] return false, sol.u[end] end
    Bf = sol.u[end]
    return all(Bf .> EXTINCTION_THRESHOLD), Bf
end

function make_A(
    A::AbstractMatrix{<:Real},
    R::Int,
    conn::Float64,
    scenario::Symbol;
    pareto_exponent::Float64 = 1.75,
    mod_gamma::Float64       = 5.0,
    B_term::Bool             = false
)
    S = size(A,1)
    C = S - R
    fill!(A, 0.0)

    # define the allowed prey‐pool for each consumer i
    function prey_indices(i)
        base = 1:R
        if B_term
            # include other consumers but not self
            return vcat(base, setdiff(R+1:S, i))
        else
            return base
        end
    end

    if scenario == :ER
        # Erdős–Rényi: each candidate link with probability conn
        for i in (R+1):S
            for j in prey_indices(i)
                if rand() < conn
                    A[i,j] = abs(rand(Normal()))
                    A[j,i] = -abs(rand(Normal()))
                end
            end
        end

    elseif scenario == :PL
        # Power‐law consumer out‐degree
        x_m = 1.0
        raw_degrees = rand(Pareto(x_m, pareto_exponent), C)
        ks = clamp.(floor.(Int, raw_degrees), 1, length(prey_indices(R+1)))

        for (idx, k) in enumerate(ks)
            ci = R + idx
            pool = prey_indices(ci)
            chosen = sample(pool, min(k, length(pool)); replace=false)
            for j in chosen
                A[ci, j] = abs(rand(Normal()))
                A[j, ci] = -abs(rand(Normal()))
            end
        end

    elseif scenario == :MOD
        # Bipartite‐modular plus optional B–B links
        halfR = fld(R,2)
        halfC = fld(C,2)
        res_block1  = 1:halfR
        res_block2  = (halfR+1):R
        cons_block1 = (R+1):(R+halfC)
        cons_block2 = (R+halfC+1):S

        for i in (R+1):S
            for j in prey_indices(i)
                # decide module‐match only for resource‐blocks;
                # for consumer–consumer we just use base conn
                if j <= R
                    same = (i in cons_block1 && j in res_block1) ||
                           (i in cons_block2 && j in res_block2)
                    p_link = same ? conn*mod_gamma : conn/mod_gamma
                else
                    # consumer–consumer: uniform random links
                    p_link = conn
                end
                if rand() < clamp(p_link,0,1)
                    A[i,j] = abs(rand(Normal()))
                    A[j,i] = -abs(rand(Normal()))
                end
            end
        end

    else
        error("Unknown scenario: $scenario")
    end

    return A
end

function persistence_sweep(;
    S=[30, 40, 50],
    C_vals=[3,6,9,12],
    conn_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    IS_vals=[0.1,1.0],
    skew_vals=[0.01,0.1,0.5],
    abundance_types=[:Normal,:Log],
    scenarios=[:ER,:PL,:MOD],
    # leave pex and modg fixed at 1.0
    eps_scales=[0.1],
    d_vals=[0.1, 1.0],
    m_vals=[0.1, 0.3],
    abundance_mean=10.0,
    number_of_combinations=1000,
    tspan=(0.,500.0),
    t_perturb=250.0,
    delta_vals=[1.0, 3.0],
    max_calib=10,
)
    locki = ReentrantLock()
    results = NamedTuple[]

    # full cartesian product
    combos = collect(Iterators.product(
        S, C_vals, conn_vals, IS_vals, skew_vals, abundance_types, scenarios, delta_vals, d_vals, m_vals
    ))
    # global cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    @threads for (S, C, conn, IS, skew, abund_type, scenario, delta, d, m) in combos[sample(1:length(combos), min(length(combos), number_of_combinations), replace=false)]
        R = S - C
        if S == 30
            cb = cb_no_trigger30
        elseif S == 40
            cb = cb_no_trigger40
        else
            cb = cb_no_trigger50
        end
        # 1) try calibrating
        xi_cons, r_res = fill(NaN,C), fill(NaN,R)
        p_final = nothing
        R_eq, C_eq = nothing, nothing
        tries = 0

        while any(isnan, xi_cons) && tries < max_calib
            tries += 1

            # build A and eps
            A = make_A(zeros(S,S), R, conn, scenario;
                        pareto_exponent=1.0, mod_gamma=1.0) .* IS
            eps_scale = rand([eps_scales][1])
            eps = clamp.(rand(Normal(eps_scale,eps_scale), S, S), 0, 1)

            # draw equilibrium
            if abund_type == :Log
            R_eq = abs.(rand(LogNormal(log(abundance_mean)-abundance_mean^2/2,
                                        abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*skew)-
                                        (abundance_mean*skew)^2/2,
                                        abundance_mean*skew), C))
            else
            R_eq = abs.(rand(Normal(abundance_mean, abundance_mean*skew), R))
            C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
            end

            # calibrate
            xi_cons, r_res = calibrate_params(
            R_eq, C_eq,
            (R,C, fill(m,C), fill(d,R), eps, A);
            xi_threshold=0.7, constraints=true
            )

            if !any(isnan, xi_cons)
            p_final = (R, C, fill(m,C), xi_cons, r_res, fill(d,R), eps, A)
            break
            end
        end

        isnothing(p_final) && continue

        # 2) local stability + survival
        fixed = vcat(R_eq, C_eq)
        D,M = compute_jacobian(fixed,p_final)
        maximum(real, eigvals(D*M)) >= 0 && continue
        ok, B_eq = survives!(fixed,p_final; tspan=tspan, cb=cb)
        # !ok && continue

        # 3) full‐model persistence
        _,_,_, pers_full, new_equil, _ = simulate_press_perturbation(
            B_eq, p_final, tspan, t_perturb, delta;
            solver=Tsit5(), cb=cb)
        aft_full = sum(new_equil .> EXTINCTION_THRESHOLD)/S

        # 4) ladder‐step persistence
        pers_steps = Dict(i => NaN for i in 1:16)
        aft_steps = Dict(i => NaN for i in 1:16)
        for step in 1:16
            A_s, eps_s = transform_for_ladder_step(step, p_final[8], p_final[7])
            p_s = (R,C, p_final[3], p_final[4], p_final[5], p_final[6], eps_s, A_s)
            _,_,_, pstep, new_equil, _ = simulate_press_perturbation(
            B_eq, p_s, tspan, t_perturb, delta;
            solver=Tsit5(), cb=cb, full_or_simple=false)
            pers_steps[step] = pstep
            aft = sum(new_equil .> EXTINCTION_THRESHOLD)/S
            aft_steps[step] = aft
        end
        step_pairs = collect(Iterators.flatten(
            ([
                Symbol("pers_step_$i") => pers_steps[i],
                Symbol("aft_step_$i") => aft_steps[i]
            ] for i in 1:16)
        ))

        record = (
            S = S, 
            C = C, 
            conn = conn, 
            IS = IS, 
            skew = skew,
            abundance = abund_type, 
            scenario = scenario,
            delta = delta, 
            pers_full = pers_full, 
            aft_full = aft_full,
            d = d, 
            m = m,
            step_pairs...,  # Properly flattened pairs
            p_final = p_final,
            R_eq = R_eq,
            C_eq = C_eq
)

        lock(locki) do
            push!(results, record)
        end
    end

    DataFrame(results)
end

# run the sweep
dfp = persistence_sweep(;
    S=[30,40,50],
    C_vals=[3,6,9,12],
    conn_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    IS_vals=[0.1,1.0, 2.0],
    skew_vals=[0.01,0.1,0.5, 1.0],
    abundance_types=[:Normal,:Log],
    scenarios=[:ER,:PL,:MOD],
    # leave pex and modg fixed at 1.0
    eps_scales=[0.1],
    d_vals=[0.1, 1.0],
    m_vals=[0.1, 0.3],
    abundance_mean=10.0,
    number_of_combinations=100,
    tspan=(0.,500.0),
    t_perturb=250.0,
    delta_vals=[1.0, 3.0, 5.0, 8.0],
    max_calib=10
)
# braoder_dfp = CSV.write("ThePaper/Ladder/Outputs/broader_persistence_sweep.csv", dfp)
serialize("persistence_sweep_304050_withAFT.jls", dfp)