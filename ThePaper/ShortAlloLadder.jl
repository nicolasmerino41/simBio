include("Ladder4.1_allo.jl")       # your short_transform_for_ladder_step etc.
@info "ShortAlloLadder.jl loaded."

using EcologicalNetworksDynamics
using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames
import Base.Threads: @threads
# -----------------------
# 1) Builder from END
# -----------------------
function build_from_scenario(S::Int, C::Int, conn::Real, IS_scale::Real)
    fw = Foodweb(:niche; S=S, C=conn, reject_if_disconnected=false)
    m  = default_model(fw)

    # pull out everything
    r       = m.r
    K       = m.K
    x       = m.x
    d       = m.d
    y       = m.y
    epsilon       = m.efficiency
    A_bool  = m.trophic.matrix
    w       = m.w
    B0_ref  = m.half_saturation_density
    i0      = m.intraspecific_interference
    h       = m.hill_exponent

    # build continuous‐attack matrix
    attack = zeros(Float64, S, S)
    @inbounds for i in 1:S, j in 1:S
        attack[i,j] = A_bool[i,j] ? x[i]*y[i]*w[i,j] : 0
    end

    # **here is the key**: a named tuple, so p.S, p.r, etc. work
    p = (
        S       = S,
        r       = r,
        K       = K,
        x       = x,
        d       = d,
        attack  = attack,
        epsilon = epsilon,
        B0_ref  = B0_ref,
        i0      = i0,
        h       = h
    )
    return m, p, attack, epsilon
end


# -----------------------
# 5) Short ladder pipeline
# -----------------------
function short_ComputingLadder_END(
    S::Int=50, C::Int=20;
    conn_vals=[0.05, 0.1, 0.2],
    IS_vals=[0.1, 1.0],
    deltas=[0.1, 0.5, 1.0],
    number_of_combinations=100,
    iterations=100,
)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(conn_vals, IS_vals, deltas, 1:iterations))
    @info "Running short ladder with $(length(combos)) combinations but only $(min(length(combos), number_of_combinations)) will be sampled."
    @threads for (conn, IS_scale, delta, ite) in sample(combos, min(length(combos), number_of_combinations); replace=false)
        # build
        local m, p, attack, epsilon = build_from_scenario(S, C, conn, IS_scale)
        p = (S = p.S, r = p.r, K = p.K, x = p.x, d = p.d, attack = p.attack, epsilon = epsilon, B0_ref = p.B0_ref, i0 = p.i0, h = p.h)
        # full-model simulation
        u0 = rand(S)
        prob_full = ODEProblem(allo_ode!, u0, (0.,500.), p)
        sol_full = solve(prob_full, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
        Bstar = sol_full.u[end]

        # stability & survival
        Jfull = compute_jacobian_allo(Bstar, p)
        if !is_locally_stable(Jfull)
            continue
        end
        # surv, Bchk = survives!(Bstar, p; cb=cb)
        # if !surv
        #     continue
        # end

        # metrics full
        res_f = compute_resilience(Bstar, p)
        rea_f = compute_reactivity(Bstar, p)
        Rm_f  = median_return_rate(Jfull, Bstar; t=1.0, n=200)
        
        # - full model press -
        rt_full, before_full, after_full, _ = simulate_press_perturbation_allo!(Bstar, p, (0.,500.), 250.0, delta; cb=cb)

        # - simplified model press at each step -
        before_S = Dict(i => NaN for i in 1:3)
        after_S  = Dict(i => NaN for i in 1:3)
        
        rea_S  = Dict(i=>NaN for i in 1:3)
        res_S  = Dict(i=>NaN for i in 1:3)
        st_S    = Dict(i=>NaN for i in 1:3)
        Rm_S    = Dict(i=>NaN for i in 1:3)
        rt_S    = Dict(i=>NaN for i in 1:3)
        
        @info "Running ladder"
        for step in 1:3
            local A_s, epsilon_s = short_transform_for_ladder_step(step, copy(attack), copy(epsilon))
            p_s = (
            S=p.S, r=p.r, K=p.K, x=p.x, d=p.d,
            attack=A_s,
            epsilon=epsilon_s,
            B0_ref=p.B0_ref, i0=p.i0, h=p.h
            )
            prob_s = ODEProblem(allo_ode!, Bstar, (0.,500.), p_s)
            sol_s = solve(prob_s, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
            B_s = sol_s.u[end]

            rt_s, before_s, after_s, B_after_s =
                simulate_press_perturbation_allo!(Bstar, p_s, (0.,500.), 250.0, delta; cb=cb)
            before_S[step] = before_s
            after_S[step]  = after_s

            rt_S[step] = mean(filter(!isnan, rt_s))

            J_s = compute_jacobian_allo(B_s, p_s)
            st_S[step]  = is_locally_stable(J_s)
            res_S[step] = compute_resilience(B_s, p_s)
            rea_S[step] = compute_reactivity(B_s, p_s)
            Rm_S[step]  = median_return_rate(J_s, B_s; t=1.0, n=20)
        end

        # flatten the step‐dicts into pairs
        step_pairs = collect(Iterators.flatten(
            ([ 
                Symbol("before_S$i") => before_S[i],
                Symbol("after_S$i")  => after_S[i],
                Symbol("rt_S$i")     => rt_S[i],
                Symbol("resilience_S$i") => res_S[i],
                Symbol("reactivity_S$i") => rea_S[i],
                Symbol("Rmed_S$i")  => Rm_S[i],
                Symbol("is_stable_S$i") => st_S[i],
            ] for i in 1:3)
        ))

        rec = merge(
            (
            conn=conn, IS=IS_scale, delta=delta,
            before_full=before_full,
            after_full=after_full,
            rt_full=mean(filter(!isnan, rt_full)),
            resilience_full=res_f,
            reactivity_full=rea_f,
            Rmed_full=Rm_f,
            ),
            Dict(step_pairs...)
        )

        lock(locki) do
            push!(results, rec)
        end
    end

    return DataFrame(results)
end

# -----------------------
# 3) Run the pipeline
# -----------------------
T = short_ComputingLadder_END(
    50, 20;
    conn_vals=0.01:0.025:0.5,
    IS_vals=[1.0],
    deltas=[0.1],
    iterations=200,
    number_of_combinations=42
)

desired = [
  :conn, :IS, :delta,
  :before_full, :after_full, :rt_full, :resilience_full, :reactivity_full, :Rmed_full,
  :before_S1, :after_S1, :rt_S1, :resilience_S1, :reactivity_S1, :Rmed_S1,
  :before_S2, :after_S2, :rt_S2, :resilience_S2, :reactivity_S2, :Rmed_S2,
  :before_S3, :after_S3, :rt_S3, :resilience_S3, :reactivity_S3, :Rmed_S3
]

A = T[:, desired]

serialize("ShortAlloLadder.jls", A)