include("Ladder4.1_allo.jl")       # your short_transform_for_ladder_step etc.
@info "ladder4.1_allo.jl loaded."
using EcologicalNetworksDynamics
@info "END loaded."
using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames
@info "Second line loaded"
import Base.Threads: @threads

@info "ShortAlloLadder.jl loaded."
# -----------------------
# 1) Builder from END
# -----------------------
function build_from_scenario(S::Int, C::Int, conn::Real, IS_scale::Real)
    fw = Foodweb(:cascade; S=S, C=conn, reject_if_disconnected=false)
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
    steps=3
)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(conn_vals, IS_vals, deltas, 1:iterations))
    @info "Running short ladder with $(length(combos)) combinations but only $(min(length(combos), number_of_combinations)) will be sampled."
    @threads for (conn, IS_scale, delta, ite) in sample(combos, min(length(combos), number_of_combinations); replace=false)
        # @info "Running iteration $(ite) with conn=$(conn), IS_scale=$(IS_scale), delta=$(delta)"
        # build
        local m, p, attack, epsilon = build_from_scenario(S, C, conn, IS_scale)
        p = (S = p.S, r = p.r, K = p.K, x = p.x, d = p.d, attack = p.attack, epsilon = epsilon, B0_ref = p.B0_ref, i0 = p.i0, h = p.h)
        # full-model simulation
        u0 = rand(S)
        prob_full = ODEProblem(allo_ode!, u0, (0.,500.), p)
        sol_full = solve(prob_full, Tsit5(); callback=cb, abstol=1e-6, reltol=1e-6)
        Bstar = sol_full.u[end]

        # stability & survival
        J_full = compute_jacobian_allo(Bstar, p)
        J_diff_full = norm(J_full - J_full)
        J_full_norm = norm(J_full)
        if !is_locally_stable(J_full)
            continue
        end
        # println("breakpoint 1")
        # surv, Bchk = survives!(Bstar, p; cb=cb)
        # if !surv
        #     continue
        # end

        # metrics full
        res_f = compute_resilience(Bstar, p)
        # println("breaokpoint 1.5")
        rea_f = compute_reactivity(Bstar, p)
        # println("breakpoint 1.75")
        Rm_f  = median_return_rate(J_full, Bstar; t=1.0, n=20)
        # println("breakpoint 2")
        # - full model press -
        rt_press_full, before_full, after_full, _ = simulate_press_perturbation_allo!(Bstar, p, (0.,500.), 250.0, delta; cb=cb)
        rt_pulse_full, _, after_pulse_full, _ = simulate_pulse_perturbation_allo!(Bstar, p, (0.,500.), 250.0, delta; cb=cb)
        # println("breakpoint 3")
        # - simplified model press at each step -
        before_S = Dict(i => NaN for i in 1:steps)
        after_S  = Dict(i => NaN for i in 1:steps)
        after_pulse_S = Dict(i => NaN for i in 1:steps)
        
        rea_S  = Dict(i=>NaN for i in 1:steps)
        res_S  = Dict(i=>NaN for i in 1:steps)
        st_S    = Dict(i=>NaN for i in 1:steps)
        Rm_S    = Dict(i=>NaN for i in 1:steps)
        rt_press_S    = Dict(i=>NaN for i in 1:steps)
        rt_pulse_S    = Dict(i=>NaN for i in 1:steps)
        J_diff_S = Dict(i=>NaN for i in 1:steps)
        
        @info "Running ladder"
        for step in 1:steps
            local A_s, epsilon_s = short_transform_for_ladder_step(step, copy(attack), copy(epsilon))
            p_s = (
            S=p.S, r=p.r, K=p.K, x=p.x, d=p.d,
            attack=A_s,
            epsilon=epsilon_s,
            B0_ref=p.B0_ref, i0=p.i0, h=p.h
            )
            # println("breakpoint 4")
            prob_s = ODEProblem(allo_ode!, Bstar, (0.,500.), p_s)
            sol_s = solve(prob_s, Tsit5(); callback=cb, abstol=1e-6, reltol=1e-6)
            B_s = sol_s.u[end]
            # println("breakpoint 5")

            rt_press_s, before_s, after_s, B_after_s =
                simulate_press_perturbation_allo!(Bstar, p_s, (0.,500.), 250.0, delta; cb=cb)
            rt_pulse_s, before_pulse_s, after_pulse_s, B_after_pulse_s =
                simulate_pulse_perturbation_allo!(Bstar, p_s, (0.,500.), 250.0, delta; cb=cb)
            before_S[step] = before_s
            after_S[step]  = after_s
            after_pulse_S[step] = after_pulse_s
            # println("breakpoint 6")
            rt_press_S[step] = mean(filter(!isnan, rt_press_s))

            J_s = compute_jacobian_allo(B_s, p_s)
            # println("breakpoint 7")
            st_S[step]  = is_locally_stable(J_s)
            res_S[step] = compute_resilience(B_s, p_s)
            rea_S[step] = compute_reactivity(B_s, p_s)
            Rm_S[step]  = median_return_rate(J_s, B_s; t=1.0, n=20)
            J_diff_S[step] = norm(J_s - J_full)
        end

        # flatten the step‐dicts into pairs
        step_pairs = collect(Iterators.flatten(
            ([ 
                Symbol("before_S$i") => before_S[i],
                Symbol("after_S$i")  => after_S[i],
                Symbol("after_pulse_S$i") => after_pulse_S[i],
                Symbol("rt_press_S$i")     => rt_press_S[i],
                Symbol("rt_pulse_S$i")     => rt_pulse_S[i],
                Symbol("resilience_S$i") => res_S[i],
                Symbol("reactivity_S$i") => rea_S[i],
                Symbol("Rmed_S$i")  => Rm_S[i],
                Symbol("is_stable_S$i") => st_S[i],
                Symbol("J_diff_S$i") => J_diff_S[i],
            ] for i in 1:steps)
        ))

        rec = merge(
            (
            conn=conn, IS=IS_scale, delta=delta, J_diff_full=J_diff_full, J_full_norm=J_full_norm,
            before_full=before_full,
            after_full=after_full, after_pulse_full=after_pulse_full,
            rt_press_full=mean(filter(!isnan, rt_press_full)),
            rt_pulse_full=mean(filter(!isnan, rt_pulse_full)),
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
    deltas=[0.1, 0.01],
    iterations=200,
    number_of_combinations=48*16,
    steps=2
)

desired = [
  :conn, :IS, :delta, :J_diff_full, :J_full_norm,
  :before_full, :after_full, :after_pulse_full, :rt_press_full, rt_pulse_full, :resilience_full, :reactivity_full, :Rmed_full, 
  :before_S1, :after_S1, :after_pulse_S1, :rt_press_S1, rt_pulse_S1, :resilience_S1, :reactivity_S1, :Rmed_S1,
  :before_S2, :after_S2, :after_pulse_S2, :rt_press_S2, rt_pulse_S2, :resilience_S2, :reactivity_S2, :Rmed_S2,
#   :before_S3, :after_S3, :after_pulse_S3, :rt_press_S3, rt_pulse_S3, :resilience_S3, :reactivity_S3, :Rmed_S3,
]

A = T[:, desired]

serialize("ShortAlloLadderCascade.jls", A)