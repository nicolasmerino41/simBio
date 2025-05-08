using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads

# ────────────────────────────────────────────────────────────────
# 5) Stability & survival checks
# ────────────────────────────────────────────────────────────────────────────────
is_locally_stable(J) = maximum(real(eigvals(J))) < 0.0

function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false,sol.u[end]) : (all(sol.u[end] .> EXTINCTION_THRESHOLD),sol.u[end])
end

function calibrate_with_margin(A, ε, R; margins=[0.05,0.2,0.5,1.0,2.0])
    for deltamargin in margins
      xi_cons, K_res = compute_default_thresholds(A, ε, R; margin=deltamargin)
      eq = calibrate_from_K_xi(xi_cons, K_res, ε, A)
      if eq !== nothing
        R_eq, C_eq = eq
        return R_eq, C_eq, xi_cons, K_res, deltamargin
      end
    end
    return nothing
  end

# ────────────────────────────────────────────────────────────────────────────────
# Main sweep: fix S & C, vary structure, collect full & ladder pers.
# ────────────────────────────────────────────────────────────────────────────────
function ComputingLadder(
    S::Int=30, C::Int=6;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[0.1,1.0],
    skew_vals=[0.01,0.1],
    abund_types=[:Normal,:Log],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[0.1],
    delta_vals=[1.0,3.0],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100
)
    R = S - C
    A = zeros(S,S)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S,EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        conn_vals,IS_vals,skew_vals,abund_types,scenarios,delta_vals
    ))

    for (conn,IS,skew,atype,scen,delta) in 
        sample(combos, min(length(combos), number_of_combinations); replace=false)

        # 1) build A & ε
        A = make_A(A,R,conn,scen).*IS
        ε = clamp.(randn(S,S).*eps_scales[1],0,1)

        # 2) calibrate ξ,K → get (R_eq,C_eq,xi_cons,K_res)
        calib = calibrate_with_margin(A, ε, R)

        calib === nothing && (@warn "Calibration failed for (conn,IS,skew,atype,scen,delta) = ($conn,$IS,$skew,$atype,$scen,$delta)"; continue)
        R_eq, C_eq, xi_cons, K_res, used_margin = calib

        # set m_cons=ξ, d_res=1, r_res=K_res
        m_cons = xi_cons
        d_res  = ones(R)
        r_res  = K_res

        # 3) package p and check stability & survival
        p      = (R, C, m_cons, xi_cons, r_res, d_res, ε, A)
        fixed  = vcat(R_eq, C_eq)
        D, M    = compute_jacobian(fixed,p)
        !is_locally_stable(D*M) && continue
        ok, B0  = survives!(fixed,p; cb=cb)
        !ok && continue

        # 4) full‐model persistence
        rt, os, ire, before_full, aft_full, new_equil, _ = simulate_press_perturbation(
            B0, p, tspan, tpert, delta;
            solver=Tsit5(),
            plot=true,
            show_warnings=true,
            full_or_simple=true,
            cb=cb,
            species_specific_perturbation=false
        )
        rt_full = mean(filter(!isnan, rt))

        # 5) ladder persistence
        pers_steps = Dict(i => NaN for i in 1:16)
        aft_steps  = Dict(i => NaN for i in 1:16)
        rt_steps   = Dict(i => NaN for i in 1:16)
        for step in 1:16
            # build simplified A and ε
            A_s, ε_s = transform_for_ladder_step(step, A, ε)

            # original equilibrium abundances
            R_eq_full, C_eq_full = R_eq, C_eq

            # 5a) Recompute xi_hat
            xi_hat = zeros(C)
            for k in 1:C
                i = R + k
                # 1) feeding gains (A>0)
                gain = sum( ε_s[i,j]*A_s[i,j]*R_eq_full[j]
                            for j in 1:R if A_s[i,j] > 0; init=0.0 )
                # 2) predation losses: A_s[i,j]<0, but we need “+(-A)B”
                loss = sum( -A_s[i,j]*C_eq_full[j-R]
                            for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                # consumer eq: xi = B_i – gain – loss
                xi_hat[k] = C_eq_full[k] - gain - loss
            end

            # 5b) Recompute K_hat
            K_hat = zeros(R)
            for i in 1:R
                # resource eq uses A[i,j] (j consumer) directly:
                drain = sum( A_s[i,j]*C_eq_full[j-R]
                                for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                # K_i = B_i + ∑ A[i,j] B_j
                K_hat[i] = R_eq_full[i] + drain
            end

            # 5c) Solve for new equilibrium
            eq = calibrate_from_K_xi(xi_hat, K_hat, ε_s, A_s)
            if eq === nothing
                error("Ladder step $step infeasible: xi_hat=$(xi_hat), K_hat=$(K_hat), step=$step")
            end
            R_eq_s, C_eq_s = eq

            # 5d) simulate simplified model
            p_s = (R, C, xi_hat, xi_hat, K_hat, ones(R), ε_s, A_s)
            rt2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                B0, p_s, tspan, tpert, delta;
                solver=Tsit5(),
                plot=false,
                cb=cb,
                full_or_simple=false
            )
            pers_steps[step] = before_s
            aft_steps[step]  = after_s
            rt_steps[step]   = mean(filter(!isnan, rt2))
        end

        step_pairs = collect(Iterators.flatten(
            ([
                Symbol("pers_step_$i") => pers_steps[i],
                Symbol("aft_step_$i") => aft_steps[i],
                Symbol("rt_step_$i") => rt_steps[i]
            ] for i in 1:16)
        ))

        rec = (
            conn=conn, IS=IS, skew=skew, scen=scen, δ=δ,
            pers_full=before_full, aft_full=aft_full, rt_full=rt_full,
            step_pairs...,  # Properly flattened pairs
            p_final = p,
            R_eq = R_eq,
            C_eq = C_eq
            )

        lock(locki) do
            push!(results, rec)
        end
    end

    return DataFrame(results)
end


# ────────────────────────────────────────────────────────────────────────────────
# 9) Run it
# ────────────────────────────────────────────────────────────────────────────────
A = ComputingLadder(
    2, 1;
    conn_vals=0.05:0.05:0.3,
    IS_vals=[1.0],
    skew_vals=[0.01,0.1,0.5],
    scenarios=[:ER,:PL,:MOD],
    delta_vals=[1.0,3.0],
    abund_types=[:Normal,:Log],
    eps_scales=[0.1],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100
)
