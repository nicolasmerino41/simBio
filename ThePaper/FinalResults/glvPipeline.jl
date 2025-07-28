using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")
# --------------------------------------------------------------------------------
# 1) gLV dynamics
# --------------------------------------------------------------------------------
# du[i] = u[i] * (K[i] - u[i] - (A*u)[i])
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for i in eachindex(u)
        du[i] = u[i] * (K[i] - u[i] + Au[i])
    end
end

# --------------------------------------------------------------------------------
# 2) Compute SL = B ./ K
# --------------------------------------------------------------------------------
function compute_SL(A::Matrix{Float64}, K::Vector{Float64})
    B = (I + A) 
            K    # equilibrium: (I + A) * B = K
    return B ./ K
end

# --------------------------------------------------------------------------------
# 3) Extract sigma/min(d) for Jacobian J
# --------------------------------------------------------------------------------
function sigma_over_min_d(A, J)
    d = -diag(J)
    if isempty(d)
        return NaN
    end
    min_d = minimum(d)
    offs = [A[i,j] for i in 1:size(A,1), j in 1:size(A,1) if i != j]
    if isempty(offs)
        return NaN
    end
    sigma = std(offs)
    return sigma / min_d
end

# --------------------------------------------------------------------------------
# 4) rewire in place
# --------------------------------------------------------------------------------
function rewire_A!(A, mask, sigma)
    A[mask] .= randn(sum(mask)) .* sigma
end

# --------------------------------------------------------------------------------
# 5) Equilibrium and feasibility for gLV
# --------------------------------------------------------------------------------
function calibrate_from_K_A(K::Vector{<:Real}, A::AbstractMatrix)
    # Solve (I + A) * u = K
    u = (I - A) \ K
    return u
end

function generate_feasible_thresholds(A::AbstractMatrix; margins=[1.0])
    S = size(A,1)
    out = NamedTuple[]
    for marg in margins
        # propose K = abs(randn(S)) * marg
        K = abs.(rand(Normal(2.0, 1.0), S) .* marg)
        K[31:50] .= 0.01
        u_eq = try
            calibrate_from_K_A(K, A)
        catch
            continue
        end
        if any(!isfinite, u_eq) || any(u_eq .<= 0)
            continue
        end
        # stability check
        J = compute_jacobian_glv(u_eq, (K, A))
        # J = D * M
        if is_locally_stable(J)
            push!(out, (K=K, margin=marg, u_eq=copy(u_eq)))
        end
    end
    return out
end

# --------------------------------------------------------------------------------
# 6) Stability check
# --------------------------------------------------------------------------------
function is_locally_stable(J::AbstractMatrix)
    if any(!isfinite, J)
        return false
    end
    mu = eigvals(J)
    maximum(real.(mu)) <= 0
end

# --------------------------------------------------------------------------------
# 7) Survival simulation
# --------------------------------------------------------------------------------
function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(gLV_rhs!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false, sol.u[end]) : (all(sol.u[end] .> 1e-6), sol.u[end])
end

# --------------------------------------------------------------------------------
# 8) Jacobian, resilience, reactivity
# --------------------------------------------------------------------------------
function compute_jacobian(u, p)
    K, A = p
    S = length(u)
    # D = diag(u)
    D = Diagonal(u)
    # Mstar = -I - A
    Mstar = -I(S) - A
    return D, Mstar
end

function compute_jacobian_glv(Bstar::AbstractVector, p)
    # g(u) returns f(u)the time‐derivative at state u
    g(u) = begin
        du = similar(u)
        gLV_rhs!(du, u, p, 0.0)
        return du
    end

    # allocate output
    S = length(Bstar)
    J = Matrix{Float64}(undef, S, S)

    # fill J in‐place
    ForwardDiff.jacobian!(J, g, Bstar)

    return J
end

function compute_resilience(u, p; extinct_species = false)
    if !extinct_species 
        J = compute_jacobian_glv(u, p)
        ev = eigvals(J)
        return maximum(real.(ev))
    else
        extant = findall(bi -> bi > EXTINCTION_THRESHOLD, u)
        J = compute_jacobian_glv(u, p)
        Jsub = J[extant, extant]
        ev = eigvals(Jsub)
        return maximum(real.(ev))
    end
end

# Reactivity: maximum eigenvalue of the symmetric part of the Jacobian.
function compute_reactivity(B, p; extinct_species = false)
    if !extinct_species 
        J = compute_jacobian_glv(B, p)
        J_sym = (J + J') / 2
        ev_sym = eigvals(J_sym)
        return maximum(real.(ev_sym))
    else
        extant = findall(bi -> bi > EXTINCTION_THRESHOLD, B)
        J = compute_jacobian_glv(B, p)
        Jsub = J[extant, extant]
        J_sym = (Jsub + Jsub') / 2
        ev_sym = eigvals(J_sym)
        return maximum(real.(ev_sym))
    end
end

# --------------------------------------------------------------------------------
# 9) Compute collectivity
# --------------------------------------------------------------------------------
function compute_collectivity(A::AbstractMatrix)
    vals = eigvals(A)
    return maximum(abs, vals)
end

# --------------------------------------------------------------------------------
# 10) Press and pulse perturbation functions remain, replacing trophic_ode! by gLV_rhs!
#     (update p and ODEProblem accordingly)
# --------------------------------------------------------------------------------
##########################
# 10) Press perturbation
##########################
function simulate_press_perturbation_glv(
    u0, p, tspan, t_perturb, delta;
    solver=Tsit5(), plot=false, cb=nothing
)
    # Phase 1: pre-perturb
    prob1 = ODEProblem(gLV_rhs!, u0, (tspan[1], t_perturb), p)
    sol1  = solve(prob1, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    pre_state = sol1.u[end]
    before_persistence = count(x -> x > 1e-6, pre_state) / length(pre_state)
    
    # Phase 2: perturb K by (1 - delta)
    K, A = p
    K_press = vcat(K[1:30] .* (1 .- delta), K[31:end])
    p_press = (K_press, A)
    prob2 = ODEProblem(gLV_rhs!, pre_state, (t_perturb, tspan[2]), p_press)
    sol2  = solve(prob2, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    new_equil = sol2.u[end]
    after_persistence = count(x -> x > 1e-6, new_equil) / length(new_equil)

    # Return times
    n = length(new_equil)
    return_times = fill(NaN, n)
    for i in 1:n
        target = new_equil[i]
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target) / (abs(target) + 1e-8) < 0.10
                return_times[i] = t - t_perturb
                break
            end
        end
    end

    return return_times, before_persistence, after_persistence, new_equil
end

##########################
# 11) Pulse perturbation
##########################
function simulate_pulse_perturbation_glv(
    u0, p, tspan, t_pulse, delta;
    solver=Tsit5(), plot=false, cb=nothing
)
    # Phase 1: pre-pulse
    prob1 = ODEProblem(gLV_rhs!, u0, (tspan[1], t_pulse), p)
    sol1  = solve(prob1, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    pre_state = sol1.u[end]
    before_persistence = count(x -> x > 1e-6, pre_state) / length(pre_state)

    # Apply pulse: u -> u * (1 + delta)
    pulsed = pre_state .* (1 .- delta)
    prob2 = ODEProblem(gLV_rhs!, pulsed, (t_pulse, tspan[2]), p)
    sol2  = solve(prob2, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    eq_state = sol2.u[end]
    after_persistence = count(x -> x > 1e-6, eq_state) / length(eq_state)

    # Return times
    n = length(eq_state)
    return_times = fill(NaN, n)
    for i in 1:n
        target = eq_state[i]
        for (t, u) in zip(sol2.t, sol2.u)
            if abs(u[i] - target) / (abs(target) + 1e-8) < 0.10
                return_times[i] = t - t_pulse
                break
            end
        end
    end

    return return_times, before_persistence, after_persistence, eq_state
end

"""
    analytical_median_return_rate(J; t=0.01)

Computes the median across species of your fully-analytical
species-level return-rates, i.e. the direct analogue of
median_return_rate but without Monte Carlo.
"""
function analytical_median_return_rate(J::AbstractMatrix; t::Real=0.01)
    rates = analytical_species_return_rates(J; t=t)
    return median(rates)
end

function checking_recalculating_demography_glv(
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
        thr_sets = generate_feasible_thresholds(A; margins=[marg])
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
            simulate_press_perturbation_glv(u0, (K,A), tspan, tpert, delta; cb=cb)
        rt_press_full   = mean(skipmissing(rt_press_vec))

        # 5b) pulse perturbation
        rt_pulse_vec, before_pulse, after_pulse, _ =
            simulate_pulse_perturbation_glv(u0, (K,A), tspan, tpert, delta; cb=cb)
        rt_pulse_full   = mean(skipmissing(rt_pulse_vec))

        rmed_full = analytical_median_return_rate(J; t=1.0)

        # 6) ladder persistence (5 simplified steps)
        after_press_S   = Dict(i=>NaN for i in 1:5)
        after_pulse_S   = Dict(i=>NaN for i in 1:5)
        rt_press_S     = Dict(i=>NaN for i in 1:5)
        rt_pulse_S     = Dict(i=>NaN for i in 1:5)
        S_S                   = Dict(i=>NaN for i in 1:5)
        collectivity_S        = Dict(i=>NaN for i in 1:5)
        resilience_S          = Dict(i=>NaN for i in 1:5)
        reactivity_S          = Dict(i=>NaN for i in 1:5)
        sigma_over_min_d_S    = Dict(i=>NaN for i in 1:5)
        SL_S                  = Dict(i=>Float64[] for i in 1:5)
        mean_SL_S             = Dict(i=>NaN for i in 1:5)
        rmed_S                = Dict(i=>NaN for i in 1:5)

        @info "Running ladder steps"
        for step in 1:5
            A_s = copy(A)
            if step == 1
                A_s = make_A(A_s, R, conn, scen; IS=IS)
            elseif step == 2
                new_conn = rand()
                while abs(new_conn - conn) < 0.4
                    new_conn = rand()
                end
                A_s = make_A(A_s, R, new_conn, scen; IS=IS)
            elseif step == 3
                A_s = make_A(A_s, R, conn, scen; IS=IS*10)
            elseif step == 4
                new_conn = rand()
                while abs(new_conn - conn) < 0.4
                    new_conn = rand()
                end
                A_s = make_A(A_s, R, new_conn, scen; IS=IS*10)
            elseif step == 5
                A_s = make_A(A_s, R-5, conn, scen; IS=IS)
            end
            if step == 5
               K[R-5:end] .= 0.0
            end
            # recalibrate equilibrium under simplified A_s
            u_eq_s = try
                calibrate_from_K_A(K, A_s)
            catch
                continue
            end

            ok2, u0_s = survives!(u_eq_s, (K,A_s); cb=cb)
            # if !ok2
            #     continue
            # end

            # record persistence & richness
            S_S[step]                = count(x->x>1e-6, u0_s)
            rt_s_press, _, after_press_s, _=
                simulate_press_perturbation_glv(u0_s, (K,A_s), tspan, tpert, delta; cb=cb)
            after_press_S[step]= after_press_s
            rt_press_S[step]  = mean(skipmissing(rt_s_press))
            rt_s_pulse, _, after_pulse_s, _=
                simulate_pulse_perturbation_glv(u0_s, (K,A_s), tspan, tpert, delta; cb=cb)
            after_pulse_S[step]= after_pulse_s
            rt_pulse_S[step]  = mean(skipmissing(rt_s_pulse))

            # simplified metrics
            collectivity_S[step]     = compute_collectivity(A_s)
            resilience_S[step]       = compute_resilience(u0_s, (K,A_s); extinct_species=false)
            reactivity_S[step]       = compute_reactivity(u0_s, (K,A_s); extinct_species=false)
            J_s_sub                 = compute_jacobian_glv(u0_s, (K,A_s))

            rmed_S[step]             = analytical_median_return_rate(J_s_sub; t=1.0)
            # J_s_sub                 = (D_s*M_s)[findall(x->x>1e-6,u0_s),
            #                                      findall(x->x>1e-6,u0_s)]
            sigma_over_min_d_S[step] = sigma_over_min_d(A_s, J_s_sub)
            # SL_S[step]               = mean(compute_SL(A_s, K))
            SL_S[step]               = diag(J_s_sub) .|> x -> x == 0.0 ? 0.0 : -1 / x
            mean_SL_S[step]          = mean(SL_S[step])
        end

        # flatten ladder dictionaries
        step_pairs = collect(Iterators.flatten(
            ([ 
                Symbol("after_press_S$i") => after_press_S[i],
                Symbol("after_pulse_S$i")       => after_pulse_S[i],
                Symbol("rt_press_S$i")          => rt_press_S[i],
                Symbol("rt_pulse_S$i")          => rt_pulse_S[i],
                Symbol("S_S$i")               => S_S[i],
                Symbol("collectivity_S$i")    => collectivity_S[i],
                Symbol("resilience_S$i")      => resilience_S[i],
                Symbol("reactivity_S$i")      => reactivity_S[i],
                Symbol("sigma_over_min_d_S$i")=> sigma_over_min_d_S[i],
                Symbol("SL_S$i")              => SL_S[i],
                Symbol("mean_SL_S$i")         => mean_SL_S[i],
                Symbol("rmed_S$i")            => rmed_S[i],
             ] for i in 1:5)
        ))

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
            step_pairs...,
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

# Invocation & serialization
R = checking_recalculating_demography_glv(
    50, 20;
    conn_vals=0.01:0.01:1.0,
    scenarios=[:ER, :PL, :MOD],
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    delta_vals=[0.1, 0.9, 1.1, 1.5, 2.0, 3.0, 4.0, 5.0, -1.0, -2.0, -3.0, -4.0, -5.0],
    margins=[1.0, 2.0, 3.0, 4.0, 5.0, 0.01],
    number_of_combinations=10000,
    iterations=5,
    pareto_exponents=[1.0, 1.25, 1.75, 2.0, 3.0, 4.0, 5.0],
    pareto_minimum_degrees=[5.0, 10.0, 15.0, 20.0],
    mod_gammas=[1.0,2.0,3.0,5.0,10.0]
)
serialize("checking_glv_10000ALL.jls", R)

# include("orderColumns.jl")
# R = reorder_columns(R)
R = deserialize("checking_glv_50000ALL.jls")
R = deserialize("checking_10000ALL.jls")