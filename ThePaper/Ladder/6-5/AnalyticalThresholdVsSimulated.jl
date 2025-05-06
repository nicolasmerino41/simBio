# loosened signature: accept any 8‐element tuple
# ----------------------------------------------------------------------------
# 1) A little helper to perturb p in various ways
# ----------------------------------------------------------------------------
function perturb_tuple(p::Tuple, δ; which::Symbol = :consumers, consumer_idx::Int = 0)
    @assert length(p) == 8 "Expected an 8‐element parameter tuple"
    R, C, m_cons, xi_cons, r_res, d_res, ε, A = p

    xi2 = copy(xi_cons)
    r2  = copy(r_res)

    if which === :consumers
        # bump all consumers
        xi2 .*= 1 .+ δ
    elseif which === :single
        @assert 1 ≤ consumer_idx ≤ C "consumer_idx must be in 1:C"
        xi2[consumer_idx] *= 1 + δ
    elseif which === :resources
        # bump all resources’ carrying capacities
        r2 .*= 1 .+ δ
    elseif which === :all
        xi2 .*= 1 .+ δ
        r2  .*= 1 .+ δ
    else
        error("Unknown perturbation target: $which")
    end

    return (R, C, m_cons, xi2, r2, d_res, ε, A)
end


# ----------------------------------------------------------------------------
# 2) Find the critical press via bisection
# ----------------------------------------------------------------------------
function find_critical_press(
    p, B_eq;
    which::Symbol      = :consumers,
    consumer_idx::Int  = 0,
    δ_lo::Float64      = 0.0,
    δ_hi::Float64      = 1.0,
    pers_target::Float64=1.0,
    tol::Float64       = 1e-3,
    tspan              = (0., 500.),
    t_perturb::Float64 = 250.0,
    cb
)
    R, C = p[1], p[2]

    # helper: run a press of size δ and return persistence
    simulate_pers = δ -> begin
        # properly build the perturbed parameter tuple
        p2 = perturb_tuple(p, δ; which=which, consumer_idx=consumer_idx)
        _, _, _, pers, _, _ = simulate_press_perturbation(
            B_eq, p2, tspan, t_perturb, 0.0;
            full_or_simple=true, cb=cb
        )
        return pers
    end

    # make sure δ_lo is safe
    if simulate_pers(δ_lo) < pers_target
        error("Persistence at δ_lo=$(δ_lo) is already below target")
    end
    # grow δ_hi until we drop below
    while simulate_pers(δ_hi) >= pers_target
        δ_hi *= 2
        if δ_hi > 1e2
            error("Could not bracket critical δ; δ_hi went beyond $(δ_hi)")
        end
    end

    # now bisect
    while δ_hi - δ_lo > tol
        δ_mid = (δ_lo + δ_hi)/2
        simulate_pers(δ_mid) >= pers_target ? (δ_lo = δ_mid) : (δ_hi = δ_mid)
    end

    return δ_lo
end

# 3) Find the critical press for each individual consumer
function find_critical_press_per_consumer(p, B_eq;
                                          δ_lo=0.0, δ_hi=1.0,
                                          pers_target=1.0, tol=1e-3,
                                          tspan=(0.0,500.0), t_perturb=250.0,
                                          cb)
    R, C = p[1], p[2]
    crits = zeros(C)
    for j in 1:C
        # persistence function for consumer j only
        f(δ) = begin
            p2 = perturb_tuple(p, δ; which=:single, consumer_idx=j)
            _,_,_, pers, _, _ = simulate_press_perturbation(
                B_eq, p2, tspan, t_perturb, δ;
                cb=cb, full_or_simple=true)
            return pers
        end

        lo, hi = δ_lo, δ_hi
        @assert f(lo) >= pers_target
        while f(hi) >= pers_target
            hi *= 2
        end

        while hi - lo > tol
            mid = (lo+hi)/2
            if f(mid) >= pers_target
                lo = mid
            else
                hi = mid
            end
        end

        crits[j] = lo
    end
    return crits
end

# === Example usage ===
# (assumes you have `p`, `B_eq`, and `cb` already defined)
δ_consumers = find_critical_press(p, B_eq; which=:consumers, cb=cb_no_trigger36)
δ_resources = find_critical_press(p, B_eq; which=:resources, cb=cb_no_trigger36)
for i in 1:p[2]
    δ_single = find_critical_press(p, B_eq; which=:single, cb=cb_no_trigger36, consumer_idx = i)
    println(δ_single)
end
δ_all = find_critical_press(p, B_eq; which=:all, cb=cb_no_trigger36)

using LinearAlgebra

"""
    analytical_critical_press(A, ε, B_eq; xi_cons=nothing)

Compute the largest relative press δ such that no species goes below zero,
via the analytic sensitivity matrix V = −(I − ε.*A)⁻¹.

If you want to weight each consumer by its equilibrium threshold ξᵢ,
pass xi_cons as a length-C vector; otherwise all presses are equal (–1).
"""
"""
    analytical_critical_press(A, ε, B_eq, xi_cons)

Compute the critical **negative** relative press δ* such that
applying `p = δ* ⋅ press` (with `press = [0…0, 1…1]`) to consumer
thresholds just drives one species to zero.

# Arguments
- `A`      : S×S interaction matrix (losses are negative A[i,j], gains positive A[i,j])
- `ε`      : S×S efficiency matrix (only applied to positive A entries)
- `B_eq`   : length-S equilibrium biomasses
- `xi_cons`: length-C vector of consumer thresholds (only used to size press)

# Returns
- `δ*` : a negative Float64, the largest (most negative) fractional
         press magnitude that still kills nothing; beyond |δ*| the
         first extinction occurs.
"""
"""
    acp_threshold(A, ε, B_eq, xi_cons; thresh=EXTINCTION_THRESHOLD)

Compute the analytic critical press δ* at which any species abundance
would first fall _to_ the extinction threshold `thresh`, rather than zero.

Arguments:
- `A`            : interaction matrix
- `ε`            : scaling matrix (only positive A entries are used)
- `B_eq`         : equilibrium abundances
- `xi_cons`      : consumer thresholds (vector of length C)
- `thresh`       : extinction threshold (default: global EXTINCTION_THRESHOLD)


"""
function acp_threshold(
    A::AbstractMatrix,
    ε::AbstractMatrix,
    B_eq::AbstractVector,
    xi_cons::AbstractVector;
    thresh = EXTINCTION_THRESHOLD
)
    S = size(A,1)
    C = length(xi_cons)
    R = S - C

    # 1) Effective A*: apply ε only to positive entries of A
    eps_eff = copy(ε)
    eps_eff[A .<= 0.0] .= 0.0
    A_star = eps_eff .* A

    # 2) True sensitivity operator V = -inv(J) where J = D * Mstar
    #    but here for simplicity we use V = -(I - A*)^-1
    V = -inv(I(S) .- A_star)

    # 3) Build unit‐press vector: 0 for resources, –1 for each consumer
    press = vcat(zeros(R), ones(C))

    # 4) Analytic per‐unit response ΔB_ana
    ΔB_ana = V * press

    # 5) For each species that declines (ΔB_ana[i]<0),
    #    find δᵢ solving:  B_eq[i] + δ ΔB_ana[i] = thresh
    δs = Float64[]
    for i in 1:S
        if ΔB_ana[i] < 0
            # δ = (thresh - B_eq[i]) / ΔB_ana[i]
            push!(δs, (thresh - B_eq[i]) / ΔB_ana[i])
        end
    end

    return minimum(δs)  # the smallest positive δ that drives someone to `thresh`
end

# ─── USAGE EXAMPLE ─────────────────────────────────────────────────────────────

# your interaction matrix A, scaling matrix ε, and equilibrium B_eq
# …build or load them here…
# let’s say R=20, C=6, S=R+C
# A = …; ε = …; B_eq = …; xi_cons = …  # length-6 vector of consumer thresholds

# if you want to weight consumers by their own ξ’s:
δ_crit = acp_threshold(A, ε, B_eq, xi_cons)

println("Critical press (equal):  ", δ_unweighted)
println("Critical press (by ξᵢ): ", δ_weighted)

δcrit = acp(A, ε, B_eq, xi_cons)

# 90% of that increase ⇒ should still survive:
xi_test  = xi_cons .* (1 + δ_crit*1.2)
p_test   = (p[1], p[2], p[3], xi_test, p[5], p[6], p[7], p[8])
_, _, _, pers90, _, _ = simulate_press_perturbation(B_eq, p_test, (0.0, 500.0), 250.0, 0.0; cb=cb)
@show pers90  # should be 1.0 (no extinctions)

# 110% of that increase ⇒ should start to lose species:
xi_test2 = xi_cons .* (1+ δ_crit*10.0)
p_test2   = (p[1], p[2], p[3], xi_test2, p[5], p[6], p[7], p[8])
_, _, _, pers110, _, _ = simulate_press_perturbation(B_eq, p_test2, (0.0, 500.0), 250.0, 0.0; cb=cb)
@show pers110  # should drop < 1.0

# your analytical result
δ_crit = acp(A, ε, B_eq, xi_cons)   # e.g. -0.4055

# the positive fraction by which to *increase* the consumer thresholds
δp = abs(δ_crit)

# test at  0.9·δp  and 1.1·δp
for α in 0.1:0.1:2.0
  xi_test   = xi_cons .* (1 .+ α*δp)
  p_test    = (p[1], p[2], p[3], xi_test2, p[5], p[6], p[7], p[8])
  sol       = solve(ODEProblem(trophic_ode!, B_eq, (0.,1e4), p_test),
                    Rodas5(); abstol=1e-12, reltol=1e-12)
  goes_extinct = any(sol.u[end] .< 0.0)
  println("α = ", α, " → extinction? ", goes_extinct)
end

using DifferentialEquations, LinearAlgebra, CairoMakie

function compare_sensitivity(p, B_eq; δ=1e-3)
    # unpack
    R, C, m_cons, xi_cons, r_res, d_res, ε, A = p
    S = R + C

    # 1) full Jacobian & sensitivity operator
    
    V        = -inv(I - A)

    # 2) analytic response
    press    = vcat(zeros(R), xi_cons)         # press = +δ in each consumer ξ
    ΔB_ana   = V * press                # length S

    # 3) tiny simulated press
    xi2      = xi_cons .* (1 .+ δ)
    p2       = (R, C, m_cons, xi2, r_res, d_res, ε, A)
    prob     = ODEProblem(trophic_ode!, B_eq, (0., 500.), p2)
    sol2     = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-12)
    B2       = sol2.u[end]
    ΔB_sim   = (B2 .- B_eq)/δ                     # already scaled by δ

    # 4) make the scatter + correlation
    fig = Figure(resolution=(500,500))
    ax  = Axis(fig[1,1];
               xlabel="Analytic ΔB", ylabel="Simulated ΔB",
               title = "δ = $(δ),  corr = $(round(cor(ΔB_ana, ΔB_sim), digits=3))")
    scatter!(ax, ΔB_ana, ΔB_sim; markersize=6, color=:steelblue)
    # 1:1 line
    # lines!(ax, -minimum(abs.(ΔB_ana))*1.1:0.0:maximum(abs.(ΔB_ana))*1.1,
    #           -minimum(abs.(ΔB_ana))*1.1:0.0:maximum(abs.(ΔB_ana))*1.1;
    #           color=:black, linestyle=:dash)
    display(fig)

    return ΔB_ana, ΔB_sim
end

# ─ Example usage within your code ─
# assume `p` is your calibrated parameter tuple, `B_eq` the equilibrium
Δana, Δsim = compare_sensitivity(p, B_eq; δ=1e-3)
cor(Δana, Δsim)


using DifferentialEquations

"""
    persistence_before_press(u0, p, t_perturb; solver=Tsit5(), cb)

Simulate `du/dt = trophic_ode!(u,p)` from t=0 to t=t_perturb,
and return `(pers, B_pre)` where

  • pers = fraction of species with B_pre[i] > 0  
  • B_pre = state vector at t_perturb

This lets you plug in ANY p = (R,C,m_cons,xi_cons,r_res,d_res,ε,A),
even one where you’ve already bumped xi_cons, and see exactly
what survives up to t_perturb.
"""
function persistence_before_press(
    u0, p, t_perturb;
    solver=Tsit5()
)
    cb=build_callbacks(length(u0),0.0)
    prob = ODEProblem(trophic_ode!, u0, (0., t_perturb), p)
    sol  = solve(prob, solver; callback=cb_no_trigger36, reltol=1e-8, abstol=1e-8)
    B_pre = sol.u[end]
    pers  = count(x->x>EXTINCTION_THRESHOLD, B_pre) / length(B_pre)
    return pers, B_pre
end
function persistence_before_press_thresh(u0, p, t_perturb; thresh=EXTINCTION_THRESHOLD)
    # build a callback that kills species when they go ≤ thresh
    cb = build_callbacks(length(u0), thresh)
    prob = ODEProblem(trophic_ode!, u0, (0., t_perturb), p)
    sol  = solve(prob, Tsit5(); callback=cb, reltol=1e-8, abstol=1e-8)
    B_pre = sol.u[end]
    pers  = count(x-> x>thresh, B_pre) / length(B_pre)
    return pers
end
# suppose p0 = your original (R,C,m_cons,xi_cons,r_res,d_res,ε,A)
# and you want to test your “+δ” bump by *modifying* xi_cons yourself:

# build a bumped p:
xi_bumped = xi_cons .* (1 .+ δ_crit * 0.9)  # e.g. 90% of analytic δ*
p_bumped  = (p[1], p[2], p[3], xi_bumped, p[5], p[6], p[7], p[8])

# now just see who survives up to t_perturb:
pers90, B_pre90 = persistence_before_press(vcat(R_eq,C_eq), p_bumped, 250.0)
@show pers90   # now you'll see <1.0 if some went extinct

# same for 110%:
xi_bumped2 = xi_cons .* (1 .+ δ_crit * 1.1)
p_bumped2  = (p[1], p[2], p[3], xi_bumped2, p[5], p[6], p[7], p[8])

pers110, _ = persistence_before_press(vcat(R_eq,C_eq), p_bumped2, 250.0)
@show pers110  # should drop below 1.0

