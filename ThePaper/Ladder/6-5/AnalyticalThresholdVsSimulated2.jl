using DifferentialEquations, LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# (1) Analytic critical press (additive in ξ)
function acp_threshold(
    A::AbstractMatrix,
    ε::AbstractMatrix,
    B_eq::AbstractVector,
    xi_cons::AbstractVector;
    thresh = EXTINCTION_THRESHOLD
)
    S = size(A, 1)
    C = length(xi_cons)
    R = S - C

    # 1a) zero out ε on non‐feeding links, build A*
    eps_eff = copy(ε)
    eps_eff[A .<= 0.0] .= 0.0
    A_star = eps_eff .* A

    # 1b) sensitivity operator V = −(I − A*)⁻¹
    V = -inv(I(S) .- A_star)

    # 2) unit‐press vector: +1 per consumer, 0 per resource
    press = vcat(zeros(R), ones(C))

    # 3) analytic per‐unit response
    ΔB_ana = V * press

    # 4) for each i with ΔB_ana[i] < 0, solve B_eq[i] + δ ΔB_ana[i] = thresh
    δs = Float64[]
    for i in 1:S
        if ΔB_ana[i] < 0
            push!(δs, (thresh - B_eq[i]) / ΔB_ana[i])
        end
    end

    return minimum(δs)
end

# ─────────────────────────────────────────────────────────────────────────────
# (2) Empirical critical press by bisection (additive in ξ)
function find_sim_thresholdd(
    A, ε, B_eq, xi_cons, m_cons, r_res, d_res;
    lower=0.0, upper=1.0,
    tol=1e-4, maxiter=30,
    t_perturb=250.0,
    thresh=EXTINCTION_THRESHOLD
)
    R = length(B_eq) - length(xi_cons)

    pers = δ -> begin
      xi2 = xi_cons .+ δ
      p2  = (R, length(xi_cons), m_cons, xi2, r_res, d_res, ε, A)
      return persistence_before_press_thresh(
        vcat(B_eq), p2, t_perturb; thresh=thresh
      )
    end

    @assert pers(lower) ≈ 1.0 "already below 1 at δ=lower"
    # expand upper until we actually get pers(upper)<1
    while pers(upper) ≥ 1.0 && upper < 1e3
      upper *= 2
    end

    for i in 1:maxiter
      mid = (lower + upper)/2
      if pers(mid) ≥ 1.0
        lower = mid
      else
        upper = mid
      end
      if upper - lower < tol
        break
      end
    end

    return lower
end

# ─────────────────────────────────────────────────────────────────────────────
# (3) Compare analytic vs. simulated
function compare_thresholds(
    A, ε, B_eq, xi_cons, m_cons, r_res, d_res;
    thresh=EXTINCTION_THRESHOLD, t_perturb=250.0
)
    δ_ana = acp_threshold(A, ε, B_eq, xi_cons; thresh=thresh)
    δ_sim = find_sim_thresholdd(A, ε, B_eq, xi_cons, m_cons, r_res, d_res;
                               lower=0.0, upper=δ_ana*2,
                               tol=1e-3, maxiter=30,
                               t_perturb=t_perturb)
    return δ_ana, δ_sim
end

# ─────────────────────────────────────────────────────────────────────────────
# USAGE EXAMPLE (you must have defined persistence_before_press earlier)
# 
δ_ana, δ_sim = compare_thresholds(
    A, ε, B_eq, xi_cons,
    p[3], p[5], p[6];
    thresh=EXTINCTION_THRESHOLD,
    t_perturb=250.0
)
println("Analytic δ* = $δ_ana")
println("Simulated δ* = $δ_sim")
