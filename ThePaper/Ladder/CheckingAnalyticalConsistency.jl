"""
    check_analytical_vs_sim(S, R, C, A, B_eq, p; δ=0.1, tspan=(0.0,500.0), t_press=250.0)

Compute and compare the analytical press response (via inverting J) vs. the simulated one.
Returns δB_ana, δB_sim, RMSE and correlation ρ.
"""
function check_analytical_vs_sim(S, R, C, A, B_eq, p;
                                  δ=0.1, tspan=(0.0,500.0), t_press=250.0)

    # 1) Build Jacobian J at B_eq
    J = compute_jacobian(B_eq, p)

    # 2) Build the "press-vector" D·b
    total = R + C
    Db = zeros(total)
    for k in 1:C
        Db[R + k] = (p[3][k] / p[4][k]) * δ   # m_cons[k]/xi_cons[k] times δ
    end

    # 3) Analytical press response: δB_ana = J^{-1}·(D·b)
    δB_ana = J \ Db

    # 4) Simulated press
    # 4a) run up to t_press
    sol1 = solve(ODEProblem(trophic_ode!, B_eq, (0.0, t_press), p),
                 Tsit5(); reltol=1e-8, abstol=1e-8)
    B_pre = sol1.u[end]

    # 4b) apply press (scale xi_cons by 1-δ)
    R_, C_, m_cons, xi_cons, r_res, d_res, ε, A_mat = p
    xi2 = xi_cons .* (1 .- δ)
    p2  = (R, C, m_cons, xi2, r_res, d_res, ε, A_mat)

    sol2 = solve(ODEProblem(trophic_ode!, B_pre, (t_press, tspan[2]), p2),
                 Tsit5(); reltol=1e-8, abstol=1e-8)
    B_post = sol2.u[end]

    δB_sim = (B_post .- B_eq) ./ δ

    # 5) compare
    err = δB_sim .- δB_ana
    return (
      δB_ana = δB_ana,
      δB_sim = δB_sim,
      rmse   = sqrt(mean(err .^ 2)),
      ρ      = cor(δB_sim, δB_ana)
    )
end

# -------------------------------
# Example of driving it once:
# -------------------------------
function find_and_test(R, C;   
  conn=0.5, α=5.0, δ=0.1,
  max_tries=1_000)

  S = R+C
  for attempt in 1:max_tries
    # 1) build A (here power-law PL(α))
    A = zeros(S,S)
    ks = clamp.(round.(Int, rand(Pareto(1,α), C)), 1, clamp(round(Int,conn*R),1,R))
    for (idx,k) in enumerate(ks)
      ci = R + idx
      js = sample(1:R, k; replace=false)
      for j in js
        w = abs(randn())
        A[ci,j] = w; A[j,ci] = -w
      end
    end

    # 2) draw equilibria & calibrate
    R_eq = abs.(rand(LogNormal(log(1.0)-0.2^2/2,0.2), R))
    C_eq = abs.(rand(LogNormal(log(0.1)-0.2^2/2,0.2), C))
    ε    = clamp.(rand(LogNormal(0.5,0.1), S,S),0,1)
    xi, r_res = calibrate_params(R_eq, C_eq, (R,C, fill(0.1,C), fill(1.0,R), ε, A);
                                 xi_threshold=0.3, constraints=true)

    if any(isnan, xi) || any(isnan, r_res)
      continue   # bad calibration → retry
    end

    # 3) build the Jacobian and test stability
    B_eq = vcat(R_eq,C_eq)
    J    = compute_jacobian(B_eq, (R,C, fill(0.1,C), xi, r_res, fill(1.0,R), ε, A))
    if maximum(real.(eigvals(J))) ≥ 0
      continue   # unstable → retry
    end

    # 4) we have a valid config! run the analytical vs   sim check
    res = check_analytical_vs_sim(S, R, C, A, B_eq,
                                  (R,C, fill(0.1,C), xi, r_res, fill(1.0,R), ε, A);
                                  δ=δ, t_press=250.0)
    return res
  end

  error("Could not find a valid stable, NaN-free configuration in $max_tries tries")
end

# Example usage:
R, C = 10, 5
result = find_and_test(R, C; conn=0.5, α=5.0, δ=0.01)
println("ρ = ", result.ρ, ", RMSE = ", result.rmse)