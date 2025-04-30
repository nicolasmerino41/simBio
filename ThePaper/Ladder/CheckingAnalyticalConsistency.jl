# — assumes you have already defined:
#     make_A!(A,R,conn,:PL;pareto_exponent)
#     trophic_ode!(du,u,p,t)
#     calibrate_params(R_eq,C_eq,R,C,ε,A)
#     compute_jacobian_factors(B_eq,p)  # returns (D,Mstar)

function check_analytical_vs_sim(
  B_eq, p; δ=1e-3, tpress=250., tspan=(0.,500.0)
)
    R,C,m_cons,ξ_cons,r_res,d_res,ε,A = p
    S = R+C

    # 1) analytical: D, M*  &  b
    D, Mstar = compute_jacobian(B_eq, p)
    b        = vcat(zeros(R), ones(C))*δ
    δB_ana   = Mstar \ b

    # 2) simulate up to tpress
    sol1 = solve(ODEProblem(trophic_ode!, B_eq, (0.,tpress), p),
                 Tsit5(); reltol=1e-8, abstol=1e-8)
    B_pre = sol1.u[end]

    # 3) apply press to ξ_cons and re-simulate
    xi2 = ξ_cons .* (1 .- δ)
    p2  = (R,C,m_cons, xi2, r_res, d_res, ε, A)
    sol2 = solve(ODEProblem(trophic_ode!, B_pre, (tpress,tspan[2]), p2),
                 Tsit5(); reltol=1e-8, abstol=1e-8)
    B_post = sol2.u[end]

    δB_sim = (B_post .- B_eq)./δ

    # 4) error metrics
    err  = δB_sim .- δB_ana
    rmse = sqrt(mean(err.^2))
    ρ    = cor(δB_sim, δB_ana)

    return δB_ana, δB_sim, rmse, ρ
end

function run_and_plot(R, C; conn=0.3, σA=2.0, δ=1e-3)
    S = R+C
    A = zeros(S,S)

    # 1) find one feasible, stable equilibrium
    B_eq = nothing; p = nothing
    while true
      make_A(A, R, conn, :ER; pareto_exponent=σA)

      # sample a target equilibrium
      R_eq = fill(10.0, R)
      C_eq = fill(1.0, C)
      # R_eq = abs.(rand(LogNormal(log(1.0)-0.2^2/2,0.2), R))
      # C_eq = abs.(rand(LogNormal(log(0.1)-0.2^2/2,0.2), C))
      ξ, r = calibrate_params(R_eq, C_eq, (R, C, m_cons, d_res, ε, A); xi_threshold=0.7, constraints=true)
      if any(isnan.(ξ)) || any(isnan.(r))
        continue
      end

      B_eq = vcat(R_eq, C_eq)
      p    = (R, C, fill(0.1,C), ξ, r, fill(1.0,R), ones(S,S), A)
      # local stability check
      D, Mstar = compute_jacobian(B_eq, p)   # or D*Mstar
      J = D * Mstar
      if maximum(real(eigvals(J))) < 0
        break
      end
    end

    # 2) run the analytical vs sim comparison
    δana, δsim, rmse, ρ = 
      check_analytical_vs_sim(B_eq, p; δ=δ)

    # 3) plot
    fig = Figure(resolution=(600,600))
    ax  = Axis(fig[1,1],
      xlabel="δBₐₙₐ", ylabel="δBₛᵢₘ",
      title="ρ=$(round(ρ,digits=3)), RMSE=$(round(rmse,digits=3))"
    )
    scatter!(ax, δana, δsim; color=:blue, markersize=8)
    lines!(ax,
      [minimum(δana), maximum(δana)],
      [minimum(δana), maximum(δana)];
      color=:black, linestyle=:dash
    )
    display(fig)

    return (δana=δana, δsim=δsim, rmse=rmse, ρ=ρ)
end

# — finally, call it:
res = run_and_plot(2, 1; conn=1.0, σA=3.0, δ=0.001)
res.δana