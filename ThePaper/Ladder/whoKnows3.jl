
# 0. IMPORT YOUR MODEL FUNCTIONS (make_A, calibrate_params, compute_jacobian, trophic_ode!)
# ————————————————————————————————————————————————————————————————
# fix system size
R = 20
C = 4

# target equilibrium
R_eq = fill(10.0, R)
C_eq = fill(1.0, C)

# baseline mortality & self–regulation
m_cons₀ = fill(0.1, C)
d_res₀  = fill(1.0, R)

# sweep grids
connectances   = [0.05, 0.1, 0.2, 0.5, 0.8]
scale_As       = [0.01, 0.1, 0.5, 1.0, 2.0]
scale_epsilons = [0.01, 0.1, 0.2, 0.4]

found = false
tol   = 1.0
cb = build_callbacks(R+C, EXTINCTION_THRESHOLD)
combos = collect(Iterators.product(connectances, scale_As, scale_epsilons))
for u in combos
    conn, sA, sε = u
    # 1) Build consumption matrix (ER model) and scale it
    A = zeros(R+C, R+C)
    A = make_A(A, R, conn, :ER; pareto_exponent=5.0)
    A .*= sA       # amplify or weaken interactions

    # 2) Build efficiency matrix and scale
    ε = fill(sε, R+C, R+C)

    # 3) Calibrate thresholds ξ and resource growth r
    pcal = (R, C, m_cons₀, d_res₀, ε, A)
    xi_cons, r_res = calibrate_params(R_eq, C_eq, pcal;
                                      xi_threshold = 0.7,
                                      constraints   = true)
    # skip if calibration failed
    if any(isnan.(xi_cons)) || any(isnan.(r_res))
        continue
    end

    # 4) Pack params & simulate un-perturbed ODE
    p        = (R, C, m_cons₀, xi_cons, r_res, d_res₀, ε, A)
    B_eq          = vcat(R_eq, C_eq)
    prob0         = ODEProblem(trophic_ode!, B_eq, (0.0, 1e4), p)
    sol0          = solve(prob0, Rodas5();  callback = cb, abstol = 1e-12, reltol = 1e-12)
    B_post0       = sol0.u[end]

    D, Mstar = compute_jacobian(B_post0, p)
    J         = D*Mstar

    ev = eigvals(J)
    if maximum(real(ev)) ≥ 0
        println("Eigenvalues not all negative: ", ev)
        continue
    else
        println("All eigenvalues negative: ", ev)
    end

    begin
        fig = Figure(; size = (600, 400))
        ax  = Axis(fig[1, 1],
                   xlabel = "Time",
                   ylabel = "Species abundances",
                   title  = "Time evolution of species abundances")
        for i in 1:R
            lines!(ax, sol0.t, sol0[i, :], label = "Resources", color = :blue)
        end
        for i in R+1:R+C
            lines!(ax, sol0.t, sol0[i, :], label = "Consumers", color = :red)
        end
        lines!(ax, sol0.t, fill(0.0, length(sol0.t)); label = "Equilibrium", color = :black, linewidth = 2)
        display(fig)
    end

    # 5) Check return‐to‐equilibrium
    err = all(B_post0 .> 0.5)
    err = maximum(abs.(B_post0 .- B_eq)) < tol
    
    if err
        println("Stable system found: connectance=$conn, A×$sA, ε×$sε")
        println("error: ", round(err, digits=3), " | connectance=$conn, A×$sA, ε×$sε")
        println(A[end, :])
        global found = true
        global B_post0 = B_post0
        global A = A
        global p = p
        global xi_cons = xi_cons
        global r_res = r_res
        global d_res₀ = d_res₀
        global ε = ε
        global m_cons₀ = m_cons₀
        global sol0 = sol0
        break
    end
end

if !found
    error("No stable combination found in the given grid.")
end

# ————————————————————————————————————————————————————————————————
# 6. ONCE STABLE, COMPUTE ANALYTICAL VS SIMULATED RESPONSE
# ————————————————————————————————————————————————————————————————
B_eq = B_post0
# rebuild Jacobian pieces
D, Mstar = compute_jacobian(B_eq, p)
n        = size(Mstar,1)
I_mat    = Matrix{Float64}(I, n, n)        # true Float64 identity
A_star   = Mstar .+ I_mat                  # nondimensional A*

# sensitivity matrix V = –(I – A*)⁻¹
V        = -inv(I_mat .- A_star)

# analytical ΔB for uniform press on consumers
δξ       = 0.5
press    = vcat(zeros(R), ones(C))
ΔB_ana   = V * press * δξ

# simulate pressed system
xi2      = xi_cons .+ δξ
params2  = (R, C, m_cons₀, xi2, r_res, d_res₀, ε, A)
prob2    = ODEProblem(trophic_ode!, B_eq, (0.0, 1e4), params2)
sol2     = solve(prob2, Rodas5(); callback = cb, abstol = 1e-12, reltol = 1e-12)
B_post2  = sol2.u[end]
ΔB_sim   = (B_post2 .- B_eq) ./ δξ

# compute correlation
rval = cor(ΔB_ana, ΔB_sim)
println("Correlation (analytic vs sim): ", round(rval, digits=3))

begin 
    fig = Figure(; size = (600, 400))
    
    prob = ODEProblem(trophic_ode!, B_eq, (0.0, 1e4), p)
    sol0 = solve(prob, saveat = 1e3)
    
    ax  = Axis(
        fig[1, 1],
        xlabel = "Time",
        ylabel = "Biomass"
    )
    for i in 1:R
        lines!(ax, sol0.t, sol0[i, :], label = "Resources", color = :blue)
    end
    for i in R+1:R+C
        lines!(ax, sol0.t, sol0[i, :], label = "Consumers", color = :red)
    end
    ax2 = Axis(
        fig[1, 2],
        xlabel = "Time",
        ylabel = "Biomass"
    )
    for i in 1:R
        lines!(ax2, sol2.t, sol2[i, :], label = "Resources", color = :blue)
    end
    for i in R+1:R+C
        lines!(ax2, sol2.t, sol2[i, :], label = "Consumers", color = :red)
    end
    lines!(ax2, sol2.t, fill(0.0, length(sol2.t)); label = "Equilibrium", color = :black, linewidth = 2)
    lines!(ax, sol0.t, fill(0.0, length(sol0.t)); label = "Equilibrium", color = :black, linewidth = 2)
    display(fig)
end

# ————————————————————————————————————————————————————————————————
# 7. PLOT WITH MAKIE (all inside begin…end)
# ————————————————————————————————————————————————————————————————
begin
    mn  = min(minimum(ΔB_ana), minimum(ΔB_sim))
    mx  = max(maximum(ΔB_ana), maximum(ΔB_sim))

    fig = Figure(; size = (600,400))
    ax  = Axis(fig[1,1];
               xlabel = "Analytical ΔB",
               ylabel = "Simulated ΔB",
               title  = "Analytic vs Sim Biomass Shifts")

    for i in 1:R
        scatter!(ax, ΔB_ana[i], ΔB_sim[i]; markersize = 8, color = :blue)
    end
    for i in R+1:R+C
        scatter!(ax, ΔB_ana[i], ΔB_sim[i]; markersize = 8, color = :red)
    end
    lines!(ax, [mn,mx], [mn,mx]; linestyle = :dash, linewidth = 2)
    display(fig)
end
