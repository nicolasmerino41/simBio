# 0. IMPORT YOUR MODEL FUNCTIONS (make_A, calibrate_params, compute_jacobian, trophic_ode!)
# ————————————————————————————————————————————————————————————————
# fix system size
R = 20
C = 4

# target equilibrium
# R_eq = fill(100.0, R)
# C_eq = fill(10.0, C)
abundance_mean = 100.0
R_eq = abs.(rand(Normal(abundance_mean, abundance_mean*0.1), R))
C_eq = abs.(rand(Normal(abundance_mean*0.1, abundance_mean*0.01), C))
# C_eq[4] = 0.0
# baseline mortality & self–regulation
m_cons₀ = fill(0.1, C)
d_res₀  = fill(1.0, R)

# sweep grids
connectances   = [0.05, 0.1, 0.2, 0.5, 0.8]
scale_As       = [0.01, 0.1, 0.5, 1.0, 2.0]
scale_epsilons = [0.01, 0.1, 0.2, 0.4]

found = false
tol   = 10.0
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
    prob0         = ODEProblem(trophic_ode!, B_eq, (0.0, 1000), p)
    sol0          = solve(prob0, Rodas5();  callback = cb, abstol = 1e-12, reltol = 1e-12)
    B_post0       = sol0.u[end]

    D, Mstar = compute_jacobian(B_post0, p)
    J         = D*Mstar

    ev = eigvals(J)
    if maximum(real(ev)) > 0
        println("Eigenvalues not all negative: ", ev)
        continue
    else
        println("All eigenvalues negative: ", ev)
    end

    begin
        fig = Figure(; size = (1100, 580))
        ax  = Axis(fig[1, 1],
                   xlabel = "Time",
                   ylabel = "Species abundances",
                   title  = "Time evolution of species abundances",
                   yticks = 0:maximum(R_eq)/10:maximum(R_eq)
                   )
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
        global A_star = A
        break
    end
end

if !found
    error("No stable combination found in the given grid.")
end

# ————————————————————————————————————————————————————————————————
# 6. ONCE STABLE, COMPUTE ANALYTICAL VS SIMULATED RESPONSE
# ————————————————————————————————————————————————————————————————
B_eq = sol0.u[end]
# 1) get the Jacobian at the “true” equilibrium B_post0
D, Mstar   = compute_jacobian(B_eq, p)
J          = D * Mstar

I_mat    = I(R+C)  # identity matrix (R+C)×(R+C)

# 2) build the press‐vector properly: ∂f_C/∂ξ_i = –B_C* at equilibrium
press      = zeros(R+C)
press[R+1:R+C] .= 1.0

# 3) analytic per‐unit sensitivity (for δξ = +1)
# ΔB_ana_unit = -inv(J) * press    # no further division needed, since δξ=1

# press_vec is length R+C, carrying the +1’s for consumers and 0’s for resources
press_vec = (zeros(R+C))    # length R+C
press_vec[R+1:R+C] .= 1.0
# press_vec = vcat(zeros(R), -B_eq[R+1:end])

# pick a delta ξ
δξ           = -1.0

# analytic per-unit sensitivity
V            = -inv(I_mat .- A_star)    # (R+C)×(R+C)
ΔB_ana_unit  = V * press_vec          # length R+C

# to simulate, we only update the C-vector ξ_cons:
xi2          = xi_cons .+ δξ           # length C

# pack parameters (ξ_cons is always length C!)
p2           = (R, C, m_cons₀, xi2, r_res, d_res₀, ε, A_star)

# simulate and get the new full-S equilibrium
sol2         = solve(ODEProblem(trophic_ode!, B_eq, (0.0,1e4), p2), Rodas5();
                     callback=cb, abstol=1e-12, reltol=1e-12)
B_post2      = sol2.u[end]
if any(B_post2 .== 0.0)
    @warn "Some species went extinct, this will fuck up the sensitivity analysis."
end

# simulated per-unit response (length R+C)
ΔB_sim_unit = (B_post2 .- B_eq) ./ δξ 

zero_idx = findall(iszero.(ΔB_sim_unit))
new_ΔB_sim_unit = copy(ΔB_sim_unit)
new_ΔB_ana_unit = copy(ΔB_ana_unit)                 
new_ΔB_sim_unit[zero_idx] .= NaN
new_ΔB_ana_unit[zero_idx] .= NaN
filter!(x -> !isnan(x), new_ΔB_sim_unit)
filter!(x -> !isnan(x), new_ΔB_ana_unit)

# now both ΔB_ana_unit and ΔB_sim_unit are length R+C and directly comparable
println("r   = ", round(cor(new_ΔB_ana_unit, new_ΔB_sim_unit), digits=3))
# println("slope = ", round(fit(Line, ΔB_ana_unit, ΔB_sim_unit).coeffs[2], digits=3))
begin 
    fig = Figure(; size = (1100, 580))
    
    prob = ODEProblem(trophic_ode!, B_eq, (0.0, 1000), p)
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
# 7. PLOT PER‐UNIT SENSITIVITIES WITH MAKIE
# ————————————————————————————————————————————————————————————————
begin
    lo = min(minimum(ΔB_ana_unit), minimum(ΔB_sim_unit))
    hi = max(maximum(ΔB_ana_unit), maximum(ΔB_sim_unit))

    fig = Figure(; size = (600,400))
    ax  = Axis(fig[1,1];
               xlabel = "Analytic ΔB per unit ξ",
               ylabel = "Simulated ΔB per unit ξ",
               title  = "Per‐unit Sensitivity: Analytic vs Sim")

    for i in 1:R
        scatter!(ax, ΔB_ana_unit[i], ΔB_sim_unit[i]; color = :blue, markersize = 5)
    end
    for i in R+1:R+C
        scatter!(ax, ΔB_ana_unit[i], ΔB_sim_unit[i]; color = :red, markersize = 5)
    end

    lines!(ax, [lo,hi], [lo,hi];
           linestyle = :dash, linewidth = 1, color = :black)

    display(fig)
end
