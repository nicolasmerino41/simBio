# 1) pick the failing cases
fails = df[df.sens_corr .< 0.9, :]
println("Found $(nrow(fails)) low–correlation runs; plotting up to 6 of them")

# 2) helper to plot one case
function plot_case(row)
    # reconstruct all the pieces
    S, C = row.S, row.C
    R    = S - C
    # pull parameters
    p = (
      row.R, row.C,
      fill(row.m, C),
      row.xi_cons,   # you’ll need to stash xi_cons & r_res in df
      row.r_res,
      fill(row.d, R),
      row.ε_mat,  # likewise, store ε & A in df
      row.A_mat
    )
    # equilibrium
    B_eq = vcat(row.R_eq, row.C_eq) # = row.B_eq

    p1 = p
    prob = ODEProblem(trophic_ode!, B_eq, (0.,1000.), p1)
    sol1  = solve(prob, Rodas5(); callback=cb, abstol=1e-9, reltol=1e-9)
    println("sol_end1 = $(sol1.u[end])")


    # 1) integrate under +1 bump in ξ
    δξ = -1.0
    press_vec = vcat(zeros(R), ones(C))
    # analytic
    Astar = p[7] .* p[8] .- transpose(p[8])
    V     = -inv(I(S) .- Astar)
    ΔB_ana = V * press_vec * δξ

    # simulation
    xi2 = p[4] .+ δξ
    p2  = (p[1],p[2],p[3], xi2, p[5], p[6], p[7], p[8])
    prob = ODEProblem(trophic_ode!, B_eq, (0.,1000.), p2)
    sol2  = solve(prob, Rodas5(); callback=cb, abstol=1e-9, reltol=1e-9)
    println("sol_end2 = $(sol2.u[end])")
    # plot
    fig = Figure(; size=(800,400))
    ax2  = Axis(fig[1,2], title="Run #$(row.scr)  –  ρ=$(round(row.sens_corr, digits=2))",
               xlabel="Time", ylabel="Abundance")
    ax1 = Axis(fig[1,1], title="Biomass trajectories", xlabel="Time", ylabel="Abundance")

    # trajectories
    
    for i in 1:S
        color = i ≤ R ? :blue : :red
        lines!(ax1, sol1.t[1:end], sol1[i, 1:end]; color=color, alpha=0.6, linewidth=1.0) # lines!(ax1, sol1.t[1:10], getindex.(sol1.u, i)[1:10]; color=color, alpha=0.6, linewidth=0.01)
        lines!(ax2, sol2.t[1:end], sol2[i, 1:end]; color=color, alpha=0.6, linewidth=1.0) # lines!(ax, sol2.t[1:10], getindex.(sol2.u, i)[1:10]; color=color, alpha=0.6, linewidth=0.01)
        # horizontal analytic prediction
        # hlines!(ax1, [B_eq[i], B_eq[i]+ΔB_ana[i]]; color=color, linestyle=:dash)
        # hlines!(ax, [B_eq[i], B_eq[i]+ΔB_ana[i]]; color=color, linestyle=:dash)
    end
    display(fig)
end

# 3) plot up to 6 of them
for row in eachrow(first(fails, 10))
    plot_case(row)
end