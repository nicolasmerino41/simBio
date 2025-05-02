# ----------------------------------------------------------------
# Precompute analytic and simulated responses for each species
# ----------------------------------------------------------------
δr  = 1.0                # resource growth perturbation step
δξ  = 1.0                # consumer threshold perturbation step
S   = R + C              # total species count
T   = 500.0

# 1) build the linear‐response operator V = (I - A*)^-1
I_mat = I(R+C)
V     = -inv(I_mat .- A_star)

# 2) Allocate storage
ΔB_ana = zeros(S, S)
ΔB_sim = zeros(S, S)

# 3) Loop over each “parameter index” i = 1…S
for i in 1:S
    press = zeros(S)

    if i ≤ R
        # resource‐growth perturbation on species i
        press[i] = 1.0
        V     = inv(I_mat .- A_star)
        ΔB_ana[:, i] = V * press * δr

        # simulate by bumping r_i
        r2 = copy(r_res)
        r2[i] += δr
        p2 = (R, C, m_cons₀, xi_cons, r2, d_res₀, ε, A)
    else
        # consumer‐threshold perturbation on consumer k = i−R
        k = i - R
        press[i] = 1.0
        V     = -inv(I_mat .- A_star)
        ΔB_ana[:, i] = V * press * δξ

        # simulate by bumping xi_cons[k]
        xi2 = copy(xi_cons)
        xi2[k] += δξ
        p2 = (R, C, m_cons₀, xi2, r_res, d_res₀, ε, A)
    end

    sol2 = solve(
        ODEProblem(trophic_ode!, B_eq, (0.0, T), p2),
        Rodas5(); callback=cb, abstol=1e-12, reltol=1e-12
    )
    B_post = sol2.u[end]
    ΔB_sim[:, i] = (B_post .- B_eq) ./ (i ≤ R ? δr : δξ)
    s_corr = cor(ΔB_sim[:, i], ΔB_ana[:, i])
    println("For species $i, the correlation is $s_corr and the s")
end

# ----------------------------------------------------------------
# Plot analytic vs simulated for each perturbation
# ----------------------------------------------------------------
begin
    cols = C
    rows = ceil(Int, S/cols)
    fig  = Figure(; size = (1100, 750))

    for i in 1:24 #S
        row = ceil(Int, i/cols)
        col = i - (row-1)*cols
        ax  = Axis(
            fig[row, col];
            title = i ≤ R ? "Δr on Res $i" : "Δξ on Con $(i-R)",
            xlabel = "Analytic ΔB",
            ylabel = "Simulated ΔB"
        )

        # mask out zero (extinct) responses
        valid = .!iszero.(ΔB_sim[:, i])
        idxs  = findall(valid)

        # separate resource vs consumer points
        res_idx  = filter(j -> j ≤ R, idxs)
        cons_idx = filter(j -> j >  R, idxs)

        # resources in blue
        scatter!(
            ax,
            ΔB_ana[res_idx,  i],
            ΔB_sim[res_idx,  i];
            color      = :blue,
            markersize = 5,
            label      = "Resources"
        )

        # consumers in red
        scatter!(
            ax,
            ΔB_ana[cons_idx, i],
            ΔB_sim[cons_idx, i];
            color      = :red,
            markersize = 7,
            label      = "Consumers"
        )

        # identity line
        lo = minimum(vcat(ΔB_ana[valid, i], ΔB_sim[valid, i]))
        hi = maximum(vcat(ΔB_ana[valid, i], ΔB_sim[valid, i]))
        lines!(
            ax,
            [lo, hi],
            [lo, hi];
            linestyle = :dash,
            color     = :black
        )

        # axislegend(ax)
    end

    display(fig)
end

begin
    # baseline unperturbed solution
    prob0 = ODEProblem(trophic_ode!, B_eq, (0.0, 1000), p)
    sol0  = solve(prob0, Rodas5(); saveat=1.0,
                  callback=cb, abstol=1e-12, reltol=1e-12)

    # set up figure with C rows and 2 columns
    fig = Figure(; size = (800, 650))
    for k in 1:C
        # 1) simulate the press on consumer k
        xi2 = copy(xi_cons)
        xi2[k] += 1.0                # one‐unit bump
        p2  = (R, C, m_cons₀, xi2, r_res, d_res₀, ε, A)
        sol2 = solve(
            ODEProblem(trophic_ode!, B_eq, (0.0, 1000), p2),
            Rodas5(); saveat=1.0,
            callback=cb, abstol=1e-12, reltol=1e-12
        )

        # 2) plot before (left) and after (right)
        ax1 = Axis(fig[k, 1],
            title = "Consumer $k: pre‐press",
            xlabel = "Time", ylabel = "Biomass")
        ax2 = Axis(fig[k, 2],
            title = "Consumer $k: post‐press",
            xlabel = "Time", ylabel = "Biomass")

        # resources in blue
        # for i in 1:R
        #     lines!(ax1, sol0.t, sol0[i, :], color=:blue)
        #     lines!(ax2, sol2.t, sol2[i, :], color=:blue)
        # end

        # consumers in red
        for i in R+1:R+C
            lines!(ax1, sol0.t, sol0[i, :], color=:red)
            lines!(ax2, sol2.t, sol2[i, :], color=:red)
        end

        # equilibrium line
        lines!(ax1, sol0.t, fill(0.0, length(sol0.t));
               color=:black, linestyle=:dash)
        lines!(ax2, sol2.t, fill(0.0, length(sol2.t));
               color=:black, linestyle=:dash)
    end

    display(fig)
end