using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads, nthreads, threadid

# ----------------------------------------------------------------------------
# 0) Assumed helper functions:
#    make_A, calibrate_params, compute_jacobian, trophic_ode!,
#    build_callbacks(S, threshold), EXTINCTION_THRESHOLD
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 1) Check local stability: max real part of Jacobian < 0
# ----------------------------------------------------------------------------
is_locally_stable(J) = maximum(real(eigvals(J))) < 0

# ----------------------------------------------------------------------------
# 2) Simulate unperturbed dynamics to check 100% persistence
# ----------------------------------------------------------------------------
function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    if sol.t[end] < tspan[2] return false, sol.u[end] end
    Bf = sol.u[end]
    return all(Bf .> EXTINCTION_THRESHOLD), Bf
end

# ----------------------------------------------------------------------------
# 3) Compute analytic vs simulated sensitivity correlation (consumers only)
# ----------------------------------------------------------------------------
function sensitivity_correlation(A, ε, p, fixed; δξ=1.0, cb, plot=false)
    R, C = p[1], p[2]
    S = R + C
    # analytic operator
    Astar = ε .* A .- transpose(A)
    V     = -inv(I(S) .- Astar)
    press = vcat(zeros(R), ones(C))
    ΔB_ana = V * press * δξ
    # simulate bump ξ_cons
    xi2 = copy(p[4]) .+ δξ
    p2  = (R, C, p[3], xi2, p[5], p[6], ε, A)
    prob = ODEProblem(trophic_ode!, fixed, (0.,500.), p2)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    B2   = sol.u[end]
    extinct = any(B2 .< EXTINCTION_THRESHOLD)
    ΔB_sim = (B2 .- fixed) ./ δξ
    # mask zeros
    idx = findall(!iszero, ΔB_sim)
    if extinct
        return NaN, extinct
    else
        if plot
            lo = min(minimum(ΔB_ana), minimum(ΔB_sim))
            hi = max(maximum(ΔB_ana), maximum(ΔB_sim))

            fig = Figure(; size = (600,400))
            ax  = Axis(fig[1,1];
                    xlabel = "Analytic ΔB per unit ξ",
                    ylabel = "Simulated ΔB per unit ξ",
                    title  = "Per‐unit Sensitivity: Analytic vs Sim")

            for i in 1:R
                scatter!(ax, ΔB_ana[i], ΔB_sim[i]; color = :blue, markersize = 5)
            end
            for i in R+1:R+C
                scatter!(ax, ΔB_ana[i], ΔB_sim[i]; color = :red, markersize = 5)
            end

            lines!(ax, [lo,hi], [lo,hi];
                linestyle = :dash, linewidth = 1, color = :black)

            display(fig)
        end
        return cor(ΔB_ana[idx], ΔB_sim[idx]), extinct
    end
end

function analytical_vs_sim(p, B_eq; δξ=1.0, cb, plot=false)
    
    R, C, m_cons, xi_cons, r_res, d_res, ε, A = p
    D, Mstar   = compute_jacobian(B_eq, p)
    J          = D * Mstar

    I_mat    = I(R+C)  # identity matrix (R+C)×(R+C)

    # press_vec is length R+C, carrying the +1’s for consumers and 0’s for resources
    press_vec = (zeros(R+C))    # length R+C
    press_vec[R+1:R+C] .= 1.0
    
    # analytic per-unit sensitivity
    # 1a) zero out ε on non‐feeding links, build A*
    eps_eff = copy(ε)
    eps_eff[A .<= 0.0] .= 0.0
    A_star = eps_eff .* A

    V            = -inv(I_mat .- A_star)    # (R+C)×(R+C)
    ΔB_ana_unit  = V * press_vec          # length R+C

    # to simulate, we only update the C-vector ξ_cons:
    xi2          = xi_cons .+ δξ           # length C

    # pack parameters (ξ_cons is always length C!)
    p2           = (R, C, m_cons, xi2, r_res, d_res, ε, A_star)

    # simulate and get the new full-S equilibrium
    sol2         = solve(ODEProblem(trophic_ode!, B_eq, (0.0,1e4), p2), Rodas5();
                        callback=cb, abstol=1e-12, reltol=1e-12)
    B_post2      = sol2.u[end]
    
    if any(B_post2 .== 0.0)
        return NaN, true
    else
        # simulated per-unit response (length R+C)
        ΔB_sim_unit = (B_post2 .- B_eq) ./ δξ 
        if plot
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
        return cor(ΔB_ana_unit, ΔB_sim_unit), false
    end
end

# ----------------------------------------------------------------------------
# 4) Main pipeline
# ----------------------------------------------------------------------------
function community_pipeline(; 
    S_vals       = [30, 40, 50],
    conn_vals    = 0.05:0.05:0.3,
    C_ratios     = [0.1,0.2,0.3],
    IS_vals      = [0.01,0.1,1.0],
    d_vals       = [0.1, 1.0],
    m_vals       = [0.25],
    eps_vals     = [0.1, 0.5],
    skew_vals    = [1.0, 0.1],
    abundance_mean = 10.0,
    abundance_types = [:Normal, :Log],
    scenarios    = [:ER, :PL, :MOD],
    pareto_exps  = [1.5, 3.0],
    mod_gammas   = [5.0, 10.0],
    single_combin_repeats = 1,
    max_communities = 1000,
    max_calib = 10,
    plot_sensitivity = false
)
    # build parameter tuples
    combos = Tuple[]
    for S in S_vals, conn in conn_vals, C_ratio in C_ratios,
        IS in IS_vals, d in d_vals, m in m_vals,
        eps in eps_vals, skew in skew_vals,
        abund in abundance_types, scenario in scenarios,
        scr in 1:single_combin_repeats

        if scenario == :PL
            for px in pareto_exps
                push!(combos, (S,conn,C_ratio,IS,d,m,eps,skew,abund,scr,scenario,px,0.0))
            end
        elseif scenario == :MOD
            for γ in mod_gammas
                push!(combos, (S,conn,C_ratio,IS,d,m,eps,skew,abund,scr,scenario,0.0,γ))
            end
        else
            push!(combos, (S,conn,C_ratio,IS,d,m,eps,skew,abund,scr,scenario,0.0,0.0))
        end
    end
    println("$(length(combos)) parameter combinations")
    results = Vector{NamedTuple}()
    locki = ReentrantLock()

    @threads for idx in sample(1:length(combos), min(length(combos), max_communities); replace=false)
        S,conn,C_ratio,IS,d,m,eps_mean,skew,abund,scr,scenario,pex,modg = combos[idx]
        C = clamp(round(Int, S*C_ratio), 1, S-1)
        R = S - C
        cb = build_callbacks(S, EXTINCTION_THRESHOLD)

        # 3) calibrate
        m_cons = fill(m, C)
        d_res  = fill(d, R)
        
        # prepare placeholders
        xi_cons, r_res = fill(NaN, C), fill(NaN, R)
        p = nothing           # will hold the final good parameter tuple
        R_eq_final = nothing  # placeholders for the equilibrium draws
        C_eq_final = nothing

        tries = 0
        while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
            tries += 1

            # (re)build A, ε, draw R_eq, C_eq …
            A = make_A(zeros(S,S), R, conn, scenario; pareto_exponent=pex, mod_gamma=modg) .* IS
            ε = clamp.(rand(Normal(eps_mean, eps_mean), S, S), 0.0, 1.0)

            if abund == :Log
                R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2,
                                        abundance_mean), R))
                C_eq = abs.(rand(LogNormal(log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2,
                                        abundance_mean*skew), C))
            else
                R_eq = abs.(rand(Normal(abundance_mean, abundance_mean*0.1), R))
                C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew*0.1), C))
            end

            # calibrate
            xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C, m_cons, d_res, ε, A);
                                            xi_threshold=0.7,
                                            constraints=true)

            if !any(isnan, xi_cons)
                # capture the working values *into the outer scope*
                p           = (R, C, m_cons, xi_cons, r_res, d_res, ε, A)
                R_eq_final  = R_eq
                C_eq_final  = C_eq
                break       # exit the while
            end
        end

        # if we never got a good calibration, skip
        if isnothing(p)
            continue
        end

        # now use p, R_eq_final, C_eq_final with confidence…
        fixed = vcat(R_eq_final, C_eq_final)
        R_eq, C_eq = R_eq_final, C_eq_final
        if any(isnan, xi_cons) continue end
        
        A  = p[8]
        ε  = p[7]

        # p = (R,C,m_cons,xi_cons,r_res,d_res,ε,A)

        # 4) check local stability
        D, Mstar = compute_jacobian(fixed, p)
        J_simp = D*Mstar
        if !is_locally_stable(J_simp) continue end

        # 5) simulate unperturbed persistence
        ok, B_eq = survives!(fixed, p; cb=cb)
        if !ok continue end

        # 5.5) collectivity & May‐complexity
        φ, C_may = collectivity_metrics(A, ε, R)

        # 5.7) temporal unpredictability
        unpredict = temporal_unpredictability(A, ε, p, B_eq;
                                            δ=1e-2, t_early=50.0, t_total=500.0)

        # 6) structural descriptors
        g = SimpleGraph(A .!= 0)
        degs = degree(g)
        degree_cv = std(degs)/mean(degs)
        link_count = count(!iszero, A[R+1:S, :])
        connectance = link_count / (C*(S-1))
        pyramid_ratio = mean(R_eq)/mean(C_eq)

        # 5.6) perturbation‐depth
        fw = floyd_warshall_shortest_paths(g)
        distmat = fw.dists  # S×S Int matrix (typemax = unreachable)
        depth = mean_knockout_depth(S, (0.0, 500.0), distmat, B_eq, p)

        # 7) sensitivity correlation
        # ρ, ext = sensitivity_correlation(A, ε, p, B_eq; cb=cb, plot=plot_sensitivity)
        # println("xi_cons: ", xi_cons)
        ρ, ext = analytical_vs_sim(p, B_eq; cb=cb, plot=plot_sensitivity, δξ=-1.0)
        if ext 
            @warn "extinction happened so scorr would not be correct"
            continue
        end
        # println("ρ: ", ρ, " when xi_cons: ", p[4], "r_res:", p[5], "and B_eq: ", B_eq, "\n")
        rec = (
            S=S, conn=conn, C_ratio=C_ratio, IS=IS, d=d, m=m, ε=eps_mean,
            scr=scr, C=C, R=R, R_eq=R_eq, C_eq=C_eq,
            scenario=scenario, pareto_exponent=pex, mod_gamma=modg,
            skew=skew, abundance_dist=abund,
            degree_cv=degree_cv, connectance=connectance, pyramid_ratio=pyramid_ratio,
            # ← newly added:
            φ=φ, C_may=C_may, depth=depth, unpredict=unpredict,
            sens_corr=ρ, xi_cons=xi_cons, r_res=r_res, ε_mat=ε, A_mat=A
        )

        lock(locki) do
            push!(results, rec)
        end
    end

    return DataFrame(results)
end

# === Run the pipeline ===
A4 = community_pipeline(
    S_vals = [30],
    conn_vals = 0.1:0.05:0.4,
    C_ratios = [0.2],
    IS_vals = [0.1, 1.0],
    d_vals = [1.0],
    m_vals = [0.1],
    eps_vals = [0.1],
    skew_vals = [0.001,0.01,0.05,0.1,0.5,0.85],
    abundance_types = [:Normal, :Log], # [:Normal, :Log],
    abundance_mean = 100.0,
    scenarios = [:ER,:PL,:MOD],
    pareto_exps = [1.5, 3.0],
    mod_gammas = [5.0, 1.0],
    max_communities = 1000,
    single_combin_repeats = 2.0,
    max_calib = 10,
    plot_sensitivity = false,
)

serialize("ThePaper/Ladder/Outputs/A4.jls", A4)
df_community_pipeline = A4

begin
    
    # 1) Define colors and markers
    scenario_colors = Dict(:ER => :blue, :PL => :green, :MOD => :red)
    is_markers      = Dict(1.0 => :circle, 0.1 => :utriangle)

    # 2) Set up a 1×4 panel figure
    fig = Figure(resolution = (1200, 400))
    ax1 = Axis(fig[1, 1]; title="Sens. Corr. vs Connectance",
                        xlabel="Connectance", ylabel="Sensitivity Correlation")
    ax2 = Axis(fig[1, 2]; title="Sens. Corr. vs Degree CV",
                        xlabel="Degree CV")
    ax3 = Axis(fig[1, 3]; title="Sens. Corr. vs Pyramid Ratio",
                        xlabel="Mean R / Mean C")
    ax4 = Axis(fig[1, 4]; title="Sens. Corr. vs Skew",
                        xlabel="Skew")

    # 3) Scatter by (scenario, IS)
    for scen in keys(scenario_colors)
        for isval in keys(is_markers)
            mask = (df_community_pipeline.scenario .== scen) .& (df_community_pipeline.IS .== isval)
            if !any(mask) continue end

            scatter!(
            ax1,
            df_community_pipeline.connectance[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval],
            label  = "$(scen), IS=$(isval)"
            )
            scatter!(
            ax2,
            df_community_pipeline.degree_cv[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
            scatter!(
            ax3,
            df_community_pipeline.pyramid_ratio[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
            scatter!(
            ax4,
            df_community_pipeline.skew[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
        end
    end

    # 4) Add a legend to the top-right panel
    #    We give it the collected plot handles & labels
    Legend(fig[1, 5], handles, labels; title = "Scenario, IS", 
        orientation = :vertical, position = :rt)

    fig
end

metrics = [
    (:φ,         "Collectivity ϕ"),
    (:C_may,     "May complexity Cₘₐy"),
    (:depth,     "Perturbation depth"),
    (:unpredict, "Temporal unpredictability"),
]

# plot sens_corr vs each metric, 4 panels in a 2×2 grid
begin
    
    # 1) Define colors & markers
    scenario_colors = Dict(:ER => :blue, :PL => :green, :MOD => :red)
    is_markers      = Dict(1.0 => :circle, 0.1 => :utriangle)

    # 2) Create a 2×2 figure
    fig = Figure(resolution = (800, 600))
    axφ     = Axis(fig[1, 1]; title = "Sens. Corr. vs Collectivity φ",
                            xlabel = "Collectivity φ", ylabel = "Sensitivity corr.")
    axCmay  = Axis(fig[1, 2]; title = "Sens. Corr. vs May complexity Cₘₐy",
                                xlabel = "May complexity Cₘₐy")
    axDepth = Axis(fig[2, 1]; title = "Sens. Corr. vs Perturbation depth",
                                xlabel = "Perturbation depth", ylabel = "Sensitivity corr.")
    axUnp   = Axis(fig[2, 2]; title = "Sens. Corr. vs Temporal unpredictability",
                                xlabel = "Temporal unpredictability")

    # 3) Plot each (scenario, IS) group and collect handles
    handles = Plot[]  # will store all the scatter plot objects
    labels  = String[]

    for scen in keys(scenario_colors)
        for isval in keys(is_markers)
            mask = (df_community_pipeline.scenario .== scen) .& (df_community_pipeline.IS .== isval)
            if !any(mask) continue end

            push!(handles, scatter!(
            axφ,
            df_community_pipeline.φ[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            ))
            push!(labels, "$(scen), IS=$(isval)")

            scatter!(
            axCmay,
            df_community_pipeline.C_may[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
            scatter!(
            axDepth,
            df_community_pipeline.depth[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
            scatter!(
            axUnp,
            df_community_pipeline.unpredict[mask],
            df_community_pipeline.sens_corr[mask];
            color  = scenario_colors[scen],
            marker = is_markers[isval]
            )
        end
    end

    # 4) Add a legend to the top-right panel
    #    We give it the collected plot handles & labels
    Legend(fig[1, 3], handles, labels; title = "Scenario, IS", 
        orientation = :vertical, position = :rt)

    display(fig)

end
