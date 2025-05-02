include("Ladder4.1_drago1.jl")
using DifferentialEquations, Random, Statistics, DataFrames, CSV
import Base.Threads: @threads, nthreads, threadid
include("ExploringFeasibleSpace.jl")

function make_A(
    A::AbstractMatrix{<:Real},
    R::Int,
    conn::Float64,
    scenario::Symbol;
    pareto_exponent::Float64 = 1.75,
    mod_gamma::Float64       = 5.0
)
    S = size(A,1)
    C = S - R
    # zero out
    fill!(A, 0.0)

    if scenario == :ER
        # Erdős–Rényi random bipartite consumer-resource links
        for i in (R+1):S, j in 1:R
            if rand() < conn
                A[i,j] = abs(rand(Normal()))
                A[j,i] = -abs(rand(Normal()))
            end
        end

    elseif scenario == :PL
        # Power-law consumer out-degree controlled by pareto_exponent
        x_m = 1.0  # minimum links per consumer
        raw_degrees = rand(Pareto(x_m, pareto_exponent), C)
        ks = clamp.(floor.(Int, raw_degrees), 1, R)

        for (idx, k) in enumerate(ks)
            ci = R + idx
            prey_pool = 1:R
            selected = sample(prey_pool, k; replace=false)
            for j in selected
                A[ci, j] = abs(rand(Normal()))
                A[j, ci] = -abs(rand(Normal()))
            end
        end

    elseif scenario == :MOD
        # Bipartite modular: two resource‐modules and two matching consumer‐modules
        halfR = fld(R, 2)
        halfC = fld(C, 2)
        res_block1 = 1:halfR
        res_block2 = (halfR+1):R
        cons_block1 = (R+1):(R+halfC)
        cons_block2 = (R+halfC+1):S
    
        for i in (R+1):S         # each consumer
            for j in 1:R         # only resources
                # do they share the same module?
                same_module = (i in cons_block1 && j in res_block1) ||
                              (i in cons_block2 && j in res_block2)
                # within‐module links denser by mod_gamma
                p_link = same_module ? conn * mod_gamma : conn / mod_gamma
                p_link = clamp(p_link, 0.0, 1.0)
                if rand() < p_link
                    A[i,j] = abs(rand(Normal()))
                    A[j,i] = -abs(rand(Normal()))
                end
            end
        end
    else
        error("Unknown scenario: $scenario")
    end

    return A
end

# # Usage examples:
# A = zeros(S,S)
# a = make_A(A, R, 0.05, :ER)
# b = make_A(A, R, 0.2, :PL; pareto_exponent=1.0)
# c = make_A(A, R, 0.2, :MOD; mod_gamma=1.0)

# begin
#     fig = Figure()
#     ax = Axis(fig[1,1])
#     MK.heatmap!(ax, a, colorrange = (-1.0, 1.0))
#     ax.yreversed = false
#     display(fig)
# end

# g = SimpleGraph(a .!= 0)
# degs = degree(g)
# degree_cv = std(degs) / mean(degs)
# -----------------------------------------------------------
# 3) rppp_guided: only over pre-screened param tuples
# -----------------------------------------------------------
function rppp_guided(df_params::DataFrame, delta_vals;
                     ladder_steps=1:16,
                     tspan=(0.0,50.0), t_perturb=25.0,
                     max_calib=10, atol=10.0,
                     abundance_mean=1.0, combinations=1000,
                     plot_full=false, plot_simple=false,
                     plot_ana_vs_sim_full = false, plot_ana_vs_sim_simp = false
)

    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()
    suffix(step) = step==1 ? "Full" : "S$(step-1)"

    @threads for row in eachrow(df_params)[sample(1:nrow(df_params), min(combinations, nrow(df_params)), replace=false)]
        S        = row.S
        conn     = round(row.conn, digits=3)
        C_ratio  = row.C_ratio
        IS   = row.IS
        d_val    = row.d
        m_val    = row.m
        epsilon_mean   = row.epsilon
        C        = row.C
        R        = row.R
        scenario = row.scenario
        constraints = row.constraints
        ssp = row.ssp
        pexs = row.pexs
        mod_gamma = row.mod_gamma
        abund = row.abundance_distribution
        skew = row.skew

        # pick correct callback
        cb = eval(Symbol("cb_no_trigger$(S)"))

        for delta in delta_vals
            
            try
                # build A & epsilon_full
                A = zeros(S,S)
                
                A = make_A(A, R, conn, scenario; pareto_exponent = pexs, mod_gamma = mod_gamma)
                A = A .* IS
                # epsilon_full = clamp.(rand(LogNormal(epsilon_mean, epsilon_mean), S, S), 0, 1)
                epsilon_full   = clamp.(rand(Normal(epsilon_mean, epsilon_mean), S, S), 0, 1)
                
                # target eq
                # R_eq = abs.(rand(LogNormal(log(abundance_mean)-abundance_mean^2/2, abundance_mean), R))
                # C_eq = abs.(rand(LogNormal(log(abundance_mean*0.1)-(abundance_mean*0.2)^2/2, abundance_mean*0.2), C))
                if abund == :Log
                    R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
                    C_eq = abs.(rand(LogNormal(log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
                elseif abund == :Normal
                    R_eq = abs.(rand(Normal(abundance_mean, abundance_mean), R))
                    C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
                end
                fixed = vcat(R_eq, C_eq)

                m_cons  = fill(m_val, C)
                d_res   = fill(d_val, R)

                # calibrate
                xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,epsilon_full,A); xi_threshold=0.7, constraints=constraints)
                tries=1
                while (any(isnan, xi_cons)||any(isnan, r_res)) && tries<max_calib
                    A .= 0
                    A = make_A(A, R, conn, scenario; pareto_exponent = pexs, mod_gamma = mod_gamma)
                    A = A .* IS
                    if abund == :Log
                        R_eq = abs.(rand(LogNormal(log(abundance_mean) - (abundance_mean^2)/2, abundance_mean), R))
                        C_eq = abs.(rand(LogNormal(log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
                    elseif abund == :Normal
                        R_eq = abs.(rand(Normal(abundance_mean, abundance_mean), R))
                        C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
                    end

                    fixed = vcat(R_eq, C_eq)
                    xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,epsilon_full,A); xi_threshold=0.7, constraints=constraints)
                    tries+=1
                end
                if any(isnan, xi_cons)||any(isnan, r_res)
                    continue
                end

                # full model solve & metrics
                p_full = (R,C,m_cons,xi_cons,r_res,d_res,epsilon_full,A)

                prob   = ODEProblem(trophic_ode!, fixed, tspan, p_full)
                sol    = solve(prob,Tsit5(); callback=cb, reltol=1e-8, abstol=1e-8)
                if sol.t[end] < t_perturb || any(isnan, sol.u[end]) ||
                    any(isinf, sol.u[end]) || any([!isapprox(sol.u[end][i], vcat(R_eq, C_eq)[i], atol=atol) for i in 1:S])
                    @warn "Error: solution did not finish properly"
                    continue
                end
                if any(sol.u[end] .< EXTINCTION_THRESHOLD)
                    @warn "Error: some species went extinct before even starting the perturbation"
                    continue
                end
                B_eq = sol.u[end]

                g = SimpleGraph(A .!= 0)
                degs = degree(g)
                degree_cv = std(degs) / mean(degs)

                meanB = mean(B_eq)
                varB = var(B_eq)
                RelVar = varB / (meanB^2 + 1e-8)
                
                resi = compute_resilience(B_eq, p_full)
                reac = compute_reactivity(B_eq, p_full)
                # press perturb
                rt_full, os_full, ire_full, _, B2, sp_rt_full = simulate_press_perturbation(
                    fixed, p_full, tspan, t_perturb, delta;
                    solver=Tsit5(), plot=plot_full,
                    show_warnings=true, full_or_simple=true,
                    cb = cb,
                    species_specific_perturbation=ssp
                )
                #############################################################
                # FIXING CORRELATION OF SENSITIVITIES
                # 1) get the Jacobian at the “true” equilibrium B_post0
                D, Mstar   = compute_jacobian(B_eq, p_full)
                J_full          = D * Mstar

                I_mat    = I(R+C)  # identity matrix (R+C)×(R+C)

                # 2) build the press‐vector properly: ∂f_C/∂ξ_i = –B_C* at equilibrium
                press_vec      = zeros(R+C)
                press_vec[R+1:R+C] .= 1.0

                # analytic per-unit sensitivity
                V            = -inv(I_mat .- A)    # (R+C)×(R+C)
                deltaB_ana  = V * press_vec          # length R+C

                # simulated per-unit response (length R+C)
                deltaB_sim = (B2 .- B_eq) ./ delta
                
                zero_idx = findall(iszero.(deltaB_sim))
                new_deltaB_sim = copy(deltaB_sim)
                new_deltaB_ana = copy(deltaB_ana)                 
                new_deltaB_sim[zero_idx] .= NaN
                new_deltaB_ana[zero_idx] .= NaN
                filter!(x -> !isnan(x), new_deltaB_sim)
                filter!(x -> !isnan(x), new_deltaB_ana)

                sens_corr_full = round(cor(new_deltaB_ana, new_deltaB_sim), digits=3)
                if sens_corr_full < 0.8
                    continue
                end
                if plot_ana_vs_sim_full
                    lo = min(minimum(deltaB_ana), minimum(deltaB_sim))
                    hi = max(maximum(deltaB_ana), maximum(deltaB_sim))

                    fig = Figure(; size = (600,400))
                    ax  = Axis(fig[1,1];
                            xlabel = "Analytic ΔB per unit ξ",
                            ylabel = "Simulated ΔB per unit ξ",
                            title  = "Per‐unit Sensitivity: Analytic vs Sim in FULL")

                    for i in 1:R
                        scatter!(ax, deltaB_ana[i], deltaB_sim[i]; color = :blue, markersize = 5)
                    end
                    for i in R+1:R+C
                        scatter!(ax, deltaB_ana[i], deltaB_sim[i]; color = :red, markersize = 5)
                    end

                    lines!(ax, [lo,hi], [lo,hi];
                        linestyle = :dash, linewidth = 1, color = :black)

                    display(fig)
                end
                #############################################################
                # record metrics for each step
                metrics = Dict("Full" => (
                    return_time=mean(filter(!isnan, rt_full)),
                    overshoot  =mean(filter(!isnan, os_full)),
                    ire        =mean(filter(!isnan, ire_full)),
                    resilience =resi,
                    reactivity =reac,
                    sens_corr  =sens_corr_full,
                    before_p   =mean(sol.u[end] .> EXTINCTION_THRESHOLD),
                    after_p    =mean(B2 .> EXTINCTION_THRESHOLD),
                    sp_rt_mean     =mean(filter(!isnan, sp_rt_full)),
                    sp_rt_cv       =std(filter(!isnan, sp_rt_full)) / mean(filter(!isnan, sp_rt_full))
                    ))

                # ladder steps 
                for step in ladder_steps[2:end]
                    suf = suffix(step)
                    A_s, epsilon_s = transform_for_ladder_step(step, A, epsilon_full)
                    p_simp   = (R,C,m_cons,xi_cons,r_res,d_res,epsilon_s,A_s)

                    sol2  = solve(
                        ODEProblem(
                            trophic_ode!, fixed, tspan, p_simp
                        ),
                        Tsit5(); 
                        callback=cb, reltol=1e-8, abstol=1e-8
                    )
                    if sol2.t[end] < t_perturb || any(isnan, sol2.u[end]) || any(isinf, sol2.u[end])
                        metrics[suf] = (
                            return_time=NaN, overshoot=NaN, ire=NaN,
                            resilience=NaN, reactivity=NaN,
                            sens_corr=NaN, before_p=NaN, after_p=NaN,
                            sp_rt_mean = NaN, sp_rt_cv = NaN
                            
                        )
                        continue
                    end

                    B_eq2 = sol2.u[end]
                    before_p = mean(sol2.u[end] .> EXTINCTION_THRESHOLD)
                    res2  = compute_resilience(B_eq2,p_simp)
                    rea2  = compute_reactivity(B_eq2,p_simp)
                    # get the full per-species press matrix just like in the full model
                    rt2, os2, ire2, _, B2_simp, sp_rt_simp = simulate_press_perturbation(
                            fixed, p_simp, tspan, t_perturb, delta;
                            solver=Tsit5(), plot=plot_simple,
                            show_warnings=true,
                            full_or_simple=false,
                            cb = cb,
                            species_specific_perturbation=ssp
                        )
    
                    ##########################################################################
                    # 1) get the Jacobian at the “true” equilibrium B_post0
                    D, Mstar   = compute_jacobian(B_eq2, p_simp)
                    J_simp          = D * Mstar

                    I_mat    = I(R+C)  # identity matrix (R+C)×(R+C)

                    # 2) build the press‐vector properly: ∂f_C/∂ξ_i = –B_C* at equilibrium
                    press_vec      = zeros(R+C)
                    press_vec[R+1:R+C] .= 1.0

                    # analytic per-unit sensitivity
                    V            = -inv(I_mat .- A_s)    # (R+C)×(R+C)
                    deltaB_ana  = V * press_vec          # length R+C

                    # simulated per-unit response (length R+C)
                    deltaB_sim = (B2_simp .- B_eq2) ./ delta

                    zero_idx = findall(iszero.(deltaB_sim))
                    new_deltaB_sim = copy(deltaB_sim)
                    new_deltaB_ana = copy(deltaB_ana)                 
                    new_deltaB_sim[zero_idx] .= NaN
                    new_deltaB_ana[zero_idx] .= NaN
                    filter!(x -> !isnan(x), new_deltaB_sim)
                    filter!(x -> !isnan(x), new_deltaB_ana)

                    corr_simp = round(cor(new_deltaB_ana, new_deltaB_sim), digits=3)
                    if plot_ana_vs_sim_simp
                        lo = min(minimum(deltaB_ana), minimum(deltaB_sim))
                        hi = max(maximum(deltaB_ana), maximum(deltaB_sim))

                        fig = Figure(; size = (600,400))
                        ax  = Axis(fig[1,1];
                                xlabel = "Analytic ΔB per unit ξ",
                                ylabel = "Simulated ΔB per unit ξ",
                                title  = "Per‐unit Sensitivity: Analytic vs Sim in step $suf")

                        for i in 1:R
                            scatter!(ax, deltaB_ana[i], deltaB_sim[i]; color = :blue, markersize = 5)
                        end
                        for i in R+1:R+C
                            scatter!(ax, deltaB_ana[i], deltaB_sim[i]; color = :red, markersize = 5)
                        end

                        lines!(ax, [lo,hi], [lo,hi];
                            linestyle = :dash, linewidth = 1, color = :black)

                        display(fig)
                    end
                    ##########################################################################
                    # compound error, step
                    # t_pred2 = mean(filter(!isnan, rt2))
                    # t_obs2  = mean(filter(!isnan, rt2))
                    # err_T2   = abs(t_obs2 - t_pred2) / ((t_obs2 + t_pred2)/2)
                    # os_pred2 = mean(filter(!isnan, os_full))
                    # os_obs = mean(filter(!isnan, os2))
                    # err_OS2  = abs(os_obs - os_pred2) / ((os_obs + os_pred2)/2)
                    # ire_pred2 = mean(filter(!isnan, ire_full))
                    # ire_obs = mean(filter(!isnan, ire2))
                    # err_IRE2  = abs(ire_obs - ire_pred2) / ((ire_obs + ire_pred2)/2)
                    # comp_err2 = (err_T2 + err_OS2 + err_IRE2)/3

                    metrics[suf] = (
                        # compound_error = comp_err2,
                        return_time    = mean(filter(!isnan, rt2)),
                        overshoot      = mean(filter(!isnan, os2)),
                        ire            = mean(filter(!isnan, ire2)),
                        resilience     = res2,
                        reactivity     = rea2,
                        sens_corr      = corr_simp,
                        before_p       = before_p,
                        after_p        = mean(B2_simp .> EXTINCTION_THRESHOLD),
                        sp_rt_mean         = mean(filter(!isnan, sp_rt_simp)),
                        sp_rt_cv           = std(filter(!isnan, sp_rt_simp)) / mean(filter(!isnan, sp_rt_simp))
                    )
                end

                # assemble record
                base = (
                    S=S, conn=conn, C_ratio=C_ratio,
                    IS=IS, d=d_val, m=m_val, epsilon=epsilon_mean,
                    delta=delta, R=R, C=C, degree_cv=degree_cv, RelVar=RelVar, 
                    scenario=scenario, constraints=constraints, ssp=ssp, pexs = pexs, mod_gamma = mod_gamma,
                    abundance_distribution = abund, skewness = skew
                )
                
                rec = base
                for step in ladder_steps
                    suf = suffix(step)
                    m = metrics[suf]
                    rec = merge(rec,
                        (; Symbol("rt_$suf")      => m.return_time,
                        Symbol("os_$suf")      => m.overshoot,
                        Symbol("ire_$suf")     => m.ire,
                        Symbol("res_$suf")     => m.resilience,
                        Symbol("rea_$suf")     => m.reactivity,
                        Symbol("scorr_$suf")   => m.sens_corr,
                        Symbol("bef_$suf")     => m.before_p,
                        Symbol("aft_$suf")     => m.after_p,
                        Symbol("sp_rt_$suf")   => m.sp_rt_mean,
                        Symbol("sp_rt_cv_$suf")   => m.sp_rt_cv)
                    )
                end

                lock(results_lock) do
                    push!(results, rec)
                end

            catch e
                @warn "Error in rppp_guided" row delta e
            end # try
        end  # delta
    end  # threads

    return DataFrame(results)
end
# -----------------------------------------------------------
# 4) Example invocation
# -----------------------------------------------------------
cb_no_trigger30 = build_callbacks(30, EXTINCTION_THRESHOLD)
# cb_no_trigger40 = build_callbacks(40, EXTINCTION_THRESHOLD)
# cb_no_trigger50 = build_callbacks(50, EXTINCTION_THRESHOLD)
# cb_no_trigger60 = build_callbacks(60, EXTINCTION_THRESHOLD)
# cb_no_trigger70 = build_callbacks(70, EXTINCTION_THRESHOLD)
# cb_no_trigger80 = build_callbacks(80, EXTINCTION_THRESHOLD)
# cb_no_trigger90 = build_callbacks(90, EXTINCTION_THRESHOLD)
# cb_no_trigger200 = build_callbacks(200, EXTINCTION_THRESHOLD)

begin
    S_vals = [30]
    conn_vals   = range(0.01, 0.4, length=10)
    C_ratios    = [0.1, 0.2]
    IS_vals     = [0.01, 0.1, 1.0]
    d_vals      = [0.1, 1.0, 2.0]
    mort_vals   = [0.25]
    epsilon_vals= [0.01, 0.1, 0.2, 0.5]
    delta_vals  = [0.1]

    df_feas = feasibility_search(
        S_vals, conn_vals, C_ratios, IS_vals,
        d_vals, mort_vals, epsilon_vals, delta_vals;
	    degree_distribution_types = [:ER, :PL, :MOD],
        pareto_exponents        = [1.1, 2.0, 5.0],   # only used for :PL
        mod_gammas              = [1.0, 5.0, 10.0],  # only used for :MOD
        pyramid_skewness = [0.5, 0.1, 0.01, 0.001],
        abundance_distribution = [:Log, :Normal],
        tspan=(0.0,500.0),
        t_perturb=250.0,
        max_calib=10,
        abundance_mean=10.0,
        atol=1.0,
        xi_threshold=0.7,
        number_of_combos=100000,
        calibrate_params_constraints = true,
        species_specific_perturbation = false
    )
end

# (b) Extract good combos
df_good = unique(
    df_feas[df_feas.feasible .== true,
            [:S, :conn, :C_ratio, :IS, :d, :m, :epsilon, :R, :C, :scenario, :constraints, :ssp, :pexs, :mod_gamma, :scorr, :skew, :abundance_distribution]]
)
df_good = unique(
    df_good[df_good.scorr .> 0.8, [:S, :conn, :C_ratio, :IS, :d, :m, :epsilon, :R, :C, :scenario, :constraints, :ssp, :pexs, :mod_gamma, :skew, :abundance_distribution]]
)

# (c) Guided rppp
begin
    delta_vals = [0.1]
    df_results = rppp_guided(
        df_good, delta_vals;
        ladder_steps=1:16,
        max_calib=15, # 10 is good
        tspan=(0.0,500.0), t_perturb=250.0,
        atol=10.0,
        abundance_mean=10.0,
        combinations=10000, # You want 10000
        plot_full=false,
        plot_simple=false,
        plot_ana_vs_sim_full = false,
        plot_ana_vs_sim_simp = false
    )
end

# (d) Save or inspect
CSV.write("guided_results_ERPLMOD_skew_logabund_constraint_s30_c3and6.csv", df_results)

c1and3and6 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERPLMOD_skew_logabund_constraint_s30_c1and3and6.csv") |> DataFrame
c9and12 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERPLMOD_skew_logabund_constraint_s30_c9and12.csv") |> DataFrame
df = vcat(c1and3and6, c9and12)