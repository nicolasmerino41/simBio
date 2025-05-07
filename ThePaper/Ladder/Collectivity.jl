# collectivity_pipeline.jl
# assumes: make_A, calibrate_params, compute_jacobian, trophic_ode!,
#           simulate_press_perturbation, transform_for_ladder_step,
#           build_callbacks, EXTINCTION_THRESHOLD

# ----------------------------------------------------------------------------
# 1. collectivity & May‐complexity
# ----------------------------------------------------------------------------
function collectivity_metrics(A, ε, R::Int)
    S = size(A,1)
    C = S - R

    # (1) spectral radius & spread of the “effective” interaction matrix A*
    Astar = ε .* A .- transpose(A)
    ϕ     = maximum(abs, eigvals(Astar))
    σ     = std(vec(Astar))

    # (2) empirical connectance among consumers
    #    count all nonzero outgoing links from each consumer
    link_count = count(!iszero, A[R+1:S, :])
    #    each consumer has (S-1) possible targets (no self‐loops)
    p = link_count / (C * (S-1))

    # (3) May’s composite complexity threshold
    C_may = σ * sqrt(p * (S-1))

    return ϕ, C_may
end

# ----------------------------------------------------------------------------
# 2. perturbation depth
# ----------------------------------------------------------------------------
### Helper: mean distance to “moved” species under single knock‐out ###
global function mean_knockout_depth(S, tspan, distmat, B_eq, p_full)
    depths = Float64[]
    for i in 1:S
    u0 = copy(B_eq); u0[i] = 0.0
    sol = solve(ODEProblem(trophic_ode!, u0, (0.,tspan[2]), p_full),
                Tsit5(); callback=cb_no_trigger30,
                abstol=1e-12, reltol=1e-12)
    moved = findall(j->abs(sol.u[end][j] - B_eq[j]) > 1e-6, 1:S)
    if !isempty(moved)
        ds = filter(!isinf, distmat[i, moved])
        if !isempty(ds)
        push!(depths, mean(ds))
        end
    end
    end
    return isempty(depths) ? NaN : mean(depths)
end

# ----------------------------------------------------------------------------
# 3. temporal unpredictability
# ----------------------------------------------------------------------------
function temporal_unpredictability(A, ε, p, B_eq;
                                   δ=1e-2, t_early=50., t_total=500.)
    # unpack R and C from the parameters tuple:
    R, C = p[1], p[2]
    S = R + C

    # 1) build A* = ε.*A - Aᵀ and the corresponding V = - (I - A*)⁻¹
    Astar = ε .* A .- transpose(A)
    V = -inv(Matrix{Float64}(I, S, S) .- Astar)

    # 2) analytic response to a uniform small bump δ in consumer thresholds
    #    press vector should be [0 for resources; 1 for consumers]
    press = vcat(zeros(R), ones(C))
    ΔB_ana = V * press * δ

    # 3) simulate the same bump
    #    p == (R, C, m_cons, xi_cons, r_res, d_res, ε, A)
    xi2 = copy(p[4]) .+ δ
    p2  = (R, C, p[3], xi2, p[5], p[6], ε, A)

    sol = solve(
      ODEProblem(trophic_ode!, B_eq, (0., t_total), p2),
      Tsit5(), saveat=[t_early, t_total],
      callback = cb_no_trigger30,
      abstol=1e-12, reltol=1e-12
    )

    early = sol.u[1] .- B_eq
    late  = sol.u[2] .- B_eq
    return 1 - cor(early, late)
end

# ----------------------------------------------------------------------------
# 4. biotic‐niche contribution
# ----------------------------------------------------------------------------
function biotic_niche_contribution(B_eq, R_eq, C_eq)
    K = vcat(R_eq, C_eq)
    N = B_eq ./ K
    return sum((N .- 1).^2) / sum(N.^2)
end

# ----------------------------------------------------------------------------
# 5. full ladder + collectivity
# ----------------------------------------------------------------------------
function collectivity_ladder(
    df_good::DataFrame;
    abundance_mean=100.0,
    tspan=(0.,500.), t_perturb=250.0,
    xi_threshold=0.7, max_calib=10,
    ladder_steps=1:16
)
    records = Vector{NamedTuple}()
    results_lock = ReentrantLock()

    @threads for row in eachrow(df_good)
        ### unpack parameters ###
        S,conn,C_ratio,IS_val,d_val,m_val,
            eps_mean,scenario,pexs,modg,skew,abund,B_term,IS =
            row.S, row.conn, row.C_ratio, row.IS, row.d, row.m,
            row.epsilon, row.scenario, row.pexs, row.mod_gamma,
            row.skew, row.abundance_distribution, row.B_term, row.IS
        
            C = clamp(round(Int, S*C_ratio),1,S-1)
        R = S - C

        ### 1) Build A & ε ###
        A = make_A(zeros(S,S), R, conn, scenario;
                    pareto_exponent=pexs, mod_gamma=modg, B_term=B_term) .* IS_val
        ε = clamp.(rand(Normal(eps_mean,eps_mean), S, S), 0, 1)

        ### 2) Draw target equilibrium ###
        if abund == :Log
            R_eq = abs.(rand(LogNormal(log(abundance_mean)-(abundance_mean^2)/2, abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*skew)-((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
        else
            R_eq = abs.(rand(Normal(abundance_mean,abundance_mean), R))
            C_eq = abs.(rand(Normal(abundance_mean*skew,abundance_mean*skew), C))
        end
        fixed = vcat(R_eq, C_eq)

        ### 3) Calibrate r & ξ ###
        m_cons = fill(m_val, C)
        d_res  = fill(d_val, R)
        xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,ε,A);
                                            xi_threshold=xi_threshold,
                                            constraints=true)
        tries = 1
        while (any(isnan,xi_cons)||any(isnan,r_res)) && tries<max_calib
            A .= 0
            A = make_A(A, R, conn, scenario;
                        pareto_exponent=pexs, mod_gamma=modg, B_term=B_term) .* IS_val
            # redraw eq
            if abund == :Log
            R_eq = abs.(rand(LogNormal(log(abundance_mean)-(abundance_mean^2)/2, abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*skew)-((abundance_mean*skew)^2)/2, abundance_mean*skew), C))
            else
            R_eq = abs.(rand(Normal(abundance_mean,abundance_mean), R))
            C_eq = abs.(rand(Normal(abundance_mean*skew,abundance_mean*skew), C))
            end
            fixed = vcat(R_eq, C_eq)
            xi_cons, r_res = calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,ε,A);
                                            xi_threshold=xi_threshold,
                                            constraints=true)
            tries += 1
        end
        if any(isnan, xi_cons) 
            continue
        end

        ### 4) Solve full‐model equilibrium ###
        p_full = (R, C, m_cons, xi_cons, r_res, d_res, ε, A)
        sol0 = solve(ODEProblem(trophic_ode!, fixed, tspan, p_full),
                        Tsit5(); callback=cb_no_trigger30,
                        abstol=1e-12, reltol=1e-12)
        B_eq = sol0.u[end]

        ### 5) Build graph & precompute all‐pairs shortest‐paths ONCE ###
        g = SimpleGraph(S)
        for i in 1:S, j in 1:S
            if A[i,j] != 0
            add_edge!(g, i, j)
            end
        end
        fw = floyd_warshall_shortest_paths(g)
        distmat = fw.dists  # S×S Int matrix (typemax = unreachable)

        ### 6) Collect full‐model metrics ###
        φ_full, C_may_full = collectivity_metrics(A, ε, R)
        depth_full        = mean_knockout_depth(S, tspan, distmat, B_eq, p_full)
        unpredict_full    = temporal_unpredictability(A, ε, p_full, B_eq)
        niche_full        = biotic_niche_contribution(B_eq, R_eq, C_eq)

        ### 7) Step through ladder once (or loop 1:16) ###
        for step in ladder_steps[2:end]
            A_s, ε_s = transform_for_ladder_step(step, A, ε)
            φ_s, C_may_s      = collectivity_metrics(A_s, ε_s, R)
            depth_s           = mean_knockout_depth(S, tspan, distmat, B_eq, p_full)
            unpredict_s       = temporal_unpredictability(A_s, ε_s, p_full, B_eq)
            niche_s           = niche_full  # unchanged by A‐only steps

            lock(results_lock) do
                push!(records, (
                    S           = S,
                    step        = step,
                    φ           = step==1 ? φ_full     : φ_s,
                    C_may       = step==1 ? C_may_full : C_may_s,
                    depth       = step==1 ? depth_full : depth_s,
                    unpredict   = step==1 ? unpredict_full : unpredict_s,
                    niche       = step==1 ? niche_full : niche_s,
                    B_term      = B_term,
                    conn        = conn,
                    C_ratio     = C_ratio,
                    IS          = IS,
                    epsilon     = eps_mean,
                    scenario    = scenario,
                    skew        = skew,
                    # abundance_distribution = abund
                ))
            end
        end
    end

    DataFrame(records)
end

# === Example usage ===
df = CSV.File("ThePaper/Ladder/Outputs/collectivity_df.csv") |> DataFrame
# @time collectivity_df = collectivity_ladder(df_good; abundance_mean=10.0, ladder_steps=1:16)
# CSV.write("ThePaper/Ladder/Outputs/collectivity_df.csv", collectivity_df)

