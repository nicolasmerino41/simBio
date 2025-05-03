# collectivity_analysis.jl
using LinearAlgebra, Statistics, Graphs, DifferentialEquations, Random, DataFrames

# === Assumed existing functions ===
# make_A(A, R, conn, scenario; pareto_exponent, mod_gamma)
# calibrate_params(R_eq, C_eq, (R,C,m_cons,d_res,ε,A); xi_threshold, constraints)
# compute_jacobian(B, p) -> (D, Mstar)
# trophic_ode!(du,u,p,t)
# simulate_press_perturbation(fixed, p, tspan, t_perturb, delta; ...)
# transform_for_ladder_step(step, A, ε) -> (A_s, ε_s)

# ----------------------------------------------------------------
# 1) Spectral radius φ and May‐complexity C
# ----------------------------------------------------------------
function collectivity_metrics(A::AbstractMatrix, ε::AbstractMatrix)
    S = size(A,1)
    # build non‐dimensionalized A*
    Astar = ε .* A .- transpose(A)
    ev = eigvals(Astar)
    φ = maximum(abs, ev)
    σ = std(vec(Astar))
    # empirical connectance p between consumers (rows R+1:S) and resources (cols 1:R)
    R = size(A,1) - size(A,2) + size(A,1) - size(A,1) # workaround to recover R from context
    C = S - R
    # count nonzero consumer→resource links
    link_count = count(!iszero, A[R+1:S, 1:R])
    p = link_count / (R*C)
    C_may = σ * sqrt(p*(S-1))
    return φ, C_may
end

# ----------------------------------------------------------------
# 2) Perturbation depth: average graph‐distance of knock‐out ripples
# ----------------------------------------------------------------
function perturbation_depth(A::AbstractMatrix, p, B_eq; threshold=1e-6, T=500.)
    S = length(B_eq)
    # build unweighted graph of trophic links (only consumer→resource or resource→consumer)
    g = Graph(S)
    for i in 1:S, j in 1:S
        if A[i,j] != 0
            add_edge!(g, i, j)
        end
    end

    depths = Float64[]
    for i in 1:S
        # remove species i
        fixed = copy(B_eq); fixed[i] = 0.0
        sol = solve(ODEProblem(trophic_ode!, fixed, (0.0,T), p), Tsit5();
                    callback=build_callbacks(S, EXTINCTION_THRESHOLD),
                    abstol=1e-12, reltol=1e-12)
        Bpost = sol.u[end]
        # which species move beyond threshold
        moved = findall(j -> abs(Bpost[j] - B_eq[j]) > threshold, 1:S)
        if !isempty(moved)
            ds = distances(g, i, moved)
            push!(depths, mean(ds))
        end
    end
    return mean(depths)
end

# ----------------------------------------------------------------
# 3) Temporal unpredictability: 1 − corr(early vs. late response)
# ----------------------------------------------------------------
function temporal_unpredictability(A::AbstractMatrix, ε::AbstractMatrix, p, B_eq;
                                   δ=1e-2, t_early=50., t_total=500.)
    S = length(B_eq)
    # analytic small‐delta response vector
    Astar = ε .* A .- transpose(A)
    V = -inv(I - Astar)
    press = vcat(zeros(S - size(ε,1) + size(ε,1)), ones(size(ε,1))) # zeros for resources, ones for consumers
    ΔB_ana = V * press * δ

    # simulate uniform press on all consumers
    xi2 = copy(p[4]) .+ δ          # p[4] is xi_cons in (R,C,m_cons,xi_cons,...)
    p2 = (p[1], p[2], p[3], xi2, p[5], p[6], ε, A)
    sol = solve(ODEProblem(trophic_ode!, B_eq, (0.0,t_total), p2), Tsit5();
                saveat=[t_early, t_total],
                callback=build_callbacks(length(B_eq), EXTINCTION_THRESHOLD),
                abstol=1e-12, reltol=1e-12)
    early = sol.u[1] .- B_eq
    late  = sol.u[2] .- B_eq
    return 1 - cor(early, late)
end

# ----------------------------------------------------------------
# 4) Biotic‐niche contribution
# ----------------------------------------------------------------
function biotic_niche_contribution(B_eq, R_eq, C_eq)
    S = length(B_eq)
    K = vcat(R_eq, C_eq)
    N = B_eq ./ K
    return sum((N .- 1).^2) / sum(N.^2)
end

# ----------------------------------------------------------------
# 5) Integrate into ladder pipeline
# ----------------------------------------------------------------
function ladder_with_collectivity(df_good::DataFrame, delta_vals;
                                  ladder_steps=1:16, tspan=(0.0,500.0), t_perturb=250.0)
    results = Vector{NamedTuple}()
    for row in eachrow(df_good)[sample(1:nrow(df_good), 100, replace=false)]
        # unpack pre‐screened parameters
        R, C = row.R, row.C
        S    = R + C
        pfull = (R, C, row.m, row.xi, row.r, row.d, row.ε, row.A)
        # compute unperturbed equilibrium
        fixed = row.B_eq
        # collect baseline metrics
        φ_full, C_full = collectivity_metrics(row.A, row.ε)
        depth_full     = perturbation_depth(row.A, pfull, fixed)
        unpred_full    = temporal_unpredictability(row.A, row.ε, pfull, fixed)
        niche_full     = biotic_niche_contribution(fixed, row.R_eq, row.C_eq)

        for step in ladder_steps
            # build simplified network & ε
            A_s, ε_s = transform_for_ladder_step(step, row.A, row.ε)
            φ_s, C_s     = collectivity_metrics(A_s, ε_s)
            depth_s      = perturbation_depth(A_s, pfull, fixed)
            unpred_s     = temporal_unpredictability(A_s, ε_s, pfull, fixed)
            niche_s      = biotic_niche_contribution(fixed, row.R_eq, row.C_eq)

            push!(results, (
                S            = S,
                step         = step,
                φ            = step==1 ? φ_full : φ_s,
                C_may        = step==1 ? C_full : C_s,
                depth        = step==1 ? depth_full : depth_s,
                unpredict    = step==1 ? unpred_full : unpred_s,
                niche        = step==1 ? niche_full : niche_s,
                sens_corr    = row[Symbol("scorr_$(step==1 ? "Full" : "S$(step-1)")")]
            ))
        end
    end

    return DataFrame(results)
end

# === Example invocation ===
# df_good = ...  # from feasibility_search
df_collect = ladder_with_collectivity(df_good, [0.1]; ladder_steps=1:16)
# CSV.write("collectivity_results.csv", df_collect)
