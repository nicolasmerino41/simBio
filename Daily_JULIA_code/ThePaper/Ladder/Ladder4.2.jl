using DifferentialEquations
using Random
using Statistics
using DataFrames
using Logging
using Graphs

# -----------------------------------------------------------
# Helper: pick up to max_steps integer R values between minR and maxR
# -----------------------------------------------------------
function pick_steps(minR::Int, maxR::Int, max_steps::Int)
    allR = collect(minR:maxR)
    if length(allR) <= max_steps
        return allR
    else
        picked = unique(round.(Int, LinRange(minR, maxR, max_steps)))
        push!(picked, minR, maxR)
        return sort(unique(picked))
    end
end

# -----------------------------------------------------------
# Main simulation pipeline with structural heterogeneity
# -----------------------------------------------------------
function rppp(
    S_vals, conn_vals, delta_vals;
    max_RC_steps   = 10,
    ladder_steps   = 1:16,
    Niter          = 10,
    tspan          = (0.0, 50.0),
    t_perturb      = 25.0,
    plot_full      = false,
    plot_simple    = false,
    d_value        = 1.0,
    IS_reduction   = 1.0,
    μ_R            = 1.0,
    σ_R            = 1.0,
    μ_C            = 0.1,
    σ_C            = 0.2,
    epsilon_mean   = 0.2,
    epsilon_sd     = 0.2,
    network_model  = :erdos_renyi
)
    results = Vector{NamedTuple}()
    results_lock = ReentrantLock()

    suffix(step) = step == 1 ? "Full" : "S$(step-1)"

    for S in S_vals
        for conn in conn_vals
            for delta in delta_vals
                for R in pick_steps(5, S-1, max_RC_steps)
                    C = S - R
                    total = S

                    for iter in 1:Niter
                        try
                            # --- 1) Build interaction matrix A_adj ---
                            A_adj = zeros(total, total)
                            if network_model == :erdos_renyi
                                for i in R+1:total, j in 1:total
                                    if i != j && rand() < conn
                                        w = rand() * IS_reduction
                                        A_adj[i,j] =  w
                                        A_adj[j,i] = -w
                                    end
                                end
                            else
                                error("Unsupported network_model: $network_model")
                            end

                            # --- 2) Sample skewed SAD (pyramidal) ---
                            R_eq = abs.(rand(LogNormal(log(μ_R) - σ_R^2/2, σ_R), R))
                            C_eq = abs.(rand(LogNormal(log(μ_C) - σ_C^2/2, σ_C), C))
                            fixed = vcat(R_eq, C_eq)

                            # --- 3) Core parameters & efficiencies ---
                            m_cons = fill(0.1, C)
                            d_res   = fill(d_value, R)
                            ε_full  = clamp.(rand(LogNormal(epsilon_mean, epsilon_sd), total, total), 0, 1)

                            # --- 4) Calibrate equilibrium ---
                            p_cal = (R, C, m_cons, d_res, ε_full, A_adj)
                            xi_cons, r_res = calibrate_params(R_eq, C_eq, p_cal)
                            if any(isnan, xi_cons) || any(isnan, r_res)
                                continue
                            end

                            # --- 5) Structural heterogeneity metrics ---
                            # Undirected graph for degree & community detection
                            g_undir = SimpleGraph(A_adj .!= 0)
                            deg_cons       = degree(g_undir)[R+1:end]
                            degree_cv      = std(deg_cons) / mean(deg_cons)

                            # Directed graph for trophic levels
                            g_dir = SimpleDiGraph(transpose(A_adj .> 0))
                            # multi-source BFS from all resources
                            levels = fill(typemax(Int), total)
                            queue = collect(1:R)
                            for v in 1:R
                                levels[v] = 0
                            end
                            while !isempty(queue)
                                v = popfirst!(queue)
                                for nbr in outneighbors(g_dir, v)
                                    if levels[nbr] > levels[v] + 1
                                        levels[nbr] = levels[v] + 1
                                        push!(queue, nbr)
                                    end
                                end
                            end
                            trophic_var  = var(levels[R+1:end])

                            # Community detection & modularity (label-propagation)
                            membership     = label_propagation_communities(g_undir)
                            modularity_sc  = modularity(g_undir, membership)

                            pyramid_skew = σ_R / σ_C

                            # --- 6) Full model Jacobian & sensitivity ---
                            p_full        = (R, C, m_cons, xi_cons, r_res, d_res, ε_full, A_adj)
                            J_full        = build_jacobian(fixed, p_full)
                            V_analytical  = compute_analytical_V(J_full, R, C, m_cons, xi_cons)

                            sol = solve(ODEProblem(trophic_ode!, fixed, tspan, p_full),
                                        Tsit5(); reltol=1e-8, abstol=1e-8)
                            if sol.t[end] < t_perturb || any(!isfinite, sol.u[end])
                                continue
                            end
                            Beq = sol.u[end]

                            # --- 7) Press perturbation & species metrics ---
                            rt_full, os_full, ire_full, before_p, B2_full =
                                simulate_press_perturbation(
                                    fixed, p_full, tspan, t_perturb, delta;
                                    solver=Tsit5(),
                                    plot=plot_full,
                                    show_warnings=true,
                                    full_or_simple=true
                                )
                            pred_resp_vec = vec(sum(V_analytical, dims=2))
                            obs_resp_vec  = vec((B2_full .- fixed) ./ delta)
                            sens_corr_full = cor(pred_resp_vec, obs_resp_vec)
                            pred_sens_vec  = vec(sum(abs.(V_analytical), dims=1))

                            # --- 8) Collect 'Full' metrics ---
                            metrics = Dict{String,NamedTuple}()
                            metrics["Full"] = (
                                persistence      = mean(Beq .> EXTINCTION_THRESHOLD),
                                return_time      = rt_full,
                                overshoot        = os_full,
                                ire              = ire_full,
                                resilience       = compute_resilience(Beq, p_full),
                                reactivity       = compute_reactivity(Beq, p_full),
                                sens_corr        = sens_corr_full,
                                pred_sensitivity = pred_sens_vec
                            )

                            # --- 9) Simplified steps ---
                            for step in ladder_steps[2:end]
                                suf = suffix(step)
                                A_s, ε_s = transform_for_ladder_step(step, A_adj, ε_full)
                                p_s      = (R, C, m_cons, xi_cons, r_res, d_res, ε_s, A_s)
                                J_s      = build_jacobian(fixed, p_s)
                                V_s      = compute_analytical_V(J_s, R, C, m_cons, xi_cons)

                                sol2 = solve(ODEProblem(trophic_ode!, fixed, tspan, p_s),
                                            Tsit5(); reltol=1e-8, abstol=1e-8)
                                if sol2.t[end] < t_perturb || any(!isfinite, sol2.u[end])
                                    metrics[suf] = NamedTuple()
                                    continue
                                end
                                Beq2 = sol2.u[end]

                                rt_s, os_s, ire_s, _, B2_s =
                                    simulate_press_perturbation(
                                        fixed, p_s, tspan, t_perturb, delta;
                                        solver=Tsit5(),
                                        plot=plot_simple,
                                        show_warnings=true,
                                        full_or_simple=false
                                    )
                                corr_s     = cor(vec(sum(V_s, dims=2)), vec((B2_s .- fixed) ./ delta))
                                pred_s_vec = vec(sum(abs.(V_s), dims=1))

                                metrics[suf] = (
                                    persistence      = mean(Beq2 .> EXTINCTION_THRESHOLD),
                                    return_time      = rt_s,
                                    overshoot        = os_s,
                                    ire              = ire_s,
                                    resilience       = compute_resilience(Beq2, p_s),
                                    reactivity       = compute_reactivity(Beq2, p_s),
                                    sens_corr        = corr_s,
                                    pred_sensitivity = pred_s_vec
                                )
                            end

                            # --- 10) Flatten results ---
                            base = (
                                species_count  = total,
                                resource_count = R,
                                consumer_count = C,
                                connectance    = conn,
                                perturb_delta  = delta,
                                iteration      = iter,
                                degree_CV      = degree_cv,
                                trophic_var    = trophic_var,
                                modularity     = modularity_sc,
                                pyramid_skew   = pyramid_skew
                            )

                            rec = base
                            for step in ladder_steps
                                suf = suffix(step)
                                m   = metrics[suf]
                                rec = merge(rec,
                                    (; Symbol("persistence_$suf")      => m.persistence,
                                       Symbol("return_time_$suf")      => m.return_time,
                                       Symbol("overshoot_$suf")        => m.overshoot,
                                       Symbol("ire_$suf")              => m.ire,
                                       Symbol("resilience_$suf")       => m.resilience,
                                       Symbol("reactivity_$suf")       => m.reactivity,
                                       Symbol("senscorr_$suf")         => m.sens_corr,
                                       Symbol("predsens_$suf")         => m.pred_sensitivity)
                                )
                            end

                            lock(results_lock) do
                                push!(results, rec)
                            end

                        catch e
                            @warn "Iteration error" S R C conn delta iter e
                        end
                    end
                end
            end
        end
    end

    return DataFrame(results)
end

# Example invocation
begin
    S_vals          = [10]
    connectance_list = 0.1:0.1:1.0
    delta_list      = 0.1:0.1:0.1

    new_results_7 = rppp(
        S_vals, connectance_list, delta_list;
        max_RC_steps   = 10,
        ladder_steps   = 1:16,
        Niter          = 10,
        tspan          = (0.0, 50.0),
        t_perturb      = 25.0,
        plot_full      = false,
        plot_simple    = false,
        d_value        = 1.0,
        IS_reduction   = 1.0,
        μ_R            = 1.0,
        σ_R            = 1.0,
        μ_C            = 0.1,
        σ_C            = 0.2,
        epsilon_mean   = 0.2,
        epsilon_sd     = 0.2,
        network_model  = :erdos_renyi
    )
end
