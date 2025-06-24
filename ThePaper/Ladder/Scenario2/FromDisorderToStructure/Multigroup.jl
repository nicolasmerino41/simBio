using Random, Statistics, LinearAlgebra
using DifferentialEquations
using Distributions
using QuadGK, Roots
using NLsolve
using CairoMakie

using Random, Statistics, LinearAlgebra
using DifferentialEquations
using Distributions
using QuadGK, Roots
using NLsolve
using CairoMakie

# --------------------------------------------------
# 1. Analytical L-group self-consistent φ
#    Solve for vector φ of length L:
#    φ_g = ∫ Dz [1 - Φ((∑_h μ_{g,h} φ_h + sqrt(∑_h σ_{g,h}^2 φ_h) * z - K0_g) / ζ_g )]
# --------------------------------------------------
function analytical_phi_L(K0_vec::Vector{Float64}, ζ_vec::Vector{Float64},
                          μ_mat::Matrix{Float64}, σ_mat::Matrix{Float64};
                          tol=1e-6)
    L = length(K0_vec)
    std_norm = Normal(0,1)

    function F!(Fout, φ_vec)
        # φ_vec: current guess, length L
        # Compute each Fout[g] = RHS_g - φ_g
        @inbounds for g in 1:L
            φg = clamp(φ_vec[g], 0.0, 1.0)
            # compute mean input: m_g = ∑_h μ_{g,h} φ_h
            m_g = 0.0
            v_g = 0.0
            for h in 1:L
                φh = clamp(φ_vec[h], 0.0, 1.0)
                m_g += μ_mat[g,h] * φh
                v_g += σ_mat[g,h]^2 * φh
            end
            sd_g = sqrt(max(v_g, 0.0))
            K0g = K0_vec[g]
            ζg = ζ_vec[g]
            # define integrand over z ~ N(0,1)
            integrand(z) = begin
                arg = (m_g + sd_g*z - K0g) / ζg
                (1 - cdf(std_norm, arg)) * pdf(std_norm, z)
            end
            # compute RHS_g
            RHS_g, _ = quadgk(integrand, -Inf, Inf; rtol=1e-6)
            Fout[g] = RHS_g - φg
        end
    end

    # initial guess: e.g., all 0.5
    φ0 = fill(0.5, L)
    sol = nlsolve(F!, φ0; xtol=tol, ftol=tol)
    φ_sol = sol.zero
    # clamp to [0,1]
    for g in 1:L
        φ_sol[g] = clamp(φ_sol[g], 0.0, 1.0)
    end
    return φ_sol
end

# --------------------------------------------------
# 2. Generate partial-order trophic network
#
# Inputs:
#   group_labels: vector of length S, with integer group IDs 1..L
#   connectance C (probability of potential link between any pair i<j)
#   ordering alpha ∈ [0,1]: 
#      with probability alpha, direct link respecting group order (if groups differ),
#      with probability (1-alpha), random between-group link direction.
#   If group(i)==group(j), we skip intragroup trophic link (A[i,j]=A[j,i]=0).
# Returns:
#   directed edge list: Vector{Tuple{i,j}} meaning i→j (i prey, j predator).
# --------------------------------------------------
function gen_partial_trophic_edges(group_labels::Vector{Int}, C::Float64, alpha::Float64, rng::AbstractRNG)
    S = length(group_labels)
    edges = Tuple{Int,Int}[]
    for i in 1:S-1
        gi = group_labels[i]
        for j in i+1:S
            gj = group_labels[j]
            if rand(rng) < C
                if gi != gj
                    if rand(rng) < alpha
                        # enforce direction from lower group → higher group
                        if gi < gj
                            push!(edges, (i,j))
                        else
                            push!(edges, (j,i))
                        end
                    else
                        # random direction
                        if rand(rng) < 0.5
                            push!(edges, (i,j))
                        else
                            push!(edges, (j,i))
                        end
                    end
                else
                    # same group: skip trophic link (intragroup interactions = 0 here)
                    # If you wanted intragroup competition, you could handle separately.
                end
            end
        end
    end
    return edges
end

# --------------------------------------------------
# 3. Build interaction matrix A from directed edges,
#    enforcing trophic sign structure and group-specific μ_mat, σ_mat.
#
# Inputs:
#   S, edges list from gen_partial_trophic_edges,
#   group_labels (length S, values 1..L),
#   μ_mat (L×L) mean interaction magnitude between group prey g → predator h,
#   σ_mat (L×L) sd of interaction magnitude between groups,
#   rng.
# For each directed edge i→j:
#   let g = group_labels[i], h = group_labels[j];
#   sample magnitude ~ truncated Normal(mean=μ_mat[g,h]/(size of prey group), sd=σ_mat[g,h]/sqrt(size_prey)),
#   then set A[j,i] = +magnitude, A[i,j] = -magnitude_predator_on_prey?
# Actually predator effect is negative on prey: sample magnitude for that direction 
# can use same magnitude or separate distribution? For simplicity, use same magnitude.
#
# Here we sample truncated normal for positive magnitude, then assign sign:
# --------------------------------------------------
function build_A_trophic(S::Int, edges::Vector{Tuple{Int,Int}}, 
                         group_labels::Vector{Int},
                         μ_mat::Matrix{Float64}, σ_mat::Matrix{Float64},
                         group_sizes::Vector{Int},
                         rng::AbstractRNG)
    A = zeros(Float64, S, S)
    L = length(group_sizes)
    # Pre-create truncated normal distributions for each (g,h) pair for efficiency?
    # But since group sizes differ, the denominator differs per g,h; we can sample on the fly.
    for (i,j) in edges
        g = group_labels[i]
        h = group_labels[j]
        # prey i in group g, predator j in group h
        # mean magnitude for A[j,i] is μ_mat[g,h]/(group_sizes[g]), sd = σ_mat[g,h]/√(group_sizes[g])
        mean_pos = μ_mat[g,h] / group_sizes[g]
        sd_pos = sqrt((σ_mat[g,h]^2) / group_sizes[g])
        # truncated normal ≥ 0
        # if mean_pos is many sd above zero, truncation effect small
        # but handle case mean_pos negative? Typically μ_mat[g,h]>0 for prey→pred.
        dist = truncated(Normal(mean_pos, sd_pos), 0.0, Inf)
        mag = rand(rng, dist)
        A[j, i] = mag
        # negative effect of predator on prey: use same magnitude distribution?
        # Option A: same mag (symmetric magnitude), but negative:
        A[i, j] = -mag
        # If you want different mean for predator→prey effect, 
        # you could sample from μ_mat[h,g], but trophic usually asymmetric: 
        # often only one sign matters. Here we keep magnitude symmetric for simplicity.
    end
    return A
end

# --------------------------------------------------
# 4. gLV simulation to equilibrium and measure φ per group
#    dN_i/dt = N_i (K_i - N_i - ∑_j A[i,j] N_j)
# Returns:
#   φ_g for each group g (fraction survivors in that group),
#   total φ_total.
# --------------------------------------------------
function gLV_rhs!(du, u, p, t)
    K, A = p
    Au = A * u
    @inbounds for idx in eachindex(u)
        du[idx] = u[idx] * (K[idx] - u[idx] - Au[idx])
    end
end

function simulate_phi_by_group(A::Matrix{Float64}, 
                               K0_vec::Vector{Float64}, ζ_vec::Vector{Float64},
                               group_labels::Vector{Int}, group_sizes::Vector{Int};
                               t_end=200.0, rng=GLOBAL_RNG)
    S = size(A,1)
    L = length(group_sizes)
    # Draw K_i by group: K_i ~ Normal(K0_vec[g], ζ_vec[g]^2)
    K = similar(group_labels, Float64)
    for i in 1:S
        g = group_labels[i]
        Ki = rand(rng, Normal(K0_vec[g], ζ_vec[g]))
        K[i] = max(Ki, 1e-6)
    end
    # initial u0
    u0 = rand(rng, S) .* K
    u0 .= max.(u0, 1e-6)
    prob = ODEProblem(gLV_rhs!, u0, (0.0, t_end), (K, A))
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    N_end = sol.u[end]
    N_end .= max.(N_end, 0.0)
    # survivors: N_end[i] > threshold
    thresh = 1e-3
    survivors = N_end .> thresh
    φ_g = zeros(Float64, L)
    for g in 1:L
        idxs = findall(x-> x==g, group_labels)
        if !isempty(idxs)
            φ_g[g] = sum(survivors[idxs]) / length(idxs)
        else
            φ_g[g] = 0.0
        end
    end
    φ_total = sum(survivors)/S
    return φ_g, φ_total
end

# --------------------------------------------------
# 5. Single-group prediction φ_glob
#    Collapse into one group using global moments:
#    - K0_glob = weighted mean of K0_vec
#    - ζ_glob = sqrt(weighted mean of (ζ_g^2 + K0_g^2) - K0_glob^2)
#    - μ_glob = S * mean(off-diagonal A entries)
#    - σ_glob = sqrt(S * var(off-diagonal A entries))
# --------------------------------------------------
function predicted_phi_single(A::Matrix{Float64}, 
                              K0_vec::Vector{Float64}, ζ_vec::Vector{Float64},
                              group_labels::Vector{Int})
    S = size(A,1)
    # global K0
    L = length(K0_vec)
    # weights = group_sizes / S
    # But we don't have group_sizes here; compute from labels:
    counts = zeros(Int, L)
    for g in group_labels
        counts[g] += 1
    end
    K0_glob = sum(counts[g] * K0_vec[g] for g in 1:L) / S
    # global ζ: Var(K) = E[K^2] - K0_glob^2
    E_K2 = sum(counts[g]*(ζ_vec[g]^2 + K0_vec[g]^2) for g in 1:L) / S
    varK = E_K2 - K0_glob^2
    ζ_glob = sqrt(max(varK, 1e-8))
    # global A moments
    off = Float64[]
    for i in 1:S, j in 1:S
        if i != j
            push!(off, A[i,j])
        end
    end
    meanA = mean(off)
    varA = var(off)
    μ_glob = S * meanA
    σ_glob = sqrt(S * varA)
    # predict φ_glob
    φ_glob = analytical_phi(K0_glob, ζ_glob, μ_glob, σ_glob)
    return φ_glob
end

# --------------------------------------------------
# 6. Multi-group prediction φ_vec
# --------------------------------------------------
function predicted_phi_multigroup(A::Matrix{Float64},
                                  K0_vec::Vector{Float64}, ζ_vec::Vector{Float64},
                                  group_labels::Vector{Int}, group_sizes::Vector{Int})
    # We need μ_mat and σ_mat estimated from A or from model parameters?
    # If we want to use the true μ_mat, σ_mat used in generation, pass them explicitly.
    # Alternatively, we can estimate μ_mat and σ_mat from A: 
    #   μ_mat[g,h] = S_g * mean of A[j,i] over edges i in group g, j in group h (for directed trophic networks, but A may be sparse).
    #   σ_mat[g,h]^2 = S_g * var of those A[j,i].
    # Here we assume we know μ_mat,σ_mat from generation. Better to pass them in.
    error("predicted_phi_multigroup requires μ_mat and σ_mat used to build A; call analytical_phi_L with those.")
end

# Better: we pass μ_mat,σ_mat explicitly when predicting; see usage below.

# --------------------------------------------------
# 7. Main experiment combining:
#    - partial ordering sweep over alpha_values ∈ [0,1]
#    - possibly multiple interaction strengths mu_values ∈ [...]
#    - group-specific K0_vec, ζ_vec
#    - known μ_mat, σ_mat
#    - for each combination: replicate network, simulate φ_g and φ_total,
#      predict φ_glob and φ_vec via analytical_phi_L,
#      compare errors.
# --------------------------------------------------
function experiment_partial_order_multigroup()
    rng = MersenneTwister(1234)

    # 7.1. Define groups
    # Example: L=3 (basal, intermediate, top), sizes:
    # You can generalize to more levels by changing group_sizes.
    group_sizes = [80, 80, 40]  # sum = S = 200
    L = length(group_sizes)
    S = sum(group_sizes)
    # Build group_labels: species 1..S, first group_sizes[1] in group 1, next in group 2, etc.
    group_labels = Vector{Int}(undef, S)
    idx = 1
    for g in 1:L
        for _ in 1:group_sizes[g]
            group_labels[idx] = g
            idx += 1
        end
    end

    # 7.2. Define carrying capacities per group
    # Example: distinct K0 for each group:
    K0_vec = [5.0, 3.0, 1.0]    # basal high, intermediate medium, top low
    ζ_vec  = [1.0, 1.0, 1.0]    # can vary if desired

    # 7.3. Define base interaction parameters μ_mat_base, σ_mat_base
    # μ_mat[g,h] is the baseline mean interaction magnitude for prey in group g → predator in group h.
    μ_mat_base = zeros(Float64, L, L)
    σ_mat_base = zeros(Float64, L, L)
    # For strict hierarchical: only allow interactions g→h if g < h (lower → higher).
    # Example: set a baseline mean μ0_base for adjacent levels, maybe smaller for further apart.
    # Here we choose:
    μ0_base = 5.0   # a base that will be scaled by mu_strength in outer loop
    σ0_base = 2.0   # base SD
    for g in 1:L
        for h in 1:L
            if g < h
                # e.g. same mean for any upward link; you can choose declining with h-g
                μ_mat_base[g,h] = μ0_base
                σ_mat_base[g,h] = σ0_base
            else
                μ_mat_base[g,h] = 0.0
                σ_mat_base[g,h] = 0.0
            end
        end
    end

    # 7.4. Define connectance and sweep parameters
    C = 0.1  # connectance: probability of a link between any pair i<j
    alpha_values = range(0.0, stop=1.0, length=6)  # e.g. [0,0.2,0.4,0.6,0.8,1.0]
    mu_strengths = [1.0, 2.0, 5.0, 10.0]  # scale factor for μ_mat_base; you can also sweep more finely
    # μ0 in build_A will be μ_mat_base[g,h] * scale / S_g etc.

    nrep = 10  # replicates per combination

    # Storage: results[(alpha, mu_s)] = Dict with keys :φ_sim_mean, :φ_sim_std, :φ_glob_pred, :φ_glob_err, 
    #                                            :φ_vec_pred (mean), :φ_vec_err (e.g. norm difference)
    results = Dict{Tuple{Float64,Float64}, Dict{Symbol, Any}}()

    for alpha in alpha_values, mu_s in mu_strengths
        φ_total_sims = Float64[]
        φ_glob_preds = Float64[]
        # For multi-group: track predicted φ_g for each replicate? Actually analytic prediction doesn't vary replicate-to-replicate (since μ_mat_base is fixed),
        # but because network random sampling yields small moment deviations, we can either use the theoretical μ_mat_base or estimate from each A.
        # Here we use the theoretical μ_mat_base * mu_s and σ_mat_base * mu_s, so prediction is the same across replicates.
        # Compute analytic prediction once per combination:
        # Build μ_mat = μ_mat_base * mu_s, σ_mat = σ_mat_base * mu_s
        μ_mat = μ_mat_base * mu_s
        σ_mat = σ_mat_base * mu_s
        # Solve analytic φ_vec:
        φ_vec_pred = analytical_phi_L(K0_vec, ζ_vec, μ_mat, σ_mat)
        φ_total_pred = sum(group_sizes[g]*φ_vec_pred[g] for g in 1:L) / S

        # For single-group prediction: need global A moments; but analytic single-group prediction if using theoretical global moments:
        # global mean interaction per link = (1/S) ∑_{g<h} [μ_mat[g,h] * (#possible links from g to h)/(#off-diagonals)] 
        # But simpler to compute per replicate from actual A. So we compute φ_glob_pred per replicate below.

        for rep in 1:nrep
            # 1) generate partial-order edges
            edges = gen_partial_trophic_edges(group_labels, C, alpha, rng)
            # 2) build A using μ_mat, σ_mat and group_sizes
            A = build_A_trophic(S, edges, group_labels, μ_mat, σ_mat, group_sizes, rng)
            # 3) simulate to get φ_g_sim, φ_total_sim
            φ_g_sim, φ_total_sim = simulate_phi_by_group(A, K0_vec, ζ_vec, group_labels, group_sizes; rng=rng)
            push!(φ_total_sims, φ_total_sim)
            # 4) single-group prediction from this A
            φ_glob = predicted_phi_single(A, K0_vec, ζ_vec, group_labels)
            push!(φ_glob_preds, φ_glob)
            # note: multi-group prediction is φ_vec_pred, same each replicate (given theoretical μ_mat).
        end

        # aggregate
        mean_sim = mean(φ_total_sims)
        std_sim = std(φ_total_sims)
        mean_glob_pred = mean(φ_glob_preds)
        std_glob_pred = std(φ_glob_preds)
        # record difference between single-group predicted and simulated
        err_glob = mean(abs.(φ_glob_preds .- φ_total_sims))

        # multi-group error: compare φ_vec_pred to simulated φ_g_sim? We need φ_g_sim per replicate:
        # To properly measure multi-group error, we'd need to record φ_g_sim per replicate. Let's do that above.

        # For simplicity, here we compute multi-group error against mean simulated φ:
        # First compute mean simulated φ_g across replicates:
        φ_g_sims_all = zeros(Float64, L, nrep)
        # We need to collect φ_g_sim each rep. Let's modify above to collect φ_g_sim per replicate:
        # Instead of pushing only φ_total_sims, also collect φ_g_sim in a matrix.
        # To keep code clear, we restructure slightly:

        # We redo the replicate loop with storage for φ_g:
        φ_total_sims2 = Float64[]
        φ_glob_preds2 = Float64[]
        φ_g_sims_mat = zeros(Float64, L, nrep)
        for rep in 1:nrep
            edges = gen_partial_trophic_edges(group_labels, C, alpha, rng)
            A = build_A_trophic(S, edges, group_labels, μ_mat, σ_mat, group_sizes, rng)
            φ_g_sim, φ_total_sim = simulate_phi_by_group(A, K0_vec, ζ_vec, group_labels, group_sizes; rng=rng)
            φ_glob = predicted_phi_single(A, K0_vec, ζ_vec, group_labels)
            # store
            φ_total_sims2_rep = φ_total_sim
            φ_total_sims2_push = push!(φ_total_sims2, φ_total_sim)
            push!(φ_glob_preds2, φ_glob)
            for g in 1:L
                φ_g_sims_mat[g, rep] = φ_g_sim[g]
            end
        end
        # recompute aggregates
        mean_sim = mean(φ_total_sims2)
        std_sim = std(φ_total_sims2)
        err_glob = mean(abs.(φ_glob_preds2 .- φ_total_sims2))
        mean_glob_pred = mean(φ_glob_preds2)
        std_glob_pred = std(φ_glob_preds2)
        # multi-group error: compare φ_vec_pred[g] to mean simulated φ_g
        mean_φ_g_sim = [mean(φ_g_sims_mat[g, :]) for g in 1:L]
        err_multi = mean(abs.(φ_vec_pred .- mean_φ_g_sim))

        # store results
        results[(alpha, mu_s)] = Dict(
            :mean_sim => mean_sim,
            :std_sim => std_sim,
            :mean_glob_pred => mean_glob_pred,
            :std_glob_pred => std_glob_pred,
            :err_glob => err_glob,
            :φ_vec_pred => φ_vec_pred,
            :mean_φ_g_sim => mean_φ_g_sim,
            :err_multi => err_multi
        )
    end

    # --------------------------------------------------
    # 8. Plotting results
    #    Example plots:
    #    (a) For each mu_s, plot error_glob and error_multi vs alpha.
    #    (b) For a chosen alpha, plot φ_total_sim ± std vs mu_s with predictions overlaid.
    #    (c) For multi-group: plot mean simulated φ_g vs φ_vec_pred for each group.
    # --------------------------------------------------

    # (a) error vs ordering alpha
    fig1 = Figure(resolution=(800,400))
    ax1 = Axis(fig1[1,1]; xlabel="Ordering α", ylabel="Mean |φ_pred - φ_sim|",
               title="Prediction error vs ordering α for different μ scales")
    colors = distinguishable_colors(length(mu_strengths))
    for (idx, mu_s) in enumerate(mu_strengths)
        errs_glob = Float64[]
        errs_multi = Float64[]
        for alpha in alpha_values
            res = results[(alpha, mu_s)]
            push!(errs_glob, res[:err_glob])
            push!(errs_multi, res[:err_multi])
        end
        # plot single-group error
        lines!(ax1, alpha_values, errs_glob; color=colors[idx], linestyle=:dash, label="glob μscale=$(mu_s)")
        # plot multi-group error
        lines!(ax1, alpha_values, errs_multi; color=colors[idx], linestyle=:solid, label="multi μscale=$(mu_s)")
    end
    axislegend(ax1; position=:rt)
    display(fig1)

    # (b) For a chosen alpha (e.g., alpha=1.0 or intermediate), plot φ_total vs mu_s
    chosen_alpha = 0.5
    fig2 = Figure(resolution=(800,400))
    ax2 = Axis(fig2[1,1]; xlabel="Interaction scale μ_s", ylabel="Surviving fraction φ_total",
               title="Simulated and predicted φ_total vs μ scale at α=$(chosen_alpha)")
    sim_means = Float64[]
    sim_stds = Float64[]
    glob_preds = Float64[]
    multi_preds = Float64[]
    for mu_s in mu_strengths
        res = results[(chosen_alpha, mu_s)]
        push!(sim_means, res[:mean_sim])
        push!(sim_stds, res[:std_sim])
        push!(glob_preds, res[:mean_glob_pred])
        push!(multi_preds, sum(group_sizes[g]*res[:φ_vec_pred][g] for g in 1:L)/S)
    end
    # error bars for simulated
    for i in 1:length(mu_strengths)
        x = mu_strengths[i]
        y = sim_means[i]; yerr = sim_stds[i]
        lines!(ax2, [x,x], [y-yerr, y+yerr]; color=:gray, linestyle=:dash)
    end
    scatter!(ax2, mu_strengths, sim_means; color=:black, marker=:circle, label="sim")
    lines!(ax2, mu_strengths, glob_preds; color=:blue, linestyle=:dash, linewidth=2, label="single-group pred")
    lines!(ax2, mu_strengths, multi_preds; color=:green, linestyle=:solid, linewidth=2, label="multi-group pred")
    axislegend(ax2; position=:rt)
    display(fig2)

    # (c) For each combination, one could also plot group-level φ:
    # For example at chosen_alpha and chosen mu_s:
    chosen_mu = mu_strengths[2]  # pick second
    res = results[(chosen_alpha, chosen_mu)]
    mean_φ_g = res[:mean_φ_g_sim]
    φ_g_pred = res[:φ_vec_pred]
    fig3 = Figure(resolution=(400,300))
    ax3 = Axis(fig3[1,1]; xlabel="Group index g", ylabel="φ_g",
               title="Group-level φ: sim vs pred at α=$(chosen_alpha), μscale=$(chosen_mu)",
               xticks=(1:L, ["G$(g)" for g in 1:L]))
    scatter!(ax3, 1:L, mean_φ_g; color=:black, marker=:circle, label="sim")
    lines!(ax3, 1:L, φ_g_pred; color=:red, linestyle=:solid, marker=:diamond, label="pred")
    axislegend(ax3; position=:rt)
    display(fig3)

    return results
end

# Run the experiment
results = experiment_partial_order_multigroup()
