
# ————————————————
# 2) Resilience & reactivity
# ————————————————
# Resilience ≔ –max Re(eig(J))
function compute_resilience(J::AbstractMatrix)
    λ = eigvals(J)
    return -maximum(real.(λ))
end

# Reactivity ≔ max Re(eig((J+J')/2))
function compute_reactivity(J::AbstractMatrix)
    M = (J + transpose(J)) / 2
    μ = eigvals(M)
    return maximum(real.(μ))
end

# ————————————————
# 3) Shuffle‐J routine
# ————————————————
"""
    shuffle_J!(J, step, R)

In‐place permutes entries of J according to `step`:
1 ⇒ only diag, separately for resources 1:R and consumers (R+1):S  
2 ⇒ only off‐diag, separately for resources and consumers  
3 ⇒ all entries (diag+off) but separately by groups  
4 ⇒ entire matrix at once  
"""
function shuffle_J!(J::AbstractMatrix, step::Int, R::Int)
    S = size(J,1)
    resources = 1:R
    consumers = (R+1):S

    if step == 1
        # Diagonal only, separate
        for grp in (resources, consumers)
            vals = [J[i,i] for i in grp]
            shuffle!(vals)
            for (i, v) in zip(grp, vals)
                J[i,i] = v
            end
        end

    elseif step == 2
        # Off-diagonal only, separate
        for grp in (resources, consumers)
            pos  = [(i,j) for i in grp for j in 1:S if j != i]
            vals = [J[i,j] for (i,j) in pos]
            shuffle!(vals)
            for ((i,j), v) in zip(pos, vals)
                J[i,j] = v
            end
        end

    elseif step == 3
        # All entries, separate
        for grp in (resources, consumers)
            pos  = [(i,j) for i in grp for j in 1:S]
            vals = [J[i,j] for (i,j) in pos]
            shuffle!(vals)
            for ((i,j), v) in zip(pos, vals)
                J[i,j] = v
            end
        end

    elseif step == 4
        # Diagonal only, all species together
        inds = 1:S
        vals = [J[i,i] for i in inds]
        shuffle!(vals)
        for (i, v) in zip(inds, vals)
            J[i,i] = v
        end

    elseif step == 5
        # Off-diagonal only, all species together
        pos  = [(i,j) for i in 1:S for j in 1:S if i != j]
        vals = [J[i,j] for (i,j) in pos]
        shuffle!(vals)
        for ((i,j), v) in zip(pos, vals)
            J[i,j] = v
        end

    elseif step == 6
        # All entries, full matrix shuffle
        pos  = [(i,j) for i in 1:S for j in 1:S]
        vals = [J[i,j] for (i,j) in pos]
        shuffle!(vals)
        for ((i,j), v) in zip(pos, vals)
            J[i,j] = v
        end

    else
        error("Invalid step = $step; must be 1–6")
    end
end


# ————————————————
# 4) Main experiment
# ————————————————
S = 50
R = 30
C = S - R

# Diagonal‐dominance levels: high ⇒ very stable, low ⇒ near‐critical
dominances = [0.001, 0.1, 1.0, 2.0, 5.0]
n_per_level = 100  # number of random Jacobs per level

# Prepare results table
results = DataFrame(
    dominance = Float64[],
    replicate = Int[],
    step      = Int[],   # 0 = original, 1–4 = shuffles
    resilience = Float64[],
    reactivity = Float64[],
    median_return = Float64[],
)

for (lvl, dom) in enumerate(dominances)
    for rep in 1:n_per_level
        # 4.1) Build a random Jacobian J with specified diagonal dominance
        # Off‐diagonal iid Normal(0,1):
        J = randn(S, S)
        for i in 1:S
            J[i,i] = 0.0
        end
        σ_off = mean(abs.(J))              # avg off‐diag magnitude
        dval = dom * σ_off
        for i in 1:S
            J[i,i] = -abs(rand(Normal(0.0, dval)))                # enforce diagonal dominance
        end
        for i in 1:S, j in 1:S
            if i != j && i <= R && j > R
                J[i,j] = -abs(J[i,j])
            elseif i != j && i > R && j <= R
                J[i,j] = abs(J[i,j])
            elseif i != j && i > R && j > R
                J[i,j] = 0.0
            elseif i != j && i <= R && j <= R
                J[i,j] = 0.0
            end
        end

        # 4.2) Compute original properties (step = 0)
        push!(results, (
            dominance = dom,
            replicate = rep,
            step      = 0,
            resilience   = compute_resilience(J),
            reactivity   = compute_reactivity(J),
            median_return = analytical_median_return_rate(J),
        ))

        # 4.3) Apply each shuffle‐step
        for step in 1:6
            J_s = copy(J)
            shuffle_J!(J_s, step, R)
            push!(results, (
                dominance = dom,
                replicate = rep,
                step      = step,
                resilience   = compute_resilience(J_s),
                reactivity   = compute_reactivity(J_s),
                median_return = analytical_median_return_rate(J_s),
            ))
        end
    end
end

# ————————————————
# 5) Plotting function
# ————————————————
# ————————————————
# 6) Makie‐based plotting
# ————————————————
# assume `results` exists with columns:
#   :dominance, :replicate, :step, :resilience, :reactivity, :median_return

# only look at one dominance level at a time—for example:
dom = unique(results.dominance)[1]   # pick the level you want to plot
df = filter(row -> row.dominance == dom, results)

############################################################################
############################################################################
steps = 1:3
metrics = (:resilience, :reactivity, :median_return)
dominance_levels = unique(results.dominance)

for step in steps
    fig = Figure(; size = (1000, 500))
    Label(fig[0, 1:length(dominances)], "Step $step"; fontsize=16)

    for (j, metric) in enumerate(metrics), (i, dom) in enumerate(dominance_levels)
        df = filter(row -> row.dominance == dom, results)

        ax = Axis(fig[j, i];
            title     = "Dom: $(dom), Metric: $(metric)",
            xlabel    = "Original $(metric)",
            ylabel    = "Step $step $(metric)",
            xticksvisible = true,
            yticksvisible = true,
            titlesize=9,
            xlabelsize=10, ylabelsize=10,
            xticklabelsize=9, yticklabelsize=9,
        )

        # extract original and shuffled values by replicate
        orig = df[(df.step .== 0), metric]
        shuf = df[(df.step .== step), metric]
        rep  = df[(df.step .== 0), :replicate]

        orig = orig[sortperm(rep)]
        shuf = shuf[sortperm(rep)]

        # plot
        scatter!(ax, orig, shuf; markersize = 6, color = :steelblue)

        mn = minimum(vcat(orig, shuf))
        mx = maximum(vcat(orig, shuf))
        lines!(ax, [mn, mx], [mn, mx]; linestyle = :dash, color = :black)

        # R² to 1:1 line
        ss_tot = sum((shuf .- mean(shuf)).^2)
        ss_res = sum((shuf .- orig).^2)
        r2_1to1 = 1 - ss_res / ss_tot

        text!(ax, "R²=$(round(r2_1to1, digits=3))";
            position=(mx, mn),
            align=(:right, :bottom),
            fontsize=9,
            color=:black
        )
    end

    display(fig)
end
############################################################################
############################################################################
# your definitions
step_names       = ["Diagonal", "Off-diagonal", "All"]
steps            = 1:length(step_names)
metrics          = (:resilience, :reactivity, :median_return)
metric_labels    = ("Resilience", "Reactivity", "Median-return")
dominance_levels = sort(unique(results.dominance))
dom_labels       = string.(dominance_levels)

begin
    
    fig = Figure(size = (1100, 500))

    for (r, step_name) in enumerate(step_names)
        step = steps[r]  # numeric step 1,2,3
        # build error table for this step
        errs = DataFrame(dominance=Float64[], replicate=Int[],
                        metric=Symbol[], error=Float64[])

        for dom in dominance_levels, rep in unique(results.replicate)
            orig = filter(row -> row.dominance==dom &&
                                row.replicate==rep &&
                                row.step==0, results)
            shuf = filter(row -> row.dominance==dom &&
                                row.replicate==rep &&
                                row.step==step, results)
            isempty(orig) || isempty(shuf) && continue

            for metr in metrics
                err = abs(shuf[1, metr] - orig[1, metr])
                push!(errs, (dom, rep, metr, err))
            end
        end

        for (c, metr) in enumerate(metrics)
            ax = Axis(fig[r, c];
                title  = "$(step_name) – $(metric_labels[c])",
                xlabel = "Diagonal Dominance",
                ylabel = "Absolute error",
                xticks = (1:length(dominance_levels), dom_labels),
            )

            sub = filter(row -> row.metric==metr, errs)
            positions = [ searchsortedfirst(dominance_levels, d) for d in sub.dominance ]

            boxplot!(ax, positions, sub.error)
        end
    end

    display(fig)
end

############################################################################
############################################################################
# # metrics to plot
# steps = 1:3  # make sure all 6 steps are included
# metrics = (:resilience, :reactivity, :median_return)

# for (p, dom) in enumerate(unique(results.dominance))
#     df = filter(row -> row.dominance == dom, results)
#     fig = Figure(; size = (1400, 700))
#     Label(fig[0, 1:6], "Dominance level: $dom"; fontsize=14)

#     for (j, metric) in enumerate(metrics), (i, step) in enumerate(steps)
#         ax = Axis(fig[j, i];
#             title     = "Step $step — $(metric)",
#             xlabel    = "Original $(metric)",
#             ylabel    = "Step $step $(metric)",
#             xticksvisible = true,
#             yticksvisible = true,
#             titlesize=9,
#             xlabelsize=9, ylabelsize=9,
#             xticklabelsize=9, yticklabelsize=9,
#         )

#         # extract and align by replicate
#         orig = df[(df.step .== 0), metric]
#         shuf = df[(df.step .== step), metric]
#         rep  = df[(df.step .== 0), :replicate]

#         orig = orig[sortperm(rep)]
#         shuf = shuf[sortperm(rep)]

#         scatter!(ax, orig, shuf; markersize = 6, color = :steelblue)

#         mn = minimum(vcat(orig, shuf))
#         mx = maximum(vcat(orig, shuf))
#         lines!(ax, [mn, mx], [mn, mx]; linestyle = :dash, color = :black)

#         # R² to 1:1 line
#         ss_tot = sum((shuf .- mean(shuf)).^2)
#         ss_res = sum((shuf .- orig).^2)
#         r2_1to1 = 1 - ss_res / ss_tot

#         text!(ax, "R²=$(round(r2_1to1, digits=3))";
#             position=(mx, mn),
#             align=(:right, :bottom),
#             fontsize=9,
#             color=:black
#         )
#     end

#     display(fig)
# end

