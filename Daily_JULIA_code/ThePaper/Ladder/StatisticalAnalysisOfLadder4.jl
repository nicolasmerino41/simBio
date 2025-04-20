############################################################
# 0) Load all required packages
############################################################
using DataFrames
using CategoricalArrays
using Statistics
using GLM
using RegressionTables
using CairoMakie         # for plotting
using StatsModels

############################################################
# 1) Build & clean the “wide” DataFrame
############################################################

# — your raw results, split into “left” (inputs) and “right” (metrics)
left  = new_results4[:, 1:4]
right = new_results4[:, 5:end]

# — compute consumer:total species ratio
left.C_ratio = left.consumer_count ./ left.species_count

# — concatenate into one DataFrame
df = hcat(left, right)

# — convert any NaN in _any_ Float column to `missing`
for c in names(df)
    if eltype(df[!, c]) <: AbstractFloat
        df[!, c] = ifelse.(isnan.(df[!, c]), missing, df[!, c])
    end
end

# — drop rows where we’re missing any of the vars needed for the first regressions
required_cols = [
    :return_time_Full, :resilience_Full, :reactivity_Full,
    :species_count,   :C_ratio,         :connectance
]
df = dropmissing(df)
df = filter(row -> !(row.return_time_Full == 0.0 || row.resilience_Full == 0.0), df)
############################################################
# 2) Fit & display the three “Full”‐step linear models
############################################################
# (A) Return‐time ∼ species_count + C_ratio + connectance
m_rt = lm(@formula(return_time_Full ~ species_count + C_ratio + connectance), df)

# (B) Resilience  ∼ species_count + C_ratio + connectance
m_res = lm(@formula(resilience_Full   ~ species_count + C_ratio + connectance), df)

# (C) Reactivity  ∼ species_count + C_ratio + connectance
m_rea = lm(@formula(reactivity_Full   ~ species_count + C_ratio + connectance), df)

println("\n=== Return‐time model ===")
println(coeftable(m_rt))

println("\n=== Resilience model ===")
println(coeftable(m_res))

println("\n=== Reactivity model ===")
println(coeftable(m_rea))

# — pairwise correlations (now with no missing values)
println("\nr(return_time, species_count) = ",
        cor(df.return_time_Full, df.species_count))
println("r(resilience, connectance)     = ",
        cor(df.resilience_Full, df.connectance))

# — multicollinearity check on the RT model
println("\nVariance inflation factors for the RT model:")
println(vif(m_rt))

############################################################
# 3) Melt into a “long” table of simplification‐errors
############################################################
# — prepare empty “long” DataFrame
long = DataFrame(
    step   = String[],
    ΔRT    = Float64[], δRT   = Float64[],
    ΔRes   = Float64[], δRes  = Float64[],
    ΔReac  = Float64[], δReac = Float64[],
)

# — define your simplification‐step labels
step_keys = ["S$(i)" for i in 1:15]   # adjust 15→actual number of steps

# — loop through each step, compute absolute (Δ) & relative (δ) deviations
for s in step_keys
    full_rt   = df.return_time_Full
    rt_s      = df[!, Symbol("return_time_$s")]

    full_res  = df.resilience_Full
    res_s     = df[!, Symbol("resilience_$s")]

    full_rea  = df.reactivity_Full
    rea_s     = df[!, Symbol("reactivity_$s")]

    # Δ = step – full,  δ = Δ/full
    ΔRT   = rt_s   .- full_rt
    δRT   = ΔRT    ./ full_rt

    ΔRes  = res_s  .- full_res
    δRes  = ΔRes   ./ full_res

    ΔReac = rea_s  .- full_rea
    δReac = ΔReac  ./ full_rea

    n = length(ΔRT)
    df_temp = DataFrame(
        step   = fill(s, n),
        ΔRT    = ΔRT,    δRT   = δRT,
        ΔRes   = ΔRes,   δRes  = δRes,
        ΔReac  = ΔReac,  δReac = δReac,
    )
    append!(long, df_temp)
end

# — again convert any remaining NaNs to missing
for c in [:ΔRT, :δRT, :ΔRes, :δRes, :ΔReac, :δReac]
    long[!, c] = ifelse.(isnan.(long[!, c]), missing, long[!, c])
end

# — drop any rows missing the key error‐metrics
long = dropmissing(long)
############################################################
# 4) Summarize & plot mean relative‐errors by step
############################################################
# — convert `step` into a true factor (categorical) and drop unused levels
long.step = categorical(long.step)
# droplevels! works on the CategoricalArray itself:
droplevels!(long.step)

# — compute mean relative‐error for each step
summary = combine(groupby(long, :step),
    :δRT   => mean => :mean_relRT,
    :δRes  => mean => :mean_relRes,
    :δReac => mean => :mean_relReac,
)

# — for plotting, build a numeric index 1…n_steps
n_steps = length(step_keys)
summary.step_idx = map(s->findfirst(==(s), step_keys), summary.step)

# — line plots of mean δRT, δRes, δReac vs. step
begin
    fig = Figure(; size = (1400,500))
    ax1 = Axis(fig[1,1], xlabel="Step", ylabel="Mean δRT",
            xticks=(1:n_steps, step_keys))
    lines!(ax1, summary.step_idx, summary.mean_relRT)

    ax2 = Axis(fig[1,2], xlabel="Step", ylabel="Mean δRes",
            xticks=(1:n_steps, step_keys))
    lines!(ax2, summary.step_idx, summary.mean_relRes)

    ax3 = Axis(fig[1,3], xlabel="Step", ylabel="Mean δReac",
            xticks=(1:n_steps, step_keys))
    lines!(ax3, summary.step_idx, summary.mean_relReac)

    display(fig)
end

# — boxplot of δRT by step
begin
    cats = [ idx_map[s] for s in long.step ]
    fig2 = Figure(; size = (800,400))
    ax = Axis(
        fig2[1,1], xlabel="Step", ylabel="δRT"
        )
    boxplot!(ax, cats, long.δRT)

    # — boxplot of δRes by step
    fig3 = Figure(; size = (800,400))
    ax = Axis(
        fig2[1,2], xlabel="Step", ylabel="δRes"
    )
    boxplot!(ax, cats, long.δRes)

    # — boxplot of δReac by step
    fig4 = Figure(; size = (800,400))
    ax = Axis(
        fig2[1,3], xlabel="Step", ylabel="δReac"
    )
    boxplot!(ax, cats, long.δReac)
    display(fig2)
end

############################################################
# 5) ANOVA + linear/quadratic fits on δRT vs. step
############################################################
# (A) ANOVA‐style: separate mean for each step
anova_rt = lm(@formula(δRT ~ 1 + step), long)
println("\nANOVA‐style OLS (δRT ~ step):")
println(coeftable(anova_rt))

# (B) Linear in numeric step index
long.step_idx = cats
lin_rt   = lm(@formula(δRT ~ 1 + step_idx), long)
println("\nLinear OLS (δRT ~ step_idx):")
println(coeftable(lin_rt))

# (C) Quadratic
quad_rt  = lm(@formula(δRT ~ 1 + step_idx + step_idx^2), long)
println("\nQuadratic OLS (δRT ~ step_idx + step_idx^2):")
println(coeftable(quad_rt))

# — overlay fits on the mean‐error scatter
lin_pred  = predict(lin_rt,  summary)
quad_pred = predict(quad_rt, summary)

begin
    # assume `clean` has only two columns: :step_idx and :mean_relRT
    clean = sort(summary, :step_idx)
    x = Float64.(clean.step_idx)
    y = Float64.(clean.mean_relRT)
    lin_pred  = coefs1[1] .+ coefs1[2] .* x
    quad_pred = coefs2[1] .+ coefs2[2] .* x .+ coefs2[3] .* x.^2

    fig = Figure(; size=(800,400))
    ax  = Axis(fig[1,1], xlabel="Step index", ylabel="Mean δRT",
            title="Data vs. Linear & Quadratic Fit")

    scatter!(ax, x, y, label="Data")
    lines!(ax, x, lin_pred,  linewidth=2, label="Linear fit")
    lines!(ax, x, quad_pred, linewidth=2, label="Quadratic fit")

    axislegend(ax)
    display(fig)
end

############################################################
# 6) ANOVA + linear/quadratic fits on δRes vs. step
############################################################
# (A) ANOVA‐style: separate mean for each step
anova_res = lm(@formula(δRes ~ 1 + step), long)
println("\nANOVA‐style OLS (δRes ~ step):")
println(coeftable(anova_res))

# (B) Linear in numeric step index
lin_res   = lm(@formula(δRes ~ 1 + step_idx), long)
println("\nLinear OLS (δRes ~ step_idx):")
println(coeftable(lin_res))

# (C) Quadratic
quad_res  = lm(@formula(δRes ~ 1 + step_idx + step_idx^2), long)
println("\nQuadratic OLS (δRes ~ step_idx + step_idx^2):")
println(coeftable(quad_res))

# — overlay fits on the mean‐error scatter
lin_pred  = predict(lin_res,  summary)
quad_pred = predict(quad_res, summary)

begin
    clean = sort(summary, :step_idx)
    x = Float64.(clean.step_idx)
    y = Float64.(clean.mean_relRes)
    lin_pred  = coefs1[1] .+ coefs1[2] .* x
    quad_pred = coefs2[1] .+ coefs2[2] .* x .+ coefs2[3] .* x.^2

    fig = Figure(; size=(800,400))
    ax  = Axis(fig[1,1], xlabel="Step index", ylabel="Mean δRes",
            title="Data vs. Linear & Quadratic Fit")

    scatter!(ax, x, y, label="Data")
    lines!(ax, x, lin_pred,  linewidth=2, label="Linear fit")
    lines!(ax, x, quad_pred, linewidth=2, label="Quadratic fit")

    axislegend(ax)
    display(fig)
end

############################################################
# 7) ANOVA + linear/quadratic fits on δReac vs. step
############################################################
# (A) ANOVA‐style: separate mean for each step
anova_reac = lm(@formula(δReac ~ 1 + step), long)
println("\nANOVA‐style OLS (δReac ~ step):")
println(coeftable(anova_reac))

# (B) Linear in numeric step index
lin_reac   = lm(@formula(δReac ~ 1 + step_idx), long)
println("\nLinear OLS (δReac ~ step_idx):")
println(coeftable(lin_reac))

# (C) Quadratic
quad_reac  = lm(@formula(δReac ~ 1 + step_idx + step_idx^2), long)
println("\nQuadratic OLS (δReac ~ step_idx + step_idx^2):")
println(coeftable(quad_reac))

# — overlay fits on the mean‐error scatter
lin_pred  = predict(lin_reac,  summary)
quad_pred = predict(quad_reac, summary)

begin
    clean = sort(summary, :step_idx)
    x = Float64.(clean.step_idx)
    y = Float64.(clean.mean_relReac)
    lin_pred  = coefs1[1] .+ coefs1[2] .* x
    quad_pred = coefs2[1] .+ coefs2[2] .* x .+ coefs2[3] .* x.^2

    fig = Figure(; size=(800,400))
    ax  = Axis(fig[1,1], xlabel="Step index", ylabel="Mean δReac",
            title="Data vs. Linear & Quadratic Fit")

    scatter!(ax, x, y, label="Data")
    lines!(ax, x, lin_pred,  linewidth=2, label="Linear fit")
    lines!(ax, x, quad_pred, linewidth=2, label="Quadratic fit")

    axislegend(ax)
    display(fig)
end
