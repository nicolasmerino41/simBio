left = new_results3[:, 1:3]
right = new_results3[:, 4:end]
left.total_species = left[:, 2] + left[:, 3]
left.CR_ratio = left[:, 3] ./ left[:, 2]
new_results3 = hcat(left, right)

# (A) Return time
m_rt = lm(@formula(return_time_Full ~ total_species + CR_ratio + connectance), new_results3)
println(coeftable(m_rt))

# (B) Resilience
m_res = lm(@formula(resilience_Full ~ total_species + CR_ratio + connectance), new_results3)
println(coeftable(m_res))

# (C) Reactivity
m_rea = lm(@formula(reactivity_Full ~ total_species + CR_ratio + connectance), new_results3)
println(coeftable(m_rea))

rs = cor(new_results3.return_time_Full, new_results3.total_species)
println("r(return_time, total_species) = ", rs)

rs = cor(new_results3.resilience_Full, new_results3.connectance)
println("r(resilience, connectance)     = ", rs)

# …and so on for each pair you care about…
using RegressionTables
vif(m_rt)

###############################################################
###############################################################
# 1) Prepare an empty “long” DataFrame with the right column types
long = DataFrame(
    step    = String[],
    ΔRT     = Float64[],
    δRT     = Float64[],
    ΔRes    = Float64[],
    δRes    = Float64[],
    ΔReac   = Float64[],
    δReac   = Float64[],
)

# 2) Melt in each step’s errors
step_keys = ["S$(i)" for i in 1:15]   # adjust to your actual number of steps
df = new_results3

for s in step_keys
    # pull out full vs step metrics
    full_rt = df.return_time_Full
    rt_s    = df[!, Symbol("return_time_$s")]

    full_res= df.resilience_Full
    res_s   = df[!, Symbol("resilience_$s")]

    full_rea= df.reactivity_Full
    rea_s   = df[!, Symbol("reactivity_$s")]

    # compute absolute & relative deviations
    ΔRT   = rt_s .- full_rt
    δRT   = ΔRT  ./ full_rt

    ΔRes  = res_s .- full_res
    δRes  = ΔRes ./ full_res

    ΔReac = rea_s .- full_rea
    δReac = ΔReac./ full_rea

    n = length(ΔRT)

    # build a small DataFrame *with keyword args only*
    df_temp = DataFrame(
      step   = fill(s, n),
      ΔRT    = ΔRT,
      δRT    = δRT,
      ΔRes   = ΔRes,
      δRes   = δRes,
      ΔReac  = ΔReac,
      δReac  = δReac,
    )
    append!(long, df_temp)
end

# 3) Now you can summarize by step:
summary = DF.combine(groupby(long, :step),
  :δRT   => mean => :mean_relRT,
  :δRes  => mean => :mean_relRes,
  :δReac => mean => :mean_relReac,
)

# 3) Plot mean relative‐error vs. step
begin
    fig = Figure(resolution=(1400,500))
    ax1 = Axis(fig[1,1]; xlabel="Step", ylabel="Mean δRT",  xticks=(1:length(step_keys), step_keys))
    lines!(ax1, 1:length(step_keys), summary.mean_relRT)

    ax2 = Axis(fig[1,2]; xlabel="Step", ylabel="Mean δRes", xticks=(1:length(step_keys), step_keys))
    lines!(ax2, 1:length(step_keys), summary.mean_relRes)

    ax3 = Axis(fig[1,3]; xlabel="Step", ylabel="Mean δReac",xticks=(1:length(step_keys), step_keys))
    lines!(ax3, 1:length(step_keys), summary.mean_relReac)

    display(fig)
end

# 4) Boxplots of δRT by step
begin
    cat = map(x->findfirst(isequal(x), step_keys), long.step) |> x->x./maximum(x) .* 15 .+ 1
    fig2 = Figure(resolution=(800,400))
    ax = Axis(fig2[1,1]; xlabel="Step", ylabel="δRT")
    MK.boxplot!(ax, cat, long.δRT)
    display(fig2)
end

# 5) One‐way ANOVA / OLS to test which steps differ
#    Here as an example for relative RT:
model = lm(@formula(δRT ~ 1 + step), long)
println(coeftable(model))


long.step_idx = map(x->findfirst(isequal(x), step_keys), long.step) |> x->x./maximum(x) .* 15 .+ 1
model1 = lm(@formula(δRT ~ 1 + step_idx), long)
println(coeftable(model1))

using StatsModels
model2 = lm(@formula(δRT ~ 1 + step_idx + step_idx^2), long)
println(coeftable(model2))

##################################################################
# PLOTTING THE FIT OF THE MODELS
# 1) Define your exact step order
step_keys = ["S$(i)" for i in 1:15]    # adjust to however many steps you actually have

# 2) Summarize mean δRT by step
summary = DF.combine(groupby(long, :step),
    :δRT => mean => :mean_relRT
)
# Add numeric step index
summary.step_idx = map(s->findfirst(==(s), step_keys), summary.step)

# 3) Compute model‐predicted δRT at each step
coefs1 = coef(model1)  # [Intercept, step_idx]
lin_pred  = coefs1[1] .+ coefs1[2] .* summary.step_idx

coefs2 = coef(model2)  # [Intercept, step_idx, step_idx^2]
quad_pred = coefs2[1] .+
            coefs2[2] .* summary.step_idx .+
            coefs2[3] .* (summary.step_idx .^ 2)

# 4) Plot data + fits
begin
    fig = Figure(resolution = (800, 400))
    ax  = Axis(fig[1, 1];
        xlabel = "Step index",
        ylabel = "Mean rel. RT error (δRT)",
        title  = "Data vs. Linear and Quadratic Fit"
    )

    # scatter the raw means
    scatter!(
        ax,
        summary.step_idx,
        summary.mean_relRT;
        color      = :blue,
        markersize = 8,
        label      = "Data"
    )

    # linear fit
    lines!(
        ax,
        summary.step_idx,
        lin_pred;
        color      = :red,
        linewidth  = 2,
        label      = "Linear fit"
    )

    # quadratic fit
    lines!(
        ax,
        summary.step_idx,
        quad_pred;
        color      = :green,
        linewidth  = 2,
        label      = "Quadratic fit"
    )

    axislegend(ax; position = :rt, tellwidth = false)
    display(fig)
end