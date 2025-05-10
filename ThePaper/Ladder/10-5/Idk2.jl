
using CairoMakie, Statistics

# 1) compute ratio and base‐mask
res   = A.resilience_full
ratio = A.rt_pulse_full ./ A.rt_pulse_S6
mask0 = .!ismissing.(res) .& .!ismissing.(ratio) .&
        .!isnan.(res)      .& .!isnan.(ratio)      .&
        isfinite.(res)     .& isfinite.(ratio)

# 2) find 95th percentile cutoff on the “good” subset
pct95 = quantile(ratio[mask0], 0.95)

# 3) build final mask excluding the top 5% of ratios
mask = mask0 .& (ratio .<= pct95)

# 4) recompute correlation on filtered data
r_val = cor(res[mask], ratio[mask])

# 5) plot without the big outliers
fig = Figure(resolution = (600,400))
ax  = Axis(fig[1,1];
    xlabel = "Resilience (–max Re(λ))",
    ylabel = "rt_pulse_full / rt_pulse_S6",
    title  = "Resilience vs RT‐ratio (Step 6) — 95% quantile"
)

scatter!(ax,
    res[mask], ratio[mask];
    color      = :blue,
    markersize = 6,
    alpha      = 0.7
)

# reference line at y=1
lines!(ax,
    [minimum(res[mask]), maximum(res[mask])],
    [1, 1];
    color     = :black,
    linestyle = :dash
)

# annotate correlation
text!(ax,
    "r = $(round(r_val, digits=2))";
    position = (maximum(res[mask]), maximum(ratio[mask])),
    align    = (:right, :top),
    fontsize = 14
)

display(fig)
