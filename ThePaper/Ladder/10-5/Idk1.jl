
using CairoMakie, LinearAlgebra, Statistics

# — compute gap and rel_dev as before —
A.gap = [ begin
    p    = row.p_final
    B    = vcat(row.R_eq, row.C_eq)
    D,Mstar = compute_jacobian(B, p)
    J    = D*Mstar
    ev   = sort(real(eigvals(J)))
    ev[end] - ev[end-1]
end for row in eachrow(A) ]

# overall correlation
r_S6 = cor(A.rt_pulse_full, A.rt_pulse_S6)
println("Global r(step 6) = ", round(r_S6, digits=3))

# per‐run relative deviation
rel_dev_S6 = A.rt_pulse_S6 ./ A.rt_pulse_full

# — OPTION A: filter out the worst 5% of outliers —
valid   = .!isnan.(rel_dev_S6) .& isfinite.(rel_dev_S6)
cutoff  = quantile(rel_dev_S6[valid], 0.9)   # 95th percentile
mask_A  = valid .& (rel_dev_S6 .<= cutoff)

# — OPTION B: keep all but plot on log‐scale (comment out if you prefer A) —
mask_B  = valid

# choose mask = mask_A  or  mask_B
mask = mask_A   # ← drop top 5% outliers
# mask = mask_B   # ← keep all, but we'll use log‐scale below

# — now make the scatter —
fig = Figure(resolution = (600,400))
ax  = Axis(fig[1,1];
    xlabel   = "Spectral gap (λ₂−λ₁)",
    ylabel   = "rt_S6 / rt_full",
    title    = "Step 6 deviation vs spectral gap",
    yscale   = mask===mask_B ? log10 : identity   # log‐scale if keeping all
)

# draw points
scatter!(ax, A.gap[mask], rel_dev_S6[mask];
         color      = :blue,
         markersize = 6,
         alpha      = 0.7)

# reference line at ratio = 1
lines!(ax, [minimum(A.gap), maximum(A.gap)], [1,1];
       color=:black, linestyle=:dash)

# annotate global r
text!(ax, "r = $(round(r_S6, digits=2))";
      position = (maximum(A.gap), maximum(rel_dev_S6[mask])),
      align    = (:right,:top),
      fontsize = 14)

display(fig)
