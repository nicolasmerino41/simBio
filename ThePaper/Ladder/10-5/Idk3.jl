
using LinearAlgebra, Statistics

"""
    satellite_sensitivity(B_eq, J_full; t=1.0, n=500, rare_quantile=0.1)

Given full‐model equilibrium B_eq (length S) and Jacobian J_full,
identifies the “rare” species whose biomass ≤ the `rare_quantile`‐th
quantile of B_eq.  For each such species i, removes i from the community
(reducing J_full to J_minus of size (S−1)×(S−1)), then returns a
DataFrame with columns:

  • i                – the index of the removed species  
  • Rinf_removed     – new asymptotic resilience of the S−1 community  
  • ΔRinf            – Rinf_removed − Rinf_full  
  • Rmed_removed     – finite‐time median return‐rate for J_minus  
  • ΔRmed            – Rmed_removed − Rmed_full  
"""
function satellite_sensitivity(B_eq::AbstractVector, J_full::AbstractMatrix;
                               t::Real=1.0, n::Int=500, rare_quantile::Real=0.1)

    S = length(B_eq)
    # baseline metrics
    ev_full  = sort(real(eigvals(J_full)))
    Rinf_full = -ev_full[end]
    Rmed_full = median_return_rate(J_full, B_eq; t=t, n=n)

    # find rare species
    thr = quantile(B_eq, rare_quantile)
    rares = findall(x -> x <= thr, B_eq)
    rows = Vector{NamedTuple}()

    for i in rares
        # build reduced Jacobian
        keep = setdiff(1:S, i)
        J_minus = J_full[keep, keep]
        B_minus = B_eq[keep]

        # recompute resilience
        ev_m    = sort(real(eigvals(J_minus)))
        Rinf_m  = -ev_m[end]

        # recompute finite‐time rate
        Rmed_m  = median_return_rate(J_minus, B_minus; t=t, n=n)

        push!(rows, (
            i              = i,
            Rinf_removed   = Rinf_m,
            ΔRinf          = Rinf_m - Rinf_full,
            Rmed_removed   = Rmed_m,
            ΔRmed          = Rmed_m - Rmed_full
        ))
    end

    return DataFrame(rows)
end

# 1) test satellite‐species sensitivity in the full model
sat_full = satellite_sensitivity(B_eq, J_full; t=1.0, n=500, rare_quantile=0.1)
# you could store `sat_full` or summary stats, e.g. maximum |ΔRinf| or |ΔRmed|
max_dRinf  = maximum(abs.(sat_full.ΔRinf))
max_dRmed  = maximum(abs.(sat_full.ΔRmed))

sat_step = satellite_sensitivity(B_eq, J_s; t=1.0, n=500, rare_quantile=0.1)
max_dRinf_s = maximum(abs.(sat_step.ΔRinf))
max_dRmed_s = maximum(abs.(sat_step.ΔRmed))

rec = (
  …,
  Rinf_full      = -maximum(real(eigvals(J_full))),
  Rmed_full      = Rmed_full,
  max_dRinf_full = max_dRinf,
  max_dRmed_full = max_dRmed,
  # inside step_pairs flattening, for each step:
    Symbol("max_dRinf_S$step") => max_dRinf_s,
    Symbol("max_dRmed_S$step") => max_dRmed_s,
  …
)

for row in eachrow(A)
    τ_full = row.tau_full
    τ_step = row[Symbol("tau_S5")]  # vector
    rt     = row.rt_pulse_full_vector     # however you stored the per‐species pulse‐RTs

    fig = Figure()
    ax  = Axis(fig[1,1];)
    scatter!(ax, τ_full, rt; color=:blue)
    lines!(ax, [0,maximum(rt)], [0,maximum(rt)]; color=:black, linestyle=:dash)
    display(fig)

    # species‐by‐species correlation
    println("cor(τ_full,τ_step) = ", cor(τ_full, rt))
    # println("Stable is $(row.stable_full)")
end
for (i, row) in enumerate(eachrow(A))

  τ_full = row.tau_full
  rt     = row.rt_pulse_full_vector

  # clean up (in case of NaNs)
  mask = .!isnan.(τ_full) .& .!isnan.(rt)
  x = τ_full[mask]
  y = rt[mask]
  if length(x) < 2
      @warn "Not enough points to fit regression on instance $i"
      continue
  end

  # compute Pearson r and regression coefficients
  r = cor(x, y)
  m = cov(x,y) / var(x)
  b = mean(y) - m * mean(x)

  # bounds for drawing lines and placing text
  lo, hi = minimum(x), maximum(x)
  ylo, yhi = minimum(y), maximum(y)

  # build the figure
  fig = Figure(; size = (500, 500))
  ax = Axis(fig[1,1];
      xlabel = "Theoretical τᵢ = 1 / Bᵢ*",
      ylabel = "Simulated pulse-return time",
      title  = "Case $i: r=$(round(r, digits=2)))"
  )

  # scatter & regression line
  scatter!(ax, x, y; color = :blue, markersize = 6)
  lines!(ax, [lo, hi], [b + m*lo, b + m*hi];
        color = :red, linewidth = 2)

  # add the after_persistence_full text in top-left
  text!(
    ax,
    "collectivity = $(round(row.collectivity_full, digits=2)), \n" *
    "after_persistence = $(round(row.after_persistence_full, digits=2)), \n" *
    "resilience = $(round(row.resilience_full, digits=2)), \n" * 
    "reactivity = $(round(row.reactivity_full, digits=2)), \n";
    position = (lo, yhi),
    align    = (:left, :top),
    fontsize = 12,
    color    = :black
  )

  display(fig)

  # print out fit details
  println("  Pearson r(τ, RT) = ", round(r, digits=2))
  # println("  regression: RT = $(round(b,2)) + $(round(m,2))*τ")
end

function fit_tau_rt_cor_drivers(A::DataFrame)
  # 1) compute r_i for each community
  rs = Float64[]
  for row in eachrow(A)
      τ = row.tau_full
      rt = row.rt_pulse_full_vector
      mask = .!isnan.(τ) .& .!isnan.(rt)
      if sum(mask) < 2
          push!(rs, NaN)
      else
          push!(rs, cor(τ[mask], rt[mask]))
      end
  end

  # 2) build meta‐DataFrame
  meta = DataFrame(
      r                      = rs,
      conn                   = A.conn,
      IS                     = A.IS,
      scen                   = A.scen,
      delta                  = A.delta,
      eps                    = A.eps,
      after_persistence_full = A.after_persistence_full,
      collectivity_full      = A.collectivity_full,
      resilience_full        = A.resilience_full,
      reactivity_full        = A.reactivity_full,
      Rmed_full              = A.Rmed_full,
  )

  # drop any rows where r is NaN
  dropmissing!(meta, :r)

  # make scenario a categorical predictor
  # meta.scen = CategoricalArray(meta.scen)

  # 3) fit linear model
  f = @formula(r ~ conn + IS + scen + delta + eps +
                 after_persistence_full + collectivity_full +
                 resilience_full + reactivity_full + Rmed_full)
  lm = GLM.lm(f, meta)

  println("\n=== Model summary ===")
  display(coeftable(lm))

  return lm, meta
end

# usage:
lm, meta = fit_tau_rt_cor_drivers(A)