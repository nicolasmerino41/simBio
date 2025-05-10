
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
    lines!(ax, [0,maximum(τ_full)], [0,maximum(τ_full)]; color=:black, linestyle=:dash)
    display(fig)

    # species‐by‐species correlation
    # println("Step $k  cor(τ_full,τ_step) = ", cor(τ_full, τ_step))
    println("Stable is $(row.stable_full)")
end