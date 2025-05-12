# 1) Full‐model Dii minima vs resilience_full
minD_full = Float64[]
res_full  = Float64[]

res_cols = Symbol.("resilience_" .* step_keys)
A = filter(row -> all(row[c] < 0 for c in res_cols), A)

for row in eachrow(A)
    # unpack p_final and equilibrium
    R, C, m_cons, xi_cons, r_res, d_res, ε, _ = row.p_final
    B0 = row.B0
    S  = R + C

    # build the diagonal D
    D = zeros(S)
    # resources
    for i in 1:R
        D[i] = d_res[i] * B0[i]
    end
    # consumers
    for k in 1:C
        i = R + k
        D[i] = (m_cons[k] / xi_cons[k]) * B0[i]
    end

    push!(minD_full, minimum(D))
    push!(res_full,  row.resilience_full)
end

println("Full‐model cor(min D_{ii}, R_∞) = ", cor(minD_full, res_full))

# scatterplot
fig = Figure(resolution=(500,400))
ax  = Axis(fig[1,1];
    xlabel="min D_{ii}", ylabel="resilience_full",
    title="Full model")
scatter!(ax, minD_full, res_full; markersize=6)
# lines!(ax, extremes(minD_full), extremes(minD_full); color=:black, linestyle=:dash)
display(fig)


# 2) Repeat for each ladder step
for s in 1:16
    minD_step = Float64[]
    res_step  = Float64[]

    for row in eachrow(A)
        # tau_Ss is the S‐vector of 1/B_eq at step s
        tau = row[Symbol("tau_S$s")]
        # skip if missing/NaN
        if any(isnan, tau)
            continue
        end

        # reconstruct eq = 1 ./ tau, and Dii = eq (since m_cons_hat/xi_hat = 1, d_res_hat = 1)
        D = 1.0 ./ tau

        push!(minD_step, minimum(D))
        push!(res_step,  row[Symbol("resilience_S$s")])
    end

    println("Step $s: cor(min D_{ii}, R_∞) = ", cor(minD_step, res_step))

    # optional: make a mini‐scatter for step s
    fig = Figure(resolution=(400,350))
    ax  = Axis(fig[1,1];
        xlabel="min D_{ii} (step $s)",
        ylabel="resilience_S$s",
        title="Step $s")
    scatter!(ax, minD_step, res_step; markersize=5, color=:blue)
    # lines!(ax, extremes(minD_step), extremes(minD_step); color=:black, linestyle=:dash)
    display(fig)
end
