
function compare_real_diag_effect(
    S::Int=50, C::Int=20;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[1.0],
    scenarios=[:ER],
    eps_scales=[1.0],
    mortality_vals=[1.0, 0.1, 0.2, 0.3, 0.4, 0.5],
    growth_vals=[0.03, 0.5, 1.0, 3.0, 5.0, 7.0],
    number_of_combinations = 100,
    iterations=1
)
    R = S - C
    results = Vector{NamedTuple}()
    combos = collect(Iterators.product(conn_vals,IS_vals,scenarios,eps_scales,1:iterations,mortality_vals,growth_vals))

    @info "Number of combinations: $(length(combos))"
    for (conn,IS,scen,epsi,ite,m_val,g_val) in Iterators.take(shuffle!(combos), number_of_combinations)

        # 1) build A, epsilon, get thresholds, equilibria, jacobian
        A = make_A(zeros(S,S),R,conn,scen; IS=IS,
                   pareto_exponent=1.25, pareto_minimum_degree=1.0, mod_gamma=1.0)
        epsilon = clamp.(fill(epsi,S,S),0,1)
        thr_sets = generate_feasible_thresholds(A, epsilon, R)
        isempty(thr_sets) && continue
        t = thr_sets[1]
        xi_cons, K_res = t.xi_cons, t.K_res

        B0 = calibrate_from_K_xi(xi_cons, K_res, epsilon, A)
        m_cons = abs.(rand(Normal(m_val, 0.2), C))
        r_res  = abs.(rand(Normal(g_val, 0.2), R))
        p = (R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)

        ok, _ = survives!(B0, p; cb=build_callbacks(S,1e-6))
        !ok && continue

        D0, M0 = compute_jacobian(B0, p)
        J0 = D0 * M0

        # containers
        metrics = (:Res,:Rea,:Rmed,:Pulse,:Pers,:Nonnorm)
        Δ_full = Dict(m=>Float64[] for m in metrics)
        Δ_diag = Dict(m=>Float64[] for m in metrics)
        Δ_diag_fixed = Dict(m=>Float64[] for m in metrics)

        # 2) J1: off-diag only resample of M0
        Sdim = size(M0,1)
        idx_off = findall(!=(I(Sdim)), trues(Sdim,Sdim))
        vals_off = M0[idx_off]; μ_off, σ_off = mean(vals_off), std(vals_off)
        M1 = copy(M0)
        M1[idx_off] .= randn(length(idx_off)) .* (2σ_off) .+ 2μ_off
        J1 = D0 * M1

        # 3) J2: full-matrix reshuffle of M0
        flat = vec(M0); μ_all, σ_all = mean(flat), std(flat)
        M2 = copy(M0)
        M2[:] .= randn(length(flat)) .* (2σ_all) .+ 2μ_all
        J2 = D0 * M2

        # 4) J3: diag-only resample of M0
        d_idx = LinearIndices(M0)[diagind(M0)]
        flat_diag = M0[d_idx]; μ_diag, σ_diag = mean(flat_diag), std(flat_diag)
        M3 = copy(M0)
        M3[d_idx] .= randn(length(d_idx)) .* (2σ_diag) .+ 2μ_diag
        J3 = D0 * M3

        # 5) compute deltas against J0
        for (Jfull, Jdiag, dict_full, dict_diag, label) in
            ((J2, J1, Δ_full, Δ_diag, "off-diag"),
             (J3, J1, Δ_full, Δ_diag_fixed, "diag-only"))
            
            push!(dict_full[:Res],  abs(-maximum(real(eigvals(Jfull)))   - (-maximum(real(eigvals(J0))))))
            push!(dict_full[:Rea],  abs(maximum(real(eigvals((Jfull+Jfull')/2))) -
                                          maximum(real(eigvals((J0+J0')/2)))))
            push!(dict_full[:Rmed], abs( analytical_Rmed(Jfull) - analytical_Rmed(J0) ))
            push!(dict_full[:Pulse],abs( 1/(-maximum(real(eigvals(Jfull)))) -
                                          1/(-maximum(real(eigvals(J0)))) ))
            push!(dict_full[:Pers], abs( count(<(0), real.(eigvals(Jfull))) / length(eigvals(Jfull)) -
                                         count(<(0), real.(eigvals(J0)))    / length(eigvals(J0)) ))
            d0 = -diag(J0)
            # push!(dict_full[:Coll], 0.0 )  # collectivity unchanged by these tests
            push!(dict_full[:Nonnorm], abs( norm(Jfull*Jfull' - Jfull'*Jfull) -
                                            norm(J0*J0'       - J0'*J0) ))

            push!(dict_diag[:Res],  abs(-maximum(real(eigvals(Jdiag)))   - (-maximum(real(eigvals(J0))))))
            push!(dict_diag[:Rea],  abs(maximum(real(eigvals((Jdiag+Jdiag')/2))) -
                                          maximum(real(eigvals((J0+J0')/2)))))
            push!(dict_diag[:Rmed], abs( analytical_Rmed(Jdiag) - analytical_Rmed(J0) ))
            push!(dict_diag[:Pulse],abs( 1/(-maximum(real(eigvals(Jdiag)))) -
                                          1/(-maximum(real(eigvals(J0)))) ))
            push!(dict_diag[:Pers], abs( count(<(0), real.(eigvals(Jdiag))) / length(eigvals(Jdiag)) -
                                         count(<(0), real.(eigvals(J0)))    / length(eigvals(J0)) ))
            # push!(dict_diag_fixed[:Coll], 0.0)
            push!(dict_diag[:Nonnorm], abs( norm(Jdiag*Jdiag' - Jdiag'*Jdiag) -
                                            norm(J0*J0'       - J0'*J0) ))
        end

        # 6) save row
        push!(results, (
            conn, IS, scen, epsi,
            Δ_res_full     = mean(Δ_full[:Res]),     Δ_res_offdiag     = mean(Δ_diag[:Res]),
            Δ_res_diagonly = mean(Δ_diag_fixed[:Res]),
            Δ_rea_full     = mean(Δ_full[:Rea]),     Δ_rea_offdiag     = mean(Δ_diag[:Rea]),
            Δ_rea_diagonly = mean(Δ_diag_fixed[:Rea]),
            Δ_rmed_full    = mean(Δ_full[:Rmed]),    Δ_rmed_offdiag    = mean(Δ_diag[:Rmed]),
            Δ_rmed_diagonly= mean(Δ_diag_fixed[:Rmed]),
            Δ_pulse_full   = mean(Δ_full[:Pulse]),   Δ_pulse_offdiag   = mean(Δ_diag[:Pulse]),
            Δ_pulse_diagonly=mean(Δ_diag_fixed[:Pulse]),
            Δ_pers_full    = mean(Δ_full[:Pers]),    Δ_pers_offdiag    = mean(Δ_diag[:Pers]),
            Δ_pers_diagonly=mean(Δ_diag_fixed[:Pers]),
            # Δ_coll_full    = mean(Δ_full[:Coll]),    Δ_coll_offdiag    = mean(Δ_diag[:Coll]),
            # Δ_coll_diagonly=mean(Δ_diag_fixed[:Coll]),
            Δ_nonnorm_full = mean(Δ_full[:Nonnorm]), Δ_nonnorm_offdiag = mean(Δ_diag[:Nonnorm]),
            Δ_nonnorm_diagonly=mean(Δ_diag_fixed[:Nonnorm])
        ))
    end

    df = DataFrame(results)
    if isempty(df)
        @warn "No valid runs—returning empty DataFrame"
        return df
    end

    # 7) plot
    metrics = [:Res,:Rea,:Rmed,:Pulse,:Pers,:Nonnorm]
    titles = ["Δ Res","Δ Rea","Δ Rmed","Δ Pulse","Δ Pers","Δ Nonnorm"]
    fig = Figure(; size=(1200,300))
    Label(fig[0, 1:3], "σ = $(IS_vals[1])"; fontsize=14)
    for (i, m) in enumerate(metrics)
        ax = Axis(
            fig[1,i];
            xlabelsize=6,
            xticklabelrotation=π/4,
            xticksize=6, ylabelsize=10
            )
        xs = 1:3
        ys = [ mean(df[!, Symbol("Δ_$(lowercase(string(m)))_full")]),
               mean(df[!, Symbol("Δ_$(lowercase(string(m)))_offdiag")]),
               mean(df[!, Symbol("Δ_$(lowercase(string(m)))_diagonly")]) ]
        barplot!(ax, xs, ys)
        ax.xticks = (xs, ["full","off-diag","diag-only"])
        ax.title = titles[i]
    end
    display(fig)

    return df
end

# run across σ’s
for sigma in [0.01, 0.1, 1.0, 2.0, 5.0]
    df = compare_real_diag_effect(
      50, 20;
      conn_vals=0.01:0.02:1.0,
      IS_vals=[sigma],
      scenarios=[:ER],
      eps_scales=[1.0],
      number_of_combinations=1000,
      iterations=10
    )
    println("Done σ=$sigma, got $(nrow(df)) rows")
end
