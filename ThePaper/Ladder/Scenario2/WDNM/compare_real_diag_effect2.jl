# ─────────────────────────────────────────────────────────────────────────────
# 1) Core Δ‐metric functions
# ─────────────────────────────────────────────────────────────────────────────
measure_resilience(J) = -maximum(real(eigvals(J)))
measure_reactivity(J) = maximum(real(eigvals((J+J')/2)))
function analytical_Rmed(J; t=0.01)
    E = exp(t*J); M = E*E'
    rates = -((diag(M).-1)./(2t))
    sorted = sort(rates)
    lo = ceil(Int, length(sorted)*0.05)+1
    hi = floor(Int, length(sorted)*0.95)
    mean(sorted[lo:hi])
end
pulse_return_time(J) = 1/measure_resilience(J)
measure_persistence(J) = count(<(0), real.(eigvals(J))) / size(J,1)
nonnormality(J) = norm(J*J' - J'*J)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Build one Jacobian via your full pipeline; ensure local stability
# ─────────────────────────────────────────────────────────────────────────────
function build_stable_J(R, C; conn, IS, scen, epsi, m_val, g_val, cb, MAX_TRIES=20)
    S = R + C
    for _ in 1:MAX_TRIES
        # 1) assemble A, ε
        A = make_A(zeros(S,S), R, conn, scen; IS=IS)
        ε = clamp.(fill(epsi,S,S), 0, 1)
        thr = generate_feasible_thresholds(A, ε, R)
        isempty(thr) && continue
        t = thr[1]
        xi_cons, K_res = t.xi_cons, t.K_res
        B0 = calibrate_from_K_xi(xi_cons, K_res, ε, A)
        m_cons = abs.(rand(Normal(m_val,0.2), C))
        r_res  = abs.(rand(Normal(g_val,0.2), R))
        p = (R, C, m_cons, xi_cons, r_res, K_res, ε, A)
        ok, _ = survives!(B0, p; cb=cb)
        ok || continue
        D, M = compute_jacobian(B0, p)
        J = D*M
        is_locally_stable(J) || continue
        return D, M, J, A
    end
    return nothing, nothing, nothing, nothing
    # error("Failed to build a stable J in $MAX_TRIES tries")
end

function trimmed(v::AbstractVector{<:Real}; p::Float64=0.05)
    n = length(v)
    if n == 0
        return NaN
    end
    s = sort(v)
    lo = floor(Int,  p*n) + 1
    hi = ceil(Int,   (1-p)*n)
    return mean(s[lo:hi])
end
# ─────────────────────────────────────────────────────────────────────────────
# 3) Main routine: sample grid, build refs, bin by IS and σ/min(d), reshuffle
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# 4) The full comparison pipeline
# ─────────────────────────────────────────────────────────────────────────────
function compare_real_diag_effect(
    ; S=50, C=20,
      conn_vals=[0.05,0.1,0.2],
      IS_vals=[0.01,0.1,1.0,2.0],
      scenarios=[:ER],
      eps_scales=[1.0],
      mortality_vals=[0.1,0.2,0.3,0.4,0.5],
      growth_vals=[0.5,1.0,3.0,5.0,7.0],
      number_of_combinations=500,
      iterations=1,
      reps_shuffle=100,
      nbins=8
)
    R = S - C
    # 1) sample parameter combos
    combos = collect(Iterators.product(conn_vals, IS_vals, scenarios,
                                       eps_scales, mortality_vals,
                                       growth_vals, 1:iterations))
    Random.seed!(123)
    combos = sample(combos, min(length(combos), number_of_combinations); replace=false)

    cb = build_callbacks(S, EXTINCTION_THRESHOLD)
    locki = ReentrantLock()

    raw = Vector{NamedTuple{(:combo_id, :IS, :sigma_md, :type, :metric, :Δ, :degree_cv),
                           Tuple{Int,Float64,Float64,Symbol,Symbol,Float64,Float64}}}()

    mets = [:Res, :Rea, :Rmed, :Pulse, :Pers, :Nonnorm]

    # raw = NamedTuple[]  # or use your existing Vector{NamedTuple}
    combo_id = 0
    # 2) threaded sweep
    @threads for tup in combos
        combo_id += 1
        (conn, IS, scen, epsi, m_val, g_val, _) = tup
        D0, M0, J0, A = build_stable_J(R,C;
            conn=conn, IS=IS, scen=scen,
            epsi=epsi, m_val=m_val, g_val=g_val, cb=cb)

        D0 === nothing && continue

        # compute sigma/min d
        d0       = abs.(-diag(J0))
        sigma_md = IS / minimum(d0)
        
        g = SimpleGraph(A .!= 0)         # unweighted graph
        degs = degree(g)
        degree_cv = std(degs) / mean(degs)
        
        # base metrics
        base = (
          Res    = measure_resilience(J0),
          Rea    = measure_reactivity(J0),
          Rmed   = analytical_Rmed(J0),
          Pulse  = pulse_return_time(J0),
          Pers   = measure_persistence(J0),
          Nonnorm= nonnormality(J0),
        )

        # stats of M0
        flat    = vec(M0)
        μ_all, σ_all = mean(flat), std(flat)
        off_idx = findall(!=(I(size(M0,1))), trues(size(M0)))
        μ_off, σ_off = mean(M0[off_idx]), std(M0[off_idx])
        d_idx       = LinearIndices(D0)[diagind(D0)]
        μ_diag, σ_diag = mean(D0[d_idx]), std(D0[d_idx])

        flatJ   = vec(J0)
        off_idx = findall(!=(I(size(J0,1)), trues(size(J0))))
        diag_idx= LinearIndices(J0)[diagind(J0)]

        # reshuffles
        for _ in 1:reps_shuffle
            # Mf = reshape(randn(length(flat)).*(σ_all).+μ_all, size(M0))
            # Df = copy(D0); Df[d_idx]   .= randn(length(d_idx)).*(σ_diag).+μ_diag
            # Jo = copy(M0); Jo[off_idx] .= randn(length(off_idx)).*(σ_off).+μ_off
            # Jd = copy(D0); Jd[d_idx]   .= randn(length(d_idx)).*(σ_diag).+μ_diag
            # Jf, Jo, Jd = Df*Mf, D0*Jo, Jd*M0

            # 1) full reshuffle of M
            flatM = vec(M0)
            flatD = vec(D0)
            # Mf = reshape(shuffle!(copy(flatM)), size(M0))
            # Df = reshape(shuffle!(copy(flatD)), size(D0))     
            Mf = reshape(copy(flatM).*(rand(length(flatM)).*(σ_all).+μ_all), size(M0))
            Df = reshape(copy(flatD).*(rand(length(flatD)).*(σ_diag).+μ_diag), size(D0))
            Jf    = Df * Mf
            
            # 2) off‐diag only: permute off‐diagonals of M
            M_off = copy(M0)
            off_idx = findall(!=(I(size(M0,1))), trues(size(M0)))
            # M_off[off_idx] = shuffle!(copy(M0[off_idx]))
            M_off[off_idx] = copy(M0[off_idx]).*(rand(length(off_idx)).*(σ_off).+μ_off)
            Jo = D0 * M_off

            # 3) diag‐only: permute diagonal of D, keep M fixed
            D_diag = copy(D0)
            d_idx  = LinearIndices(D0)[diagind(D0)]
            # D_diag[d_idx] = shuffle!(copy(D0[d_idx]))
            D_diag[d_idx] = copy(D0[d_idx]).*(rand(length(d_idx)).*(σ_diag).+μ_diag)
            Jd = D_diag * M0

            # full permute
            Jf = reshape(shuffle!(copy(flatJ)), size(J0))

            # off‐diag only
            Jo = copy(J0)
            Jo[off_idx] .= shuffle!(copy(J0[off_idx]))

            # diag‐only
            Jd = copy(J0)
            Jd[diag_idx] .= shuffle!(copy(J0[diag_idx]))

            vals = Dict(
              :full     => (
                Res=measure_resilience(Jf),
                Rea=measure_reactivity(Jf),
                Rmed=analytical_Rmed(Jf),
                Pulse=pulse_return_time(Jf),
                Pers=measure_persistence(Jf),
                Nonnorm=nonnormality(Jf),
              ),
              :offdiag  => (
                Res=measure_resilience(Jo),
                Rea=measure_reactivity(Jo),
                Rmed=analytical_Rmed(Jo),
                Pulse=pulse_return_time(Jo),
                Pers=measure_persistence(Jo),
                Nonnorm=nonnormality(Jo),
              ),
              :diagonly => (
                Res=measure_resilience(Jd),
                Rea=measure_reactivity(Jd),
                Rmed=analytical_Rmed(Jd),
                Pulse=pulse_return_time(Jd),
                Pers=measure_persistence(Jd),
                Nonnorm=nonnormality(Jd),
              ),
            )

            lock(locki) do
                for typ in (:full,:offdiag,:diagonly), met in mets
                    Δ = abs(getfield(vals[typ],met) - getfield(base,met))
                    push!(raw, (combo_id, IS, sigma_md, typ, met, Δ, degree_cv))
                end
            end
        end
    end

    # 3) build DataFrame & compute bin‐indices manually
    df = DataFrame(raw)
    # define the bin edges on the **log10** scale
    is_edges = range(minimum(df.IS),  maximum(df.IS),  length=nbins+1)
    sm_edges = 10 .^ range(log10(minimum(df.sigma_md)), log10(maximum(df.sigma_md)), length=nbins+1)
    # for each row, find the bin (1…nbins)
    df.IS_bin  = [ clamp(searchsortedfirst(is_edges,  x) - 1, 1, nbins) for x in df.IS ]
    df.SM_bin  = [ clamp(searchsortedfirst(sm_edges, x) - 1, 1, nbins) for x in df.sigma_md ]

    # 4) aggregate trimmed‐means for each metric & reshuffle type
    aggIS = combine(groupby(df, [:IS_bin, :type, :metric]),
                    :Δ => (v-> trimmed(v)) => :Δ_trim)
    aggSM = combine(groupby(df, [:SM_bin, :type, :metric]),
                    :Δ => (v-> trimmed(v)) => :Δ_trim)

    # 5) pivot for ease of plotting
    pivotIS = unstack(aggIS,  [:IS_bin, :metric], :type, :Δ_trim)
    pivotSM = unstack(aggSM,  [:SM_bin, :metric], :type, :Δ_trim)

    # 6) plot
    fig = Figure(; size=(1100,450))
    titles = ["Δ Res","Δ Rea","Δ Rₘₑd","Δ PulseRT","Δ Pers","Δ Nonnorm"]

    for (k,(pivot, edges, label)) in enumerate(((pivotIS, is_edges, "IS"), (pivotSM, sm_edges, "σ/min d")))
        for (i, met) in enumerate(mets)
            row, col = divrem(i-1,3)
            ax = Axis(fig[row+1, col*2 + k],
                    xlabelsize=8, ylabelsize=8,
                    xticklabelsize=8, yticklabelsize=8)
            mask = pivot.metric .== met    # use the column `metric`
            full    = pivot[mask, :full]
            offdiag = pivot[mask, :offdiag]
            diagonly= pivot[mask, :diagonly]

            xs = edges[1:end-1] .+ diff(edges)/2  # midpoints
            xx = xs[1:length(full)]
            lines!(ax, xx, full,    label="full")
            lines!(ax, xx, offdiag, label="off‐diag")
            lines!(ax, xx, diagonly,label="diag‐only")
            ax.xlabel = label
            ax.ylabel = titles[i]
            k==2 && (ax.xscale = log10)
            # row==0 && col==0 && axislegend(ax)
        end
    end
    display(fig)

    # ─────────────────────────────────────────────────────────────────────────
    # 7) RETURN the raw and pivoted tables so you can inspect counts, etc.
    return (
      raw     = df,        # every single Δ with its IS_bin & SM_bin
      summary = (is=pivotIS, σmin=pivotSM),
      bins    = (IS_edges=is_edges, SM_edges=sm_edges)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 7) Run it!
# ─────────────────────────────────────────────────────────────────────────────
out1 = compare_real_diag_effect(
  S=50, C=20;
  conn_vals=0.01:0.02:1.0,
  IS_vals=0.01:0.01:2.0,
  scenarios=[:ER],
  eps_scales=[1.0],
  mortality_vals=[0.1,0.2,0.3,0.4,0.5],
  growth_vals=[0.5,1.0,3.0,5.0,7.0],
  number_of_combinations=1000,
  iterations=5,
  reps_shuffle=50,
  nbins=8
)

# how many Δ's fell in each IS‐bin?
count_by_ISbin = combine(groupby(out.raw, :IS_bin), nrow => :count)
# counts in each combination of IS‐bin × reshuffle‐type
xt = combine(groupby(out.raw, [:IS_bin, :type]), :Δ => length => :n_obs)

# How many Δ’s fell in each IS‐bin?
count_by_ISbin = combine(groupby(outER.raw, :IS_bin),
                        nrow => :n_obs)
count_by_ISbin = combine(groupby(outPL.raw, :IS_bin),
                        nrow => :n_obs)
count_by_ISbin = combine(groupby(outMOD.raw, :IS_bin),
                        nrow => :n_obs)
count_by_ISbin = combine(groupby(outALL.raw, :IS_bin),
                        nrow => :n_obs)

# How many _distinct_ communities in each IS‐bin?
uniq_by_ISbin  = combine(groupby(outER.raw, :IS_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_ISbin  = combine(groupby(outPL.raw, :IS_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_ISbin  = combine(groupby(outMOD.raw, :IS_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_ISbin  = combine(groupby(outALL.raw, :IS_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)

# Similarly for σ/min d
uniq_by_SMbin  = combine(groupby(outER.raw, :SM_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_SMbin  = combine(groupby(outPL.raw, :SM_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_SMbin  = combine(groupby(outMOD.raw, :SM_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
uniq_by_SMbin  = combine(groupby(outALL.raw, :SM_bin),
                         :combo_id => (x-> length(unique(x))) => :n_communities)
