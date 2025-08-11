using DataFrames
using LinearAlgebra
using Statistics
using CairoMakie

# =============== Structural metrics (your definitions) =================
function structural_from_A(A::AbstractMatrix)
    S = size(A, 1)
    Adj = A .!= 0.0
    @inbounds for i in 1:S
        Adj[i,i] = false
    end

    conn = sum(Adj) / (S * (S - 1))

    degs = sum(Adj, dims=2)
    mdeg = mean(degs)
    degree_cv = mdeg == 0 ? NaN : std(degs) / mdeg

    k = degs[:]
    m = sum(k) / 2
    modu = NaN
    if m > 0
        B = zeros(Float64, S, S)
        @inbounds for u in 1:S, v in 1:S
            B[u,v] = (Adj[u,v] ? 1.0 : 0.0) - (k[u]*k[v])/(2m)
        end
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        svec = map(x -> x >= 0 ? 1.0 : -1.0, v1)
        modu = (svec' * (B * svec)) / (4m)
    end

    return (modularity=modu, degree_cv=degree_cv, connectance=conn)
end

function add_structural_columns!(R::DataFrame)
    n = nrow(R)
    modularity  = Vector{Float64}(undef, n)
    degree_cv   = Vector{Float64}(undef, n)
    connectance = Vector{Float64}(undef, n)

    @inbounds for (i, row) in enumerate(eachrow(R))
        A = row.p_final[2]       # tuple; A at index 2
        mets = structural_from_A(A)
        modularity[i]  = mets.modularity
        degree_cv[i]   = mets.degree_cv
        connectance[i] = mets.connectance
    end

    R.modularity  = modularity
    R.degree_cv   = degree_cv
    R.connectance = connectance
    return R
end

# ===================== Binning utilities ======================
"""
    bin_xy(xs, ys; n_bins=50, method=:equalwidth, qclip=(nothing, nothing),
           agg_error=:std, min_bin=5)

Return (mx, my, ey, nx) where:
- mx: bin centers (mean x in bin)
- my: mean y per bin
- ey: error per bin (std or sem)
- nx: counts per bin
`method` = :equalwidth or :quantile
`qclip` = (qx_lo, qx_hi) to drop x outliers before binning (e.g., (0.01, 0.99))
`agg_error` = :std or :sem
`min_bin` = minimum points required to keep a bin
"""
function bin_xy(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
                n_bins::Int=50,
                method::Symbol=:equalwidth,
                qclip::Tuple{Union{Nothing,Real},Union{Nothing,Real}}=(nothing, nothing),
                agg_error::Symbol=:std,
                min_bin::Int=5)

    @assert length(xs) == length(ys)
    x = collect(skipmissing(xs))
    y = collect(skipmissing(ys))

    # finite mask (fresh, not reused after slicing)
    keep0 = map(i -> isfinite(x[i]) && isfinite(y[i]), eachindex(x))
    x = x[keep0];  y = y[keep0]

    # optional quantile clipping — use fresh masks each time
    qlo, qhi = qclip
    if qlo !== nothing
        xlo = quantile(x, Float64(qlo))
        keep1 = x .>= xlo
        x = x[keep1];  y = y[keep1]
    end
    if qhi !== nothing
        xhi = quantile(x, Float64(qhi))
        keep2 = x .<= xhi
        x = x[keep2];  y = y[keep2]
    end

    if isempty(x)
        return Float64[], Float64[], Float64[], Int[]
    end

    # build edges
    edges = if method == :equalwidth
        range(minimum(x), maximum(x), length=n_bins+1)
    elseif method == :quantile
        qs = collect(range(0.0, 1.0; length=n_bins+1))
        ed = quantile(x, qs)
        # remove non-increasing duplicates; ensure at least 2 edges
        uniq = [ed[1]]
        for i in 2:length(ed)
            if ed[i] > last(uniq)
                push!(uniq, ed[i])
            end
        end
        if length(uniq) < 2
            # fallback: equal-width on the degenerate range
            uniq = collect(range(minimum(x), maximum(x), length=min(n_bins+1, max(2, length(x)))))
        end
        uniq
    else
        error("Unknown binning method: $method")
    end

    nb = max(1, length(edges) - 1)
    # digitize
    bix = searchsortedlast.(Ref(edges), x)        # in [1, nb+1]
    @inbounds for i in eachindex(bix)
        bix[i] = clamp(bix[i], 1, nb)
    end

    mx = Float64[]; my = Float64[]; ey = Float64[]; nx = Int[]
    for b in 1:nb
        idxs = findall(==(b), bix)
        if length(idxs) >= min_bin
            xb = @view x[idxs]; yb = @view y[idxs]
            push!(mx, mean(xb))
            push!(my, mean(yb))
            if agg_error == :std
                push!(ey, std(yb))
            elseif agg_error == :sem
                push!(ey, std(yb) / sqrt(length(yb)))
            else
                error("agg_error must be :std or :sem")
            end
            push!(nx, length(idxs))
        end
    end

    # sort by x so lines/bands draw cleanly
    if !isempty(mx)
        ord = sortperm(mx)
        mx = mx[ord]; my = my[ord]; ey = ey[ord]; nx = nx[ord]
    end
    return mx, my, ey, nx
end

# ===================== Makie plotting ======================
"""
    plot_full_binned_lines(R; n_bins=50, error_kind=:sem, show_band=true,
                           bin_method=:equalwidth, qclip=(nothing,nothing),
                           save_plot=false)

Create a grid: rows = dynamics (resilience_full, reactivity_full, rt_pulse_full, after_press_full),
cols = structure (modularity, degree_cv, connectance). For each panel, bin x and draw mean(y) with errors.
"""
function plot_full_binned_lines(R::DataFrame;
        n_bins::Int=50,
        error_kind::Symbol=:sem,       # :std or :sem
        show_band::Bool=true,          # shaded band vs discrete errorbars
        bin_method::Symbol=:equalwidth,# or :quantile
        qclip::Tuple{Union{Nothing,Real},Union{Nothing,Real}}=(nothing, nothing),
        save_plot::Bool=false)

    dyns    = [:resilience_full, :reactivity_full, :rt_pulse_full, :after_press_full]
    dynname = Dict(
        :resilience_full  => "Resilience",
        :reactivity_full  => "Reactivity",
        :rt_pulse_full    => "Return Time (pulse)",
        :after_press_full => "Persistence (after press)",
    )
    structs = [:modularity, :degree_cv, :connectance]
    sname   = Dict(:modularity=>"Modularity", :degree_cv=>"Degree CV", :connectance=>"Connectance")

    # verify columns exist
    for c in vcat(dyns, structs)
        @assert hasproperty(R, c) "Missing column $(c). Did you run add_structural_columns! ?"
    end

    fig = Figure(; size=(1200, 1400))

    for (ri, y) in enumerate(dyns)
        for (ci, x) in enumerate(structs)
            ax = Axis(fig[ri, ci];
                title = "$(dynname[y]) vs $(sname[x])",
                xlabel = sname[x],
                ylabel = dynname[y],
                xgridvisible=false, ygridvisible=false
            )

            xs = R[!, x]
            ys = R[!, y]
            mx, my, ey, _ = bin_xy(xs, ys;
                                   n_bins=n_bins,
                                   method=bin_method,
                                   qclip=qclip,
                                   agg_error=error_kind,
                                   min_bin=5)

            # light scatter background (optional; comment out if unwanted)
            # scatter!(ax, xs, ys; color=:gray, markersize=2, transparency=true)

            # main line + error
            lines!(ax, mx, my; linewidth=2)
            if show_band
                # shaded error band
                upper = my .+ ey
                lower = my .- ey
                band!(ax, mx, lower, upper; transparency=true)
            else
                errorbars!(ax, mx, my, ey)
            end
        end
    end

    display(fig)
    if save_plot
        save("full_binned_lines.png", fig; px_per_unit=6.0)
    end
    return fig
end

# ===================== Typical pipeline ======================
add_structural_columns!(R)
fig = plot_full_binned_lines(R; n_bins=50, error_kind=:sem, show_band=true,
                             bin_method=:quantile, qclip=(0.01, 0.99), save_plot=false)
