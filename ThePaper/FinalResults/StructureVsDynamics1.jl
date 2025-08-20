using DataFrames
using LinearAlgebra
using Statistics
using CairoMakie

# ───────────────────────────────────────────────────────────────
# Structural metrics from your snippet
# ───────────────────────────────────────────────────────────────
function structural_from_A(A::AbstractMatrix)
    S = size(A, 1)
    Adj = A .!= 0.0
    # kill self-loops if present
    for i in 1:S
        Adj[i, i] = false
    end

    conn = sum(Adj) / (S * (S - 1))

    degs = sum(Adj, dims=2)
    mdeg = mean(degs)
    degree_cv = mdeg == 0 ? NaN : std(degs) / mdeg

    k = degs[:]
    m = sum(k) / 2
    B = zeros(Float64, S, S)
    for u in 1:S, v in 1:S
        B[u,v] = (Adj[u,v] ? 1.0 : 0.0) - (k[u] * k[v]) / (2m)
    end
    vals, vecs = eigen(Symmetric(B))
    v1 = vecs[:, argmax(vals)]
    svec = map(x -> x >= 0 ? 1.0 : -1.0, v1)
    modu = (svec' * (B * svec)) / (4m)

    return (modularity=modu, degree_cv=degree_cv, connectance=conn)
end

function add_structural_columns!(R::DataFrame)
    n = nrow(R)
    modularity  = Vector{Float64}(undef, n)
    degree_cv   = Vector{Float64}(undef, n)
    connectance = Vector{Float64}(undef, n)

    for (i, row) in enumerate(eachrow(R))
        A = row.p_final[2]  # tuple, A is index 2
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

# ───────────────────────────────────────────────────────────────
# Makie plotting
# ───────────────────────────────────────────────────────────────
function plot_full_scatter(R::DataFrame; save_plot=false, pixels_per_unit=3)
    dyns = [:resilience_full]
    structs = [:modularity, :degree_cv, :connectance]

    # Drop rows with missing values in any relevant column
    cols = vcat(dyns, structs)
    Rclean = dropmissing(copy(R), cols)

    fig = Figure(
        ; size=(600, 450),
        )

    for (i, y) in enumerate(dyns)
        for (j, x) in enumerate(structs)
            ax = Axis(
                fig[i, j],
                xlabel=String(x), ylabel=String(y),
                yscale=log10,
                # xscale=log10
                )
            scatter!(ax, Rclean[!, x], abs.(Rclean[!, y]), color=:blue, markersize=5, transparency=true)
        end
    end

    display(fig)
    if save_plot
        filename = "structural_vs_dynamics.png"
        save(filename, fig; px_per_unit=pixels_per_unit)
    end
    # return fig
end

# ───────────────────────────────────────────────────────────────
# Pipeline:
# ───────────────────────────────────────────────────────────────
add_structural_columns!(R_subset)
fig = plot_full_scatter(R_subset; save_plot=true)

function plot_one_scatter(R::DataFrame, xcol::Symbol, ycol::Symbol)
    Rclean = dropmissing(copy(R), [xcol, ycol])
    fig = Figure(; size=(500, 400))
    ax = Axis(fig[1, 1], xlabel=String(xcol), ylabel=String(ycol))
    scatter!(ax, Rclean[!, xcol], Rclean[!, ycol], color=:blue, markersize=6, transparency=true)
    display(fig)
    return fig
end

# Example:  
plot_one_scatter(R, :connectance, :reactivity_full)