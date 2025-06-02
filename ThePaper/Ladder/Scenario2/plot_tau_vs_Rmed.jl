using CairoMakie, DataFrames

"""
    plot_species_tau_vs_Rmed(df;
        which    = :all,         # :all, :resources or :consumers
        conn_col = :conn        # column to color‐code by
    )

Assumes each row of `df` has:
  • `:tau`       ⇒ Vector{Float64} of length 50 (species τᵢ),
  • `:Rmed_full` ⇒ scalar median return rate,
  • `conn_col`   ⇒ value to color points by.

Plots 1/τᵢ vs Rₘₑd for the selected species.
"""
function plot_tau_vs_Rmed(
    df::DataFrame;
    which::Symbol    = :all,
    conn_col::Symbol = :conn
)
    # fixed community sizes
    n_res = 30
    n_cons = 20
    N = n_res + n_cons

    # Fix: avoid line breaks in ternary and range expressions, and add spaces around ':'
    idxs = if which === :resources
        1:n_res
    elseif which === :consumers
        (n_res+1):N
    else
        1:N
    end

    xs = Float64[]
    ys = Float64[]
    cs = Float64[]

    # gather data
    for row in eachrow(df)[4:4]
        τ = row[:tau_full]                # Vector of length 50
        # Rmed = row[:ssp_rmed_full]          # scalar median return rate
        Rmed = row[:analytical_rmed_full]          # scalar median return rate
        clr = row[conn_col]          # e.g. connectivity

        for i in idxs
            push!(xs, 1/τ[i])        # 1/τᵢ
            push!(ys, Rmed[i])          # same Rmed for all species
            push!(cs, clr)
        end
    end

    # make the scatter
    fig = Figure(resolution = (800, 400))
    ax = Axis(fig[1,1],
        xlabel = "1/τᵢ (species-level return rate)",
        ylabel = "Analytical Rₘₑd (median return rate)",
        title  = which === :resources ? "Resources" :
                 which === :consumers ? "Consumers" :
                                         "All Species"
    )

    scatter!(ax, xs, ys;
        # color      = cs,
        colormap   = :viridis,
        # colorrange = extrema(cs),
        markersize = 6,
        alpha      = 0.8
    )

    # Colorbar(fig[1,2], ax.plots[1], label = string(conn_col))
    display(fig)
end

# Example usage:
plot_tau_vs_Rmed(
    df; which = :consumers, conn_col = :conn
)