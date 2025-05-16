# (Descriptor functions assumed to be defined above)
# assume add_all_descriptors! and add_biomass_aggregates! defined above

# -------------------------------------------------------------------
# 2) Species‐level stats: degree vs neighbor biomass
# -------------------------------------------------------------------
"""
compute_species_neighbor_stats(df) -> DataFrame

Returns a DataFrame with one row per species per community,
containing:
  :community        community index (row number)
  :species          global species index (1..S)
  :type             :consumer or :resource
  :degree           number of links (prey for consumers, predators for resources)
  :neighbor_biomass total biomass of those linked species at equilibrium
"""
function compute_species_neighbor_stats(df::DataFrame)
    stats = DataFrame(
        community        = Int[],
        species          = Int[],
        type             = Symbol[],
        degree           = Int[],
        neighbor_biomass = Float64[]
    )
    for (idx, row) in enumerate(eachrow(df))
        A    = row.p_final[8]
        R_eq = row.R_eq
        C_eq = row.C_eq
        R    = length(R_eq)
        C    = length(C_eq)
        S    = R + C

        for i in 1:S
            if i > R
                # consumer i: row i, positive entries are prey
                prey_cols = [j for j in 1:R if A[i,j] > 0]
                deg = length(prey_cols)
                nb  = sum(R_eq[j] for j in prey_cols)
                push!(stats, (idx, i, :consumer, deg, nb))
            else
                # resource i: row i, negative entries are predators
                pred_cols = [j for j in (R+1):S if A[i,j] < 0]
                deg = length(pred_cols)
                nb  = sum(C_eq[j - R] for j in pred_cols)
                push!(stats, (idx, i, :resource, deg, nb))
            end
        end
    end
    return stats
end

"""
plot_degree_vs_neighbor_biomass(df; color_by=nothing)

Creates a two‐panel figure:
  1) Consumer degree vs. total prey biomass
  2) Resource degree vs. total predator biomass
Optionally colors points by any df column (e.g. :connectance).
"""
function plot_degree_vs_neighbor_biomass(
    df::DataFrame;
    color_by::Union{Symbol,Nothing}=nothing
)
    stats = compute_species_neighbor_stats(df)
    # split into consumers and resources
    cons = filter(row -> row.type == :consumer, stats)
    ress = filter(row -> row.type == :resource, stats)

    fig = Figure(resolution = (1200, 400))
    ax1 = Axis(fig[1, 1];
        title = "Consumer Degree vs Prey Biomass",
        xlabel = "Consumer Degree",
        ylabel = "Total Prey Biomass"
    )
    ax2 = Axis(fig[1, 2];
        title = "Resource Degree vs Predator Biomass",
        xlabel = "Resource Degree",
        ylabel = "Total Predator Biomass"
    )

    function scatter_panel(ax, substats)
        x = substats.degree
        y = substats.neighbor_biomass
        mask = .!ismissing.(x) .& .!ismissing.(y)
        if color_by === nothing
            scatter!(ax, x[mask], y[mask]; markersize=6, alpha=0.7)
        else
            cvals = df[substats.community, color_by]
            scatter!(ax, x[mask], y[mask];
                     color = cvals[mask],
                     colormap = :viridis,
                     colorrange = (minimum(cvals), maximum(cvals)),
                     markersize = 6, alpha = 0.7)
        end
        # if sum(mask) >= 2
        #     println("Correlation ", ax.title.text, ": ",
        #             round(cor(x[mask], y[mask]), digits=2))
        # end
    end

    scatter_panel(ax1, cons)
    scatter_panel(ax2, ress)

    display(fig)
end

# -------------------------------------------------------------------
# 3) Example usage
# -------------------------------------------------------------------
plot_degree_vs_neighbor_biomass(A; color_by = :IS)

using DataFrames, Statistics, CairoMakie

"""
    compute_species_neighbor_stats(df::DataFrame)

For each community (row of `df`) and each species in it, compute:
- :type            => "consumer" or "resource"
- :degree          => number of prey (for consumers) or predators (for resources)
- :sum_biomass     => total biomass of those neighbors
- :mean_biomass    => average biomass per neighbor
- :cv_biomass      => CV (std/mean) of neighbor biomasses
Also carries over any columns you want to color by (e.g. :conn, :scen).
Assumes `row.p_final[8]` is your signed interaction matrix A,
`row.R_eq` is a Vector of length R, `row.C_eq` of length C.
"""
function compute_species_neighbor_stats(df::DataFrame)
    stats = DataFrame(
        scen         = String[],
        conn         = Float64[],
        type         = String[],
        degree       = Int[],
        sum_biomass  = Float64[],
        mean_biomass = Float64[],
        cv_biomass   = Float64[],
    )

    for row in eachrow(df)
        A    = row.p_final[8]
        R_eq = row.R_eq
        C_eq = row.C_eq
        R    = length(R_eq)
        S    = size(A,1)

        # carry‐over metadata
        scen = row.scen
        conn = row.conn

        # --- consumers: inspect row i for positive entries into columns 1:R
        for k in 1:length(C_eq)
            i = R + k
            prey_idx = [j for j in 1:R   if A[i,j] > 0]
            deg   = length(prey_idx)
            sb    = sum(R_eq[j] for j in prey_idx)
            mb    = deg > 0 ? sb/deg : 0.0
            cb    = (deg > 1 ? std(R_eq[prey_idx]) : 0.0)
            cvb   = mb > 0 ? cb/mb : 0.0
            push!(stats, (scen, conn, "consumer", deg, sb, mb, cvb))
        end

        # --- resources: inspect row i for negative entries into columns R+1:S
        for i in 1:R
            pred_idx = [j for j in (R+1):S if A[i,j] < 0]
            deg   = length(pred_idx)
            sb    = sum(C_eq[j-R] for j in pred_idx)
            mb    = deg > 0 ? sb/deg : 0.0
            cb    = (deg > 1 ? std(C_eq[pred_idx .- R]) : 0.0)
            cvb   = mb > 0 ? cb/mb : 0.0
            push!(stats, (scen, conn, "resource", deg, sb, mb, cvb))
        end
    end

    return stats
end

"""
    plot_neighbor_stats(df::DataFrame; color_by=nothing)

From `df`, computes species‐neighbor stats via `compute_species_neighbor_stats`,
then makes three side‐by‐side scatter grids for:
  1) Total neighbor biomass
  2) Mean neighbor biomass
  3) CV of neighbor biomass

Each grid has two panels: left=consumers, right=resources.
You can optionally `color_by` any column present in the stats table
(e.g. :conn, :scen, etc.).
"""
function plot_neighbor_stats(df::DataFrame; color_by::Union{Symbol,Nothing}=nothing)
    stats = compute_species_neighbor_stats(df)

    # the three metrics we want to plot
    metrics = [
        (:sum_biomass,  "Total Neighbor Biomass"),
        (:mean_biomass, "Average Neighbor Biomass"),
        (:cv_biomass,   "CV of Neighbor Biomass"),
    ]

    for (sym, ylabel) in metrics
        fig = Figure(resolution = ( 900, 400 ))

        for (col, sp_type) in enumerate(["consumer","resource"])
            ax = Axis(fig[1,col],
                title = uppercasefirst(sp_type),
                xlabel = "Degree",
                ylabel = ylabel,
            )

            sub = filter(r -> r.type == sp_type, stats)
            x = sub.degree
            y = sub[!, sym]

            if color_by === nothing
                scatter!(ax, x, y; markersize=6, alpha=0.7)
            else
                c = sub[!, color_by]
                scatter!(ax, x, y;
                    color      = c,
                    colormap   = :viridis,
                    colorrange = (minimum(c), maximum(c)),
                    markersize = 6,
                    alpha      = 0.7,
                )
            end

            if length(x) > 1
                r = cor(x, y)
                text!(ax, "r=$(round(r, digits=2))",
                    position = (maximum(x), maximum(y)),
                    align    = (:right, :top),
                    fontsize = 12,
                )
            end
        end

        display(fig)
    end

    return nothing
end

# ────────────────────────────────────────────────────────────────────────────────
# Example usage:
#   plot_neighbor_stats(A; color_by=:conn)
# Where A is your original DataFrame with columns:
#   - p_final (where p_final[8] is the signed matrix A)
#   - R_eq, C_eq (vectors of equilibrium biomasses)
#   - conn, scen, etc. for coloring / faceting
# ────────────────────────────────────────────────────────────────────────────────

plot_neighbor_stats(A; color_by=:conn)
