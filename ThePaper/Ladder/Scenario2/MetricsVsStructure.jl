
using DataFrames, Graphs, StatsBase, LinearAlgebra, CairoMakie

# -------------------------------------------------------------------
# 1) Network descriptors
# -------------------------------------------------------------------
function add_degree_cv!(df::DataFrame)
    if :degree_cv ∉ names(df)
        df.degree_cv = [ begin
            A = row.p_final[8]
            g = SimpleGraph(A .!= 0)
            degs = degree(g)
            std(degs) / mean(degs)
        end for row in eachrow(df) ]
    end
    return df
end

function add_connectance!(df::DataFrame)
    if :connectance ∉ names(df)
        df.connectance = [ begin
            A = row.p_final[8]
            S = size(A,1)
            count(!=(0), A) / (S^2)
        end for row in eachrow(df) ]
    end
    return df
end

function add_mean_degree!(df::DataFrame)
    if :mean_degree ∉ names(df)
        df.mean_degree = [ begin
            A = row.p_final[8]
            g = SimpleGraph(A .!= 0)
            mean(degree(g))
        end for row in eachrow(df) ]
    end
    return df
end

function add_strength_heterogeneity!(df::DataFrame)
    if :strength_heterogeneity ∉ names(df)
        df.strength_heterogeneity = [ begin
            A = row.p_final[8]
            vals = abs.(A[A .!= 0])
            std(vals) / mean(vals)
        end for row in eachrow(df) ]
    end
    return df
end

function add_assortativity!(df::DataFrame)
    if :assortativity ∉ names(df)
        df.assortativity = [ begin
            A = row.p_final[8]
            g = SimpleGraph(A .!= 0)
            ds = degree(g)
            uv = [(e.src,e.dst) for e in edges(g)]
            cor([ds[u] for (u,_) in uv], [ds[v] for (_,v) in uv])
        end for row in eachrow(df) ]
    end
    return df
end

function add_clustering!(df::DataFrame)
    if :clustering ∉ names(df)
        df.clustering = [ begin
            A = row.p_final[8]
            g = SimpleGraph(A .!= 0)
            global_clustering_coefficient(g)
        end for row in eachrow(df) ]
    end
    return df
end

function add_all_descriptors!(df::DataFrame)
    add_degree_cv!(df)
    add_connectance!(df)
    add_mean_degree!(df)
    add_strength_heterogeneity!(df)
    add_assortativity!(df)
    add_clustering!(df)
    return df
end

# -------------------------------------------------------------------
# 2) Plotting function
# -------------------------------------------------------------------

# """
#     plot_metrics_vs_structure(df, metric, structure;
#                               facet_by=nothing, color_by=nothing,
#                               steps=1:6, ncols=4, remove_unstable=true)

# Plots each step (S1‒S6) of the chosen `metric` against the chosen
# aggregate `structure` property, with optional faceting and coloring.
# """
# -------------------------------------------------------------------
# 2) Biomass aggregates
# -------------------------------------------------------------------
function add_biomass_aggregates!(df::DataFrame)
    if :total_resource_biomass ∉ names(df)
        df.total_resource_biomass = [ sum(row.R_eq) for row in eachrow(df) ]
        df.total_consumer_biomass = [ sum(row.C_eq) for row in eachrow(df) ]
        df.biomass_CR_ratio       = df.total_consumer_biomass ./ df.total_resource_biomass
        df.cv_resource_biomass    = [ std(row.R_eq)/mean(row.R_eq) for row in eachrow(df) ]
        df.cv_consumer_biomass    = [ std(row.C_eq)/mean(row.C_eq) for row in eachrow(df) ]
    end
    return df
end

# -------------------------------------------------------------------
# 3) Plotting function
# -------------------------------------------------------------------
function plot_metrics_vs_structure(
    df::DataFrame,
    metric::Symbol,
    structure::Symbol;
    facet_by::Union{Symbol,Nothing}=nothing,
    color_by::Union{Symbol,Nothing}=nothing,
    steps::UnitRange=1:6,
    ncols::Int=3,
    remove_unstable::Bool=false
)
    # ensure all descriptors are present
    add_all_descriptors!(df)
    add_biomass_aggregates!(df)

    # optionally filter unstable communities
    if remove_unstable
        res_cols = Symbol.(["resilience_S" * string(s) for s in 1:6])
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
    end

    # determine facets
    facets = facet_by === nothing ? [nothing] : unique(df[!, facet_by])
    step_names = ["Full Model","Global A (Global ϵ)","Global AE",
                  "Randomize m_cons","Randomize ξ","Randomize K_res",
                  "Global A (Global ϵ) Mean B","Global AE Mean B"]

    for f in facets
        sub = facet_by === nothing ? df : filter(r -> r[facet_by] == f, df)
        rows = ceil(Int, length(steps)/ncols)
        fig = Figure(; size = (ncols*300, rows*300 + (facet_by===nothing ? 0 : 40)))
        if facet_by !== nothing
            Label(fig[1,1:ncols], "$(facet_by) = $(f)"; fontsize=18)
        end

        for (idx, step) in enumerate(steps)
            row_i = div(idx-1, ncols) + (facet_by===nothing ? 1 : 2)
            col_i = mod(idx-1, ncols) + 1
            ax = Axis(fig[row_i, col_i];
                      title=step_names[step],
                      xlabel=string(structure),
                      ylabel="" * string(metric) * "_S" * string(step),
                      titlesize=14, xlabelsize=12, ylabelsize=12)

            x = sub[!, structure]
            y = sub[!, Symbol(string(metric), "_S", string(step))]
            mask = .!ismissing.(x) .& .!ismissing.(y)

            if color_by === nothing
                scatter!(ax, x[mask], y[mask]; markersize=6, alpha=0.7)
            else
                cvals = sub[!, color_by]
                scatter!(ax, x[mask], y[mask];
                         color=cvals[mask],
                         colormap=:viridis,
                         colorrange=(minimum(cvals), maximum(cvals)),
                         markersize=6, alpha=0.7)
            end

            if sum(mask) >= 2
                r_val = cor(x[mask], y[mask])
                text!(ax, "r=$(round(r_val,digits=2))";
                      position=(maximum(x[mask]), minimum(y[mask])),
                      align=(:right,:bottom), fontsize=10)
            end
        end
        display(fig)
    end
end

# -------------------------------------------------------------------
# 4) Example usage
# -------------------------------------------------------------------
# Ensure descriptors and biomass aggregates are computed
add_all_descriptors!(A)
add_biomass_aggregates!(A)

# Plot resilience vs. biomass CR ratio, colored by connectance
plot_metrics_vs_structure(A, :Rmed, :mean_tau_full; color_by = :conn, remove_unstable=true)

# -------------------------------------------------------------------
# 3) Example usage
# -------------------------------------------------------------------
for metric in [:reactivity, :resilience, :rt_press, :rt_pulse, :after_persistence, :collectivity]
    for structure in [:degree_cv, :connectance, :mean_degree, :strength_heterogeneity, :assortativity, :clustering]
        plot_metrics_vs_structure(A, metric, structure; color_by=:conn, facet_by=:scen, remove_unstable=true)
    end
end

# for metric in [:reactivity, :resilience, :rt_press, :rt_pulse, :after_persistence, :collectivity]
for metric in [:rt_pulse]
    for biomass in [:total_resource_biomass, :total_consumer_biomass, :biomass_CR_ratio, :cv_resource_biomass, :cv_consumer_biomass]
        plot_metrics_vs_structure(
            A, metric, biomass;
            color_by=:IS,
            # facet_by=:scen,
            remove_unstable=true
        )
    end
end

