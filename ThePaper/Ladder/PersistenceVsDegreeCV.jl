using CairoMakie, Graphs, Statistics

function add_degree_cv!(df::DataFrame)
    if :degree_cv ∉ names(df)
        df.degree_cv = [ begin
            A = row.p_final[8]                     # extract A from your saved params tuple
            g = SimpleGraph(A .!= 0)              # unweighted graph
            degs = degree(g)
            std(degs) / mean(degs)
        end for row in eachrow(df) ]
    end
    return df
end

function plot_persistence_vs_degree_cv(
    df::DataFrame;
    facet_by::Symbol,
    color_by::Symbol,
    steps::UnitRange=1:16,
    ncols::Int=4
)   
    df = filter(row -> row.pers_full > 0.1, df)
    df = add_degree_cv!(df)
    facets = unique(df[!, facet_by])
    cats   = unique(df[!, color_by])
    palette = distinguishable_colors(length(cats))
    color_map = Dict(cats[i] => palette[i] for i in eachindex(cats))

    for f in facets
        sub = df[df[!, facet_by] .== f, :]
        nsteps = length(steps)
        nrows  = ceil(Int, nsteps/ncols)
        fig = Figure(; size=(900, 600))
        Label(fig[0, 2], "$facet_by = $f", fontsize = 24, tellwidth = false)

        for (i, step) in enumerate(steps)
            row = 1 + div(i-1, ncols)
            col = 1 + mod(i-1, ncols)
            ax = Axis(fig[row, col];
                      title  = "step $step",
                      xlabel = "Degree CV",
                      ylabel = "Persistence")

            xs = sub.degree_cv
            ys = sub[!, Symbol("pers_step_$step")]

            # now build a color array by indexing into our Dict
            cols = [ color_map[val] for val in sub[!, color_by] ]
            mask = .!ismissing.(xs) .& .!ismissing.(ys) .& .!isnan.(xs) .& .!isnan.(ys)
            scatter!(ax, xs[mask], ys[mask]; color=cols, markersize=6, alpha=0.7)

            bad_idx = isnan.(xs) .| isnan.(ys)
            correlation = cor(xs[.!bad_idx], ys[.!bad_idx])
            text!(ax, "r=$(round(correlation, digits=2))",;
                    position = (1.9,0.3), 
                    align = (:right,:bottom), 
                    fontsize=8)

        end

        display(fig)
    end
end

subset_df = filter(row -> row.S == 50, dfp)
plot_persistence_vs_degree_cv(
    dfp;
    facet_by = :C,
    color_by = :scenario,
    steps    = 1:16,
    ncols    = 4
)