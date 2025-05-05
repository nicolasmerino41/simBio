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
    steps::UnitRange=1:16,
    ncols::Int=4
)   
    df = filter(row -> row.pers_full > 0, df)
    df = add_degree_cv!(df)
    facets = unique(df[!, facet_by])
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

            mask = .!ismissing.(xs) .& .!ismissing.(ys) .& .!isnan.(xs) .& .!isnan.(ys)
            scatter!(ax, xs[mask], ys[mask]; markersize=6, alpha=0.7)

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

plot_persistence_vs_degree_cv(
    dfp;
    facet_by = :C,
    steps    = 1:16,
    ncols    = 4
)