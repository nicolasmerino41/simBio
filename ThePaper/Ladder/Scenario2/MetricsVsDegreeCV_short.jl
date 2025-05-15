"""
    add_degree_cv!(df::DataFrame)

Ensure `:degree_cv` exists.  Builds unweighted SimpleGraph from the
8th element of `row.p_final` (the A matrix) and computes std(deg)/mean(deg).
"""
function add_degree_cv!(df::DataFrame)
    if :degree_cv ∉ names(df)
        df.degree_cv = [
            begin
                A = row.p_final[8]               # adjacency matrix
                g = SimpleGraph(A .!= 0)         # unweighted graph
                degs = degree(g)
                std(degs) / mean(degs)
            end for row in eachrow(df)
        ]
    end
    return df
end
function add_abundance_cv!(df::DataFrame)
    if :degree_cv ∉ names(df)
        df.abundance_cv = [
            begin
                degs = vcat(row.R_eq, row.C_eq)  # abundance vector
                std(degs) / mean(degs)
            end for row in eachrow(df)
        ]
    end
    return df
end

add_degree_cv!(A)

function plot_vs_degree_cv(
    df::DataFrame,
    variable::Symbol;
    facet_by::Union{Symbol,Nothing}=nothing,
    color_by::Union{Symbol,Nothing}=nothing,
    steps::UnitRange=1:6,
    ncols::Int=4
)
    step_names = [
        "Full Model", "Global A (Global ϵ)", " Global AE",
        "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻"
    ]

    # 1) filter & augment
    df = filter(row -> row.before_persistence_full > 0.0, df)
    add_degree_cv!(df)

    # 2) prepare global color_vals
    N = nrow(df)
    color_vals = if color_by === nothing
        fill(1.0, N)
    else
        Cs = df[!, color_by]
        if eltype(Cs) <: Real
            Cs
        else
            # convert categories to 1,2,3...
            levels = unique(Cs)
            map(x -> findfirst(==(x), levels), Cs)
        end
    end

    # 3) facets
    facet_keys = facet_by === nothing ? [nothing] : unique(df[!, facet_by])

    for f in facet_keys
        # boolean mask for this facet
        mask_f = facet_by === nothing ? trues(N) : (df[!, facet_by] .== f)
        sub   = df[mask_f, :]
        cvals = color_vals[mask_f]

        ns   = length(steps)
        nrows = ceil(Int, ns/ncols)
        fig = Figure(; size=(ncols*300, nrows*300 + (facet_by===nothing ? 0 : 40)))

        if facet_by !== nothing
            Label(fig[1,1:ncols], "$(facet_by) = $(f)"; fontsize=20)
        end

        for (i, step) in enumerate(steps)
            r = div(i-1, ncols) + (facet_by===nothing ? 1 : 2)
            c = mod(i-1, ncols) + 1
            ax = Axis(fig[r, c];
                title   = "$(step_names[step])",
                xlabel  = "Degree CV",
                ylabel  = string(variable, "_S", step),
                titlesize = 12, xlabelsize = 10, ylabelsize = 10
            )

            x = sub.degree_cv
            y = sub[!, Symbol(string(variable), "_S", step)]
            mask_xy = .!ismissing.(x) .& .!ismissing.(y)

            if color_by === nothing
                scatter!(ax, x[mask_xy], y[mask_xy]; color=:blue, markersize=6, alpha=0.7)
            else
                scatter!(ax, x[mask_xy], y[mask_xy];
                         color       = cvals[mask_xy],
                         colormap    = :viridis,
                         colorrange  = (minimum(color_vals), maximum(color_vals)),
                         markersize  = 6,
                         alpha       = 0.7)
            end

            if sum(mask_xy) >= 2
                r_val = cor(x[mask_xy], y[mask_xy])
                text!(ax, "r=$(round(r_val,digits=2))";
                      position = (maximum(x[mask_xy]), minimum(y[mask_xy])),
                      align    = (:right,:bottom),
                      fontsize = 10)
            end
        end

        display(fig)
    end
    return nothing
end

# Persistence vs degree_cv, coloured by connectance, faceted by scenario
plot_vs_degree_cv(
    A, :rt_pulse;
    # facet_by = :epsi,
    color_by = :IS,
)

# Return‐time vs degree_cv, no facets, no colour
step_keys  = ["S1","S2","S3","S4","S5","S6"]
res_cols = Symbol.("resilience_" .* step_keys)
B = filter(row -> all(row[c] < 0 for c in res_cols), A)
    
function plot_vs_abundance_cv(
    df::DataFrame,
    variable::Symbol;
    facet_by::Union{Symbol,Nothing}=nothing,
    color_by::Union{Symbol,Nothing}=nothing,
    steps::UnitRange=1:6,
    ncols::Int=4
)
    # 1) filter & augment
    df = filter(row -> row.before_persistence_full > 0.0, df)
    add_abundance_cv!(df)

    # 2) prepare global color_vals
    N = nrow(df)
    color_vals = if color_by === nothing
        fill(1.0, N)
    else
        Cs = df[!, color_by]
        if eltype(Cs) <: Real
            Cs
        else
            # convert categories to 1,2,3...
            levels = unique(Cs)
            map(x -> findfirst(==(x), levels), Cs)
        end
    end

    # 3) facets
    facet_keys = facet_by === nothing ? [nothing] : unique(df[!, facet_by])

    for f in facet_keys
        # boolean mask for this facet
        mask_f = facet_by === nothing ? trues(N) : (df[!, facet_by] .== f)
        sub   = df[mask_f, :]
        cvals = color_vals[mask_f]

        ns   = length(steps)
        nrows = ceil(Int, ns/ncols)
        fig = Figure(; size=(ncols*300, nrows*300 + (facet_by===nothing ? 0 : 40)))

        if facet_by !== nothing
            Label(fig[1,1:ncols], "$(facet_by) = $(f)"; fontsize=20)
        end

        for (i, step) in enumerate(steps)
            r = div(i-1, ncols) + (facet_by===nothing ? 1 : 2)
            c = mod(i-1, ncols) + 1
            ax = Axis(fig[r, c];
                title   = "step $step",
                xlabel  = "Degree CV",
                ylabel  = string(variable, "_S", step),
                titlesize = 12, xlabelsize = 10, ylabelsize = 10
            )

            x = sub.abundance_cv
            y = sub[!, Symbol(string(variable), "_S", step)]
            mask_xy = .!ismissing.(x) .& .!ismissing.(y)

            if color_by === nothing
                scatter!(ax, x[mask_xy], y[mask_xy]; color=:blue, markersize=6, alpha=0.7)
            else
                scatter!(ax, x[mask_xy], y[mask_xy];
                         color       = cvals[mask_xy],
                         colormap    = :viridis,
                         colorrange  = (minimum(color_vals), maximum(color_vals)),
                         markersize  = 6,
                         alpha       = 0.7)
            end

            if sum(mask_xy) >= 2
                r_val = cor(x[mask_xy], y[mask_xy])
                text!(ax, "r=$(round(r_val,digits=2))";
                      position = (maximum(x[mask_xy]), minimum(y[mask_xy])),
                      align    = (:right,:bottom),
                      fontsize = 10)
            end
        end

        display(fig)
    end
    return nothing
end

# Persistence vs degree_cv, coloured by connectance, faceted by scenario
plot_vs_abundance_cv(
    A, :reactivity;
    # facet_by = :epsi,
    # color_by = :conn,
)

# Return‐time vs degree_cv, no facets, no colour
step_keys  = ["S1","S2","S3","S4","S5","S6"]
res_cols = Symbol.("resilience_" .* step_keys)
B = filter(row -> all(row[c] < 0 for c in res_cols), A)
    
plot_vs_degree_cv(
    B, :reactivity;
    facet_by = nothing,
    color_by = nothing,
)
