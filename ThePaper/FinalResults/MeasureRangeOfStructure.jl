using Graphs, StatsBase, DataFrames, Statistics

# ————————————————————————————————
# Power-law slope estimation (linear regression on log-log)
# ————————————————————————————————
function powerlaw_slope(degrees::Vector{Int})
    d = degrees[degrees .> 0]
    if length(d) < 10
        return NaN
    end
    counts = countmap(d)
    k_vals = collect(keys(counts))
    pk_vals = [counts[k] / length(d) for k in k_vals]
    log_k = log10.(k_vals)
    log_pk = log10.(pk_vals)
    X = hcat(ones(length(log_k)), log_k)
    β = X \ log_pk
    return β[2]
end

# ————————————————————————————————
# Gini coefficient
# ————————————————————————————————
function gini(x::AbstractVector{<:Real})
    n = length(x)
    if n == 0
        return NaN
    end
    sorted = sort(abs.(x))
    index = collect(1:n)
    return (2 * sum(index .* sorted) - (n + 1) * sum(sorted)) / (n * sum(sorted))
end

# ————————————————————————————————
# Main structure analysis function
# ————————————————————————————————
function analyze_matrix_structures(df::DataFrame)
    results = DataFrame(
        id = Int[],
        mean_degree_resource = Float64[],
        mean_degree_consumer = Float64[],
        slope_resource = Float64[],
        slope_consumer = Float64[],
        gini_total = Float64[],
        degree_cv_total = Float64[]  # ← NEW
    )

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]  # interaction matrix
        Adj = A .!= 0.0

        k_res = sum(Adj[1:30, :], dims=2)[:]
        k_con = sum(Adj[31:end, :], dims=2)[:]
        degs_all = vcat(k_res, k_con)

        slope_res = powerlaw_slope(k_res)
        slope_con = powerlaw_slope(k_con)
        if mean(k_res) < 1 || mean(k_con) < 1
            println("Skipping ID=$i: deg_res=$(mean(k_res)), deg_con=$(mean(k_con))")
            continue
        end
        g = gini(degs_all)
        cv = std(degs_all) / mean(degs_all)

        push!(results, (
            i,
            mean(k_res),
            mean(k_con),
            slope_res,
            slope_con,
            g,
            cv
        ))
    end

    return results
end

# ————————————————————————————————
# Summary reporting
# ————————————————————————————————
function summarize_structure_ranges(results::DataFrame)
    println("=== STRUCTURAL SUMMARY ACROSS MATRICES ===\n")

    metrics = [
        :mean_degree_resource,
        :mean_degree_consumer,
        :slope_resource,
        :slope_consumer,
        :gini_total,
        :degree_cv_total   # ← NEW
    ]

    for m in metrics
        col = results[!, m]
        col_clean = skipmissing(col)
        println("• $(m):")
        println("   range = $(round(minimum(col_clean), digits=3)) → $(round(maximum(col_clean), digits=3))")
        println("   mean  = $(round(mean(col_clean), digits=3)), std = $(round(std(col_clean), digits=3))\n")
    end

    sorted_by_slope = sort(results, :slope_consumer)
    println("Matrix with lowest consumer slope (most scale-free): ID = $(sorted_by_slope.id[1]), slope = $(round(sorted_by_slope.slope_consumer[1], digits=3))")
    println("Matrix with highest consumer slope (least scale-free): ID = $(sorted_by_slope.id[end]), slope = $(round(sorted_by_slope.slope_consumer[end], digits=3))\n")

    sorted_by_gini = sort(results, :gini_total, rev=true)
    println("Matrix with highest degree inequality (Gini): ID = $(sorted_by_gini.id[1]), Gini = $(round(sorted_by_gini.gini_total[1], digits=3))")
    println("Matrix with lowest degree inequality (Gini): ID = $(sorted_by_gini.id[end]), Gini = $(round(sorted_by_gini.gini_total[end], digits=3))")

    sorted_by_cv = sort(results, :degree_cv_total, rev=true)
    println("Matrix with highest degree CV: ID = $(sorted_by_cv.id[1]), CV = $(round(sorted_by_cv.degree_cv_total[1], digits=3))")
    println("Matrix with lowest degree CV: ID = $(sorted_by_cv.id[end]), CV = $(round(sorted_by_cv.degree_cv_total[end], digits=3))")

    return nothing
end

res = analyze_matrix_structures(G)
summarize_structure_ranges(res)


function plot_random_degree_distributions_makie(df::DataFrame; n::Int = 10, seed::Int = 42)
    # Random.seed!(seed)
    chosen = rand(1:nrow(df), n)

    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)

    fig = Figure(; size=(1100, 450))
    for (i, idx) in enumerate(chosen)
        A = df.p_final[idx][2]             # interaction matrix
        Adj = A .!= 0.0
        degrees = vec(sum(Adj, dims=2))[31:end]
        sorted_deg = sort(degrees; rev=true)

        ax = Axis(fig[(i - 1) ÷ ncols + 1, (i - 1) % ncols + 1],
            title="Matrix ID = $(idx)",
            xlabel="Node rank", ylabel="Degree",
            xticklabelsize=9, yticklabelsize=9,
            titlesize=10
        )
        lines!(ax, 1:length(sorted_deg), sorted_deg)
    end

    display(fig)
    return fig
end

plot_random_degree_distributions_makie(G)
