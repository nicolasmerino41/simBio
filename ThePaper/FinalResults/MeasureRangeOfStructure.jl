using Graphs, StatsBase, DataFrames, Statistics

R = filter(row -> row.scen == :PL, R)

function powerlaw_slope(degrees::Vector{Int})
    d = degrees[degrees .> 0]  # Remove zeros
    if length(d) < 10
        return NaN
    end

    counts = countmap(d)
    k_vals = collect(keys(counts))
    pk_vals = [counts[k] / length(d) for k in k_vals]

    log_k = log10.(k_vals)
    log_pk = log10.(pk_vals)

    X = hcat(ones(length(log_k)), log_k)
    β = X \ log_pk  # linear regression

    slope = β[2]    # power-law exponent ≈ slope
    return slope
end

function gini(x::AbstractVector{<:Real})
    n = length(x)
    if n == 0
        return NaN
    end
    sorted = sort(abs.(x))
    index = collect(1:n)
    return (2 * sum(index .* sorted) - (n + 1) * sum(sorted)) / (n * sum(sorted))
end

function analyze_matrix_structures(df::DataFrame)
    results = DataFrame(
        id = Int[],
        mean_degree_resource = Float64[],
        mean_degree_consumer = Float64[],
        slope_resource = Float64[],
        slope_consumer = Float64[],
        gini_total = Float64[]
    )

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]
        Adj = A .!= 0.0

        k_res = sum(Adj[1:30, :], dims=2)[:]
        k_con = sum(Adj[31:end, :], dims=2)[:]

        slope_res = powerlaw_slope(k_res)
        slope_con = powerlaw_slope(k_con)

        degs_all = vcat(k_res, k_con)
        g = gini(degs_all)

        push!(results, (
            i,
            mean(k_res),
            mean(k_con),
            slope_res,
            slope_con,
            g
        ))
    end

    return results
end


res = analyze_matrix_structures(G)

function summarize_structure_ranges(results::DataFrame)
    println("=== STRUCTURAL SUMMARY ACROSS MATRICES ===\n")

    metrics = [
        :mean_degree_resource,
        :mean_degree_consumer,
        :slope_resource,
        :slope_consumer,
        :gini_total
    ]

    for m in metrics
        col = results[!, m]
        col_clean = skipmissing(col)
        println("• $(m):")
        println("   range = $(round(minimum(col_clean), digits=3)) → $(round(maximum(col_clean), digits=3))")
        println("   mean  = $(round(mean(col_clean), digits=3)), std = $(round(std(col_clean), digits=3))\n")
    end

    # Identify most and least scale-free-like
    sorted_by_slope = sort(results, :slope_consumer)
    println("Matrix with lowest consumer slope (most scale-free): ID = $(sorted_by_slope.id[1]), slope = $(round(sorted_by_slope.slope_consumer[1], digits=3))")
    println("Matrix with highest consumer slope (least scale-free): ID = $(sorted_by_slope.id[end]), slope = $(round(sorted_by_slope.slope_consumer[end], digits=3))\n")

    sorted_by_gini = sort(results, :gini_total, rev=true)
    println("Matrix with highest degree inequality (Gini): ID = $(sorted_by_gini.id[1]), Gini = $(round(sorted_by_gini.gini_total[1], digits=3))")
    println("Matrix with lowest degree inequality (Gini): ID = $(sorted_by_gini.id[end]), Gini = $(round(sorted_by_gini.gini_total[end], digits=3))")

    return nothing
end

summarize_structure_ranges(res)
