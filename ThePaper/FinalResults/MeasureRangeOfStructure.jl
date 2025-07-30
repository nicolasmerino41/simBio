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
using DataFrames, Statistics, StatsBase, LinearAlgebra

function analyze_matrix_structures(df::DataFrame)
    results = DataFrame(
        id = Int[],
        mean_degree_resource = Float64[],
        mean_degree_consumer = Float64[],
        slope_resource = Float64[],
        slope_consumer = Float64[],
        gini_total = Float64[],
        degree_cv_total = Float64[],
        connectance = Float64[],
        modularity = Float64[],
        nestedness = Float64[],
        # clustering = Float64[]
    )

    for (i, row) in enumerate(eachrow(df))
        A = row.p_final[2]
        Adj = A .!= 0.0
        S = size(Adj, 1)

        # Degrees
        k_res = sum(Adj[1:30, :], dims=2)[:]
        k_con = sum(Adj[31:end, :], dims=2)[:]
        degs_all = vcat(k_res, k_con)

        g = SimpleGraph(A .!= 0)         # unweighted graph
        degs = degree(g)
        degree_cv = std(degs) / mean(degs)

        # if mean(k_res) < 1 || mean(k_con) < 1
        #     println("Skipping ID=$i: deg_res=$(mean(k_res)), deg_con=$(mean(k_con))")
        #     continue
        # end

        # Slopes
        slope_res = powerlaw_slope(k_res)
        slope_con = powerlaw_slope(k_con)

        # Inequality
        g = gini(degs)
        # cv = std(degs_all) / mean(degs_all)

        # Connectance
        L = count(Adj)
        conn = L / (S^2)

        # ——— spectral modularity ———
        k = sum(Adj, dims=2)[:]
        m = sum(k) / 2
        B = zeros(Float64, S, S)
        for u in 1:S, v in 1:S
            B[u, v] = (Adj[u, v] ? 1.0 : 0.0) - (k[u] * k[v]) / (2m)
        end
        vals, vecs = eigen(Symmetric(B))
        v1 = vecs[:, argmax(vals)]
        s = map(x -> x >= 0 ? 1.0 : -1.0, v1)
        Q = (s' * (B * s)) / (4m)

        # ——— nestedness ——— (mean pairwise overlap)
        nested_sum = 0.0
        nested_count = 0
        for u in 1:S, v in u+1:S
            du = sum(Adj[u, :])
            dv = sum(Adj[v, :])
            denom = max(du, dv)
            if denom > 0
                overlap = sum(Adj[u, :] .& Adj[v, :]) / denom
                nested_sum += overlap
                nested_count += 1
            end
        end
        nestedness = nested_count > 0 ? nested_sum / nested_count : NaN

        # # ——— clustering coefficient ———
        # triangles = 0
        # triplets = 0
        # for u in 1:S
        #     neighbors = findall(Adj[u, :] .| Adj[:, u]')
        #     for i in 1:length(neighbors), j in i+1:length(neighbors)
        #         v, w = neighbors[i], neighbors[j]
        #         if Adj[v, w] || Adj[w, v]
        #             triangles += 1
        #         end
        #         triplets += 1
        #     end
        # end
        # clustering = triplets > 0 ? (triangles / triplets) : NaN

        push!(results, (
            i,
            mean(k_res),
            mean(k_con),
            slope_res,
            slope_con,
            g,
            degree_cv,
            conn,
            Q,
            nestedness,
            # clustering
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
        :degree_cv_total,   # ← NEW
        :connectance,
        :modularity,
        :nestedness
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

    sorted_by_mod = sort(results, :modularity, rev=true)
    println("Matrix with highest modularity: ID = $(sorted_by_mod.id[1]), mod = $(round(sorted_by_mod.modularity[1], digits=3))")
    println("Matrix with lowest modularity: ID = $(sorted_by_mod.id[end]), mod = $(round(sorted_by_mod.modularity[end], digits=3))")

    sorted_by_conn = sort(results, :connectance, rev=true)
    println("Matrix with highest connectance: ID = $(sorted_by_conn.id[1]), conn = $(round(sorted_by_conn.connectance[1], digits=3))")
    println("Matrix with lowest connectance: ID = $(sorted_by_conn.id[end]), conn = $(round(sorted_by_conn.connectance[end], digits=3))")

    sorted_by_nested = sort(results, :nestedness, rev=true)
    println("Matrix with highest nestedness: ID = $(sorted_by_nested.id[1]), nested = $(round(sorted_by_nested.nestedness[1], digits=3))")
    println("Matrix with lowest nestedness: ID = $(sorted_by_nested.id[end]), nested = $(round(sorted_by_nested.nestedness[end], digits=3))")

    # sorted_by_clust = sort(results, :clustering_coefficient, rev=true)
    # println("Matrix with highest clustering coefficient: ID = $(sorted_by_clust.id[1]), clust = $(round(sorted_by_clust.clustering_coefficient[1], digits=3))")
    # println("Matrix with lowest clustering coefficient: ID = $(sorted_by_clust.id[end]), clust = $(round(sorted_by_clust.clustering_coefficient[end], digits=3))")

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

function plot_high_cv_distributions(df::DataFrame, res::DataFrame; top_n::Int = 5)
    top_cv = sort(res, :degree_cv_total, rev=true)[1:top_n, :]
    fig = Figure(; size=(1000, 200 * top_n))
    for (i, row) in enumerate(eachrow(top_cv))
        idx = row.id
        A = df.p_final[idx][2]
        Adj = A .!= 0.0
        degs = vec(sum(Adj, dims=2))
        sorted_deg = sort(degs; rev=true)
        ax = Axis(fig[i, 1], title="ID=$(idx), CV=$(round(row.degree_cv_total, digits=3))",
                  xlabel="Rank", ylabel="Degree", titlesize=11)
        lines!(ax, 1:length(sorted_deg), sorted_deg)
    end
    display(fig)
end

plot_high_cv_distributions(G, res)