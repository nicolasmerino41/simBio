df = deserialize("ThePaper/Ladder/Outputs/random.jls")

using DataFrames, CairoMakie, Statistics, LinearAlgebra

# Assume R has columns :J_full (Matrix{Float64}) and :d_vec (Vector{Float64})
using DataFrames, CairoMakie, Statistics, LinearAlgebra

function plot_empirical_sigma_ratio(R::DataFrame)
    # theoretical grid
    ratios = 10 .^ range(-2, 1; length=15)

    emp_ratios = Float64[]

    for row in eachrow(R)
        J = row.J_full                # full Jacobian from DataFrame
        S = size(J, 1)

        # Extract self-regulation rates: J = D * M, Mii = -1 ⇒ Dii = -Jii
        d = -diag(J)
        min_d = minimum(d)

        # Recover A: J = D(-I + A) ⇒ -D + D*A = J ⇒ A = I + D^{-1}J
        A = I(S) + Diagonal(1 ./ d) * J

        # Estimate σ as std of off-diagonal entries of A
        offs = [A[i,j] for i in 1:S, j in 1:S if i != j]
        σ_est = std(offs)

        push!(emp_ratios, σ_est / min_d)
    end

    # Plot histogram of empirical ratios
    fig = Figure(; size = (600, 400))
    ax = Axis(fig[1,1];
        xlabel = "σ / min(d)",
        ylabel = "Count",
        title = "Empirical σ/min(d) Ratios"
    )

    # Use bins at the theoretical grid points
    bins = vcat(ratios, maximum(ratios)*1.2)
    hist!(ax, emp_ratios; bins = bins, color = :skyblue)

    # Overlay vertical lines at theoretical ratios
    for r in ratios
        vlines!(ax, [r]; color = :black, linestyle = :dash)
    end

    display(fig)
end

# Example usage:
# Assuming your DataFrame is named `R` and has column :J_full:
plot_empirical_sigma_ratio(df)
