# Full clean version:
using Random, Distributions, StatsBase, Graphs, Statistics, CairoMakie

# ------------------------------
# Step 1: Generate a better synthetic power-law network
# ------------------------------
function generate_powerlaw_degrees(N::Int, pareto_exponent::Float64; k_min::Int=1, k_max::Int=1000)
    # Use inverse transform sampling for clean power-law
    u = rand(N)
    degs = floor.(Int, k_min * (1 .- u).^( -1/(pareto_exponent-1) ))
    # Clip to k_max
    degs = clamp.(degs, k_min, k_max)
    return degs
end

# ------------------------------
# Step 2: Proper MLE estimator (Clauset-Newman)
# ------------------------------
function estimate_powerlaw_exponent_MLE(degrees::Vector{Int}; k_min::Int=1)
    filtered_degrees = degrees[degrees .>= k_min]
    n = length(filtered_degrees)
    if n == 0
        error("No degrees >= k_min")
    end
    alpha = 1 + n / sum(log.(filtered_degrees ./ (k_min - 0.5)))
    return alpha
end

# ------------------------------
# Step 3: Plot PDF and CCDF
# ------------------------------
function plot_degree_distribution(degrees::Vector{Int}, alpha::Float64; k_min::Int=1)
    fig = Figure(resolution=(1400, 600))

    # Plot PDF
    ax1 = Axis(fig[1, 1], xlabel="Degree", ylabel="P(k)", yscale=log10, xscale=log10, title="PDF of Degree")
    hist!(ax1, degrees; bins=100, normalization=:pdf)
    k_vals = range(k_min, stop=maximum(degrees), length=100)
    pdf_fit = (alpha-1) * k_min^(alpha-1) * k_vals.^(-alpha)
    lines!(ax1, k_vals, pdf_fit, linewidth=2, label="Fitted power law")

    # Plot CCDF
    ax2 = Axis(fig[1, 2], xlabel="Degree", ylabel="P(K>k)", yscale=log10, xscale=log10, title="CCDF of Degree")
    sorted_degrees = sort(degrees)
    n = length(sorted_degrees)
    ccdf_y = reverse((1:n) ./ n)
    lines!(ax2, sorted_degrees, ccdf_y, label="Empirical CCDF")

    fig
end

# ------------------------------
# Step 4: MAIN PROGRAM
# ------------------------------
function main(pareto_exponent)
    Random.seed!(1234)

    S = 5000  # much larger for a better power law
    println("Generating degrees...")

    degrees = generate_powerlaw_degrees(S, pareto_exponent; k_min=1, k_max=500)

    println("Estimating alpha...")
    alpha_est = estimate_powerlaw_exponent_MLE(degrees; k_min=1)

    degree_cv = std(degrees) / mean(degrees)

    println("\nSummary:")
    println("True Pareto exponent used to generate = $pareto_exponent")
    println("Estimated power-law exponent (alpha) = $alpha_est")
    println("Degree coefficient of variation (CV) = $degree_cv")

    fig = plot_degree_distribution(degrees, alpha_est; k_min=1)
    display(fig)
end

# Run
main(1.5)