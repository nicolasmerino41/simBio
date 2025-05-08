function estimate_powerlaw_exponent(degrees::Vector{Int}; k_min::Int=1)
    degrees = filter(d -> d >= k_min, degrees)
    n = length(degrees)
    if n == 0
        return NaN
    end
    α = 1 + n / sum(log.(degrees ./ (k_min - 0.5)))
    return α
end

# ------------------ MAIN SCRIPT ------------------
begin
    # Parameters
    pareto_exponent = 2.0
    scenario = :ER
    S = 200
    R = 20
    conn = 0.1

    # Generate A matrix and degrees
    A = zeros(S, S)
    A = make_A(A, R, conn, scenario; pareto_exponent=pareto_exponent)
    graph = SimpleGraph(A .!= 0)
    degs = degree(graph)
    println(degs)
    # Compute statistics
    alpha_est = estimate_powerlaw_exponent(degs, k_min=1)
    degree_cv = std(degs) / mean(degs)

    println("For pareto_exponent = $pareto_exponent:")
    println("    Estimated power-law exponent α = $(round(alpha_est, digits=3))")
    println("    Degree coefficient of variation CV = $(round(degree_cv, digits=3))")

    # Plot degree distribution
    fig = Figure(; size = (1000, 500))
    ax = Axis(fig[1,1], xlabel = "Degree", ylabel = "Frequency", title = "Degree Distribution (PL)")
    hist!(ax, degs; bins=0:maximum(degs), normalization=:pdf)

    # Plot fitted power-law on top
    kmin = minimum(degs)
    k_fit = range(kmin, maximum(degs); length=100)
    pdf_fit = (alpha_est - 1) * (kmin)^(alpha_est - 1) * k_fit.^(-alpha_est)
    lines!(ax, k_fit, pdf_fit; linewidth=2, linestyle=:dash, label="Power-law fit", color = :red)

    axislegend(ax)
    fig
end
