using CairoMakie, Distributions, Random
begin

    # Parameters
    m = 1

    # 1) “Incorrect” log-normal: μ = m, σ = 2m
    dist1 = LogNormal(m, 2m)

    # 2) “Correct” log-normal: μ = log(m) - σ²/2, σ = m
    dist2 = LogNormal(log(m) - m^2/2, m*2)

    # Draw samples
    Random.seed!(42)
    N = 1000000
    samps1 = rand(dist1, N)
    samps2 = rand(dist2, N)

    # Plot with Makie
    fig = Figure(; size = (800, 500))
    ax = Axis(fig[1, 1];
        xlabel = "Value",
        ylabel = "Density",
        title  = "LogNormal parameterizations"
    )

    # # Overlay two histograms, normalized to PDF
    # hist!(ax, samps1;
    #     bins        = 200,
    #     # normalization = :pdf,
    #     color       = (:steelblue, 0.5),
    #     label       = "LogNormal(μ=m, σ=2m)"
    # )
    hist!(ax, samps2;
        bins        = 200,
        # normalization = :pdf,
        color       = (:tomato, 0.5),
        label       = "LogNormal(μ=log(m)-σ²/2, σ=m)"
    )

    # xlimits!(ax, 0, 20)
    # axislegend(ax; position = :rt)

    display(fig)
end