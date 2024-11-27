results = CSV.File("results.csv") |> DataFrame
results_pse_new = CSV.File("results_pse_new.csv") |> DataFrame
results_pse_new_long = CSV.File("results_pse_new_long.csv") |> DataFrame
results_ose_short = CSV.File("results_ose_short.csv") |> DataFrame
results_ose_long = CSV.File("results_ose_long.csv") |> DataFrame
using Plots

PL.scatter(
    results_pse_new_long.mu, results_pse_new_long.total_biomass,
    xlabel="mu (Competition Coefficient)",
    ylabel="Total Biomass",
    title="Relationship between mu and Total Biomass",
    legend=false,
    grid=true
)

using StatsPlots
using CairoMakie
pl = MK.heatmap(
    results_pse_new_long.mu, results_pse_new_long.NPP, results_pse_new_long.total_biomass,
    # xlabel="mu (Competition Coefficient)",
    # ylabel="NPP (Net Primary Productivity)"
    # title="Total Biomass Heatmap",
    # colormap=:thermal
)

results_ose_long = CSV.File("Daily JULIA code/Resource Model/simulation_results.csv") |> DataFrame

fig = Figure(resolution = (500, 400))

ax = Axis(fig[1, 1], xlabel = "Connectivity", ylabel = "Total Biomass",
             title = "Scatter Plot Colored by Mu")

scatter!(
    ax,
    results_ose_long.connectivity,
    results_ose_long.total_biomass,
    color = results_ose_long.mu,
    colormap = :viridis,  # Choose a colormap
    markersize = 8,
    label = "Biomass vs Connectivity"
)