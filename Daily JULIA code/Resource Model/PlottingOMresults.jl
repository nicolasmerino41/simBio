results = CSV.File("results.csv") |> DataFrame
results_pse_new = CSV.File("results_pse_new.csv") |> DataFrame
df = CSV.File("results_pse_new.csv") |> DataFrame

using Plots

PL.scatter(
    df.mu, df.total_biomass,
    xlabel="mu (Competition Coefficient)",
    ylabel="Total Biomass",
    title="Relationship between mu and Total Biomass",
    legend=false,
    grid=true
)

using StatsPlots
using CairoMakie
pl = MK.heatmap(
    df.mu, df.NPP, df.total_biomass,
    xlabel="mu (Competition Coefficient)",
    ylabel="NPP (Net Primary Productivity)",
    title="Total Biomass Heatmap",
    colormap=:thermal
)
