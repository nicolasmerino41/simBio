df = deserialize("ThePaper/Ladder/Outputs/NonFeasibility/R960.jls")

begin
    fig = Figure(; size = (1000, 460))
    for (idx, row) in enumerate(sample(eachrow(df), 5))
    abund = vcat(row.R_eq, row.C_eq)
    nonzero_abund = abund[abund .> 0]
    zero_abund = abund[abund .== 0]
    col = (idx - 1) % 3 + 1
    row_pos = (idx - 1) ÷ 3 + 1
    ax = Axis(fig[row_pos, col],
              title = "Abundance Histogram $(idx)",
              xlabel = "Abundance",
              ylabel = "Frequency")

    # Plot non-zero abundances
    hist!(ax, nonzero_abund, bins=30, color=:blue, strokewidth=0.5, label="Non-zero abundances")
    # Plot zero abundances as bars or annotations
    if !isempty(zero_abund)
        barplot!(ax, [0], [length(zero_abund)], color=:red, width=maximum(nonzero_abund)/30, label="Zero abundances")
    end
    axislegend(ax)
end
    display(fig)
end


