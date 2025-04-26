function plot_x_vs_y(
    df::DataFrame, x_value, y_value;
    color_by = :conn,
    variable_to_subset = nothing, subset = nothing
)
    step_keys = ["Full", "S1", "S2", "S3", "S4", "S5", "S6", "S7",
                 "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15"]

    titles = [
        "Full Model", "Step 1", "Step 2", "Step 3", "Step 4",
        "Step 5", "Step 6", "Step 7", "Step 8", "Step 9",
        "Step 10", "Step 11", "Step 12", "Step 13", "Step 14", "Step 15"
    ]

    if !isnothing(variable_to_subset) && !isnothing(subset)
        df = df[df[!, variable_to_subset] .== subset, :]
    elseif (!isnothing(variable_to_subset) && isnothing(subset)) || (isnothing(variable_to_subset) && !isnothing(subset))
        @error "Error: variable_to_subset is specified but subset is not"
    end

    Cs = df[!, color_by]
    is_numeric_color = eltype(Cs) <: Real
    color_vals = is_numeric_color ? Cs : map(x -> findfirst(==(x), unique(Cs)), Cs)

    fig = Figure(; size = (1600, 700))
    ncols = 4
    nrows = ceil(Int, length(step_keys) / ncols)

    label_text = isnothing(subset) ?
        "$y_value vs $x_value (colored by $color_by)" :
        "$y_value vs $x_value when $variable_to_subset = $subset (colored by $color_by)"
    Label(fig[0, :], label_text, fontsize = 24, tellwidth = false)

    for (i, step) in enumerate(step_keys)
        r, c = fldmod1(i, ncols)
        ax = Axis(
            fig[r, c],
            xlabel = String(x_value),
            ylabel = String(y_value),
            title = titles[i]
        )

        x_key = hasproperty(df, Symbol("$(x_value)_$step")) ? Symbol("$(x_value)_$step") : x_value
        y_key = Symbol("$(y_value)_$step")

        x_raw = df[:, x_key]
        y_raw = df[:, y_key]

        if y_value == :aft
            limits!(ax, (minimum(x_raw), maximum(x_raw)), (0.0, 1.05))
        end

        valid = [isa(x, Real) && isa(y, Real) && !ismissing(x) && !ismissing(y)
                 for (x, y) in zip(x_raw, y_raw)]

        x_clean = [Float64(x) for (x, v) in zip(x_raw, valid) if v]
        y_clean = [Float64(y) for (y, v) in zip(y_raw, valid) if v]
        color_clean = [Float64(c) for (c, v) in zip(color_vals, valid) if v && isa(c, Real) && !ismissing(c)]

        scatter!(ax, x_clean, y_clean;
            markersize = 6, color = color_clean,
            colormap = :viridis,
            colorrange = (minimum(color_vals), maximum(color_vals)))

        if length(x_clean) > 2
            Xmat = [ones(length(x_clean)) x_clean]
            beta = Xmat \ y_clean

            x_fit = range(minimum(x_clean), maximum(x_clean), length = 100)
            y_fit = beta[1] .+ beta[2] .* x_fit
            lines!(ax, x_fit, y_fit; color = :black, linestyle = :dash, linewidth = 1.5)

            r_val = cor(x_clean, y_clean)
            slope_val = beta[2]

            x_text = maximum(x_clean) - 0.3 * (maximum(x_clean) - minimum(x_clean))
            y_text = maximum(y_clean) - 0.2 * (maximum(y_clean) - minimum(y_clean))

            text!(ax, "r = $(round(r_val, digits = 2)), slope = $(round(slope_val, digits = 2))",
                position = (x_text, y_text), fontsize = 13)
        end
    end

    display(fig)
end

# Example usage:
plot_x_vs_y(df_results, :degree_cv, :per; color_by=:conn)
plot_x_vs_y(df_results, :res, :scorr; color_by=:scenario, variable_to_subset=:S, subset=50)
# fig = plot_x_vs_y(df_results, :res, :rt; color_by=:conn)
for i in 0.1:0.1:0.5
    fig = plot_x_vs_y(df_results, :degree_cv, :res; color_by=:conn, variable_to_subset = :C_ratio, subset = i)
end

############## FOR WHEN MAC SPITS OUT THE RESULTS ##############
mac_results = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL.csv") |> DataFrame
mac_results1 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75.csv") |> DataFrame
mac_results2 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint.csv") |> DataFrame
df_to_plot = mac_results2
name_of_the_x = [:conn, :C_ratio, :degree_cv, :RelVar, :IS, :d, :m]
name_of_the_y = [:aft, :res, :scorr, :rt, :rea]
name_of_the_col = [:conn, :C_ratio, :IS, :degree_cv, :RelVar, :scenario, :d, :m]

for i in name_of_the_x, j in name_of_the_y, k in name_of_the_col
    if i != j && i != k && j != k
        plot_x_vs_y(df_to_plot, i, j; color_by=k)
    end
end
plot_x_vs_y(df_to_plot, :degree_cv, :aft; color_by=:conn)
subs = filter(row -> row.scenario == "PL", df_to_plot)
