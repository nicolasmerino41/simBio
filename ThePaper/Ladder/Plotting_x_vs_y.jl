function plot_x_vs_yy(
    df::DataFrame, x_value, y_value;
    color_by = :conn,
    variable_to_subset = nothing, subset = nothing
)
    step_keys = ["Full", "S1", "S2", "S3", "S4", "S5", "S6", "S7",
                 "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15"]

    titles = [
        "Full Model",
        "Full A (Species‑specific ϵ)", "Full A (Global ϵ)", "Full A (Re‑randomised ϵ)",
        "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species‑specific ϵ)", "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re‑randomised ϵ)",
        "Global A (ϵ_full)", "Global A (Species‑specific ϵ)", "Global A (Global ϵ)", "Global A (Re‑randomised ϵ)",
        "Re-randomised A (ϵ_full)", "Re-randomised A (Species‑specific ϵ)", "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
    ]

    if !isnothing(variable_to_subset) && !isnothing(subset)
        df = df[df[!, variable_to_subset] .== subset, :]
    elseif (!isnothing(variable_to_subset) && isnothing(subset)) || (isnothing(variable_to_subset) && !isnothing(subset))
        @error "Error: variable_to_subset is specified but subset is not"
    end

    Cs = df[!, color_by]
    is_numeric_color = eltype(Cs) <: Real
    color_vals = is_numeric_color ? Cs : map(x -> findfirst(==(x), unique(Cs)), Cs)

    fig = Figure(; size = (1600, 750))
    ncols = 4
    nrows = ceil(Int, length(step_keys) / ncols)

    label_text = isnothing(subset) ?
        "$y_value vs $x_value (colored by $color_by)" :
        "$y_value vs $x_value when $variable_to_subset = $subset (colored by $color_by)"
    Label(fig[0, 2], label_text, fontsize = 24, tellwidth = false)

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
        # println(minimum(y_raw))
        # println(maximum(y_raw))

        valid = [isa(x, Real) && isa(y, Real) && !ismissing(x) && !ismissing(y) && !isnan(x) && !isnan(y)
                 for (x, y) in zip(x_raw, y_raw)]

        x_clean = [Float64(x) for (x, v) in zip(x_raw, valid) if v]
        y_clean = [Float64(y) for (y, v) in zip(y_raw, valid) if v]

        if y_value == :aft
            limits!(ax, (minimum(x_clean), maximum(x_clean)), (0.5, 1.1))
        end
        # elseif y_value == :rt
        #     limits!(ax, (minimum(x_clean), maximum(x_clean)), (-10, 50.0))
        # end
        if x_value == :degree_cv && y_value == :aft
            limits!(ax, (0, 4), (0.5, 1.05))
        end

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

            x_text = maximum(x_clean) - 0.5 * (maximum(x_clean) - minimum(x_clean))
            y_text = minimum(y_clean) + 0.3 * (maximum(y_clean) - minimum(y_clean))

            text!(ax, "r = $(round(r_val, digits = 2)), slope = $(round(slope_val, digits = 2))",
                position = (x_text, y_text), fontsize = 13)
        end
    end

    display(fig)
end
df_to_plot = df_to_plot[df_to_plot[!, :scenario] .== "MOD", :]
plot_x_vs_yy(df_to_plot, :C, :aft; color_by=:conn, variable_to_subset = :scenario, subset = "PL")
plot_x_vs_yy(df_to_plot, :degree_cv, :rt; color_by=:C)
# for i in 30:10:60
#     plot_x_vs_yy(df_to_plot, :degree_cv, :aft; color_by=:scenario, variable_to_subset = :C_ratio, subset = 0.1)
# end

############## FOR WHEN MAC SPITS OUT THE RESULTS ##############
mac_results = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL.csv") |> DataFrame
mac_results1 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75.csv") |> DataFrame
mac_results2 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint.csv") |> DataFrame
mac_results3 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint_speciesSpecificRTs_s708090100.csv") |> DataFrame
mac_results4 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint_speciesSpecificRTs_s30405060.csv") |> DataFrame
mac_results5 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint_speciesSpecificRTs_s30to100.csv") |> DataFrame
only_mod =  CSV.File("ThePaper/Ladder/Outputs/guided_results_noconstraint_sspRTs_s40506070_onlyMOD.csv") |> DataFrame
df_to_plot = df
yes_const = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERPLMOD_75_yesconstraint_NOsspRTs_s30405060.csv") |> DataFrame
MOD_Results1 = CSV.File("ThePaper/Ladder/Outputs/guided_results_ERplusPL1_75_noconstraint_sspRTs_s40506070_withMOD.csv") |> DataFrame
# name_of_the_x = [:conn, :C_ratio, :degree_cv, :RelVar, :IS, :d, :m]
# name_of_the_y = [:aft, :res, :scorr, :rt, :rea]
# name_of_the_col = [:conn, :C_ratio, :IS, :degree_cv, :RelVar, :scenario, :d, :m]

# for i in name_of_the_x, j in name_of_the_y, k in name_of_the_col
#     if i != j && i != k && j != k
#         plot_x_vs_yy(df_to_plot, i, j; color_by=k)
#     end
# end

mac_results5 = filter(row -> row.C_ratio == 0.1, mac_results5)


as = df[df[!, :scenario] .== :PL, :]