# This script compares your full vs. simplified results step-by-step (correlations & mean errors),
# then runs a tiny regression and scatterplots to see which traits (skew, connectance, C, mortality, etc.)
# drive those deviations.

variable = "aft"  # or "aft"
for step in 1:16
    xs = dfp[!, Symbol("$(variable)_full")]
    ys = dfp[!, Symbol("$(variable)_step_$(step)")]
    corr = cor(xs, ys)
    println("Step $step: correlation = $corr")
end

for step in 1:16
    xs = dfp[!, Symbol("$(variable)_full")]
    ys = dfp[!, Symbol("$(variable)_step_$(step)")]
    error = mean(abs.(xs .- ys))
    println("Step $step: mean absolute error = $error")
end

# for step in 1:15
#     dfp[!, Symbol("deviation_step_$(step)")] = abs.(dfp[!, :pers_full] .- dfp[!, Symbol("pers_step_$(step)")])
# end
# using GLM
dfp = copy(dfp_304050_withAFT)
dfp = filter(row -> row.S == 50, dfp)

dfp[!, :deviation] = abs.(dfp.pers_full .- dfp.pers_step_5)  # example for step 5

IS_vals = [mean(abs.(a[8])) for a in dfp.p_final]  # A matrices per row
dfp[!, :IS_computed] = IS_vals

skew_computed = sum.(dfp.C_eq) ./ sum.(dfp.R_eq)
dfp[!, :skew_computed] = skew_computed

using GLM

model = lm(@formula(deviation ~ C + conn + IS_computed + skew_computed + abundance + delta + d + m), dfp)
println(coeftable(model))

model_simple = lm(@formula(deviation_step_15 ~ conn + skew_computed + abundance + m), subs)
println(coeftable(model_simple))
# for step in 1:15
#     dfp[!, Symbol("deviation_step_$(step)")] = abs.(dfp[!, :pers_full] .- dfp[!, Symbol("pers_step_$(step)")])
# end
# using GLM
begin
    fig = Figure(resolution = (1200, 800))

    ax1 = Axis(fig[1, 1], title="Deviation vs Skew", xlabel="Skew (C_eq / R_eq)", ylabel="Deviation")
    scatter!(ax1, dfp.skew_computed, dfp.deviation)

    ax2 = Axis(fig[1, 2], title="Deviation vs Connectance", xlabel="Connectance", ylabel="Deviation")
    scatter!(ax2, dfp.conn, dfp.deviation)

    ax3 = Axis(fig[2, 1], title="Deviation vs C", xlabel="C", ylabel="Deviation")
    scatter!(ax3, dfp.C, dfp.deviation)

    ax4 = Axis(fig[2, 2], title="Deviation vs Mortality (m)", xlabel="Mortality (m)", ylabel="Deviation")
    scatter!(ax4, dfp.m, dfp.deviation)

    display(fig)
end

# Compute deviation columns if not already done
for step in 1:15
    dfp[!, Symbol("deviation_step_$(step)")] = abs.(dfp.pers_full .- dfp[!, Symbol("pers_step_$(step)")])
end

# Variable to color by — change this to any column in dfp
color_var = dfp.IS_computed  # e.g., m, conn, C, etc.
subs = filter(row -> row.C == 9 && row.S == 50, dfp)
begin
    fig = Figure(; size = (1500, 1000))
    
    for step in 1:15
        row = cld(step, 5)
        col = (step - 1) % 5 + 1
        ax = Axis(fig[row, col], title = "Step $step", xlabel = "Skew", ylabel = "Deviation")

        scatter!(
            ax,
            subs.skew_computed,
            subs[!, Symbol("deviation_step_$(step)")],
            color = color_var,
            colormap = :viridis
        )
    end

    display(fig)
end

