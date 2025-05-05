using DataFrames, GLM, StatsModels, OrderedCollections

# assume df is your collected DataFrame, e.g. df_collect_1
df = collectivity_df
# 1) make sure categorical columns are encoded as factors
categorical(df_collect_1, :scenario)
df.B_term = convert.(Bool, df.B_term)  # if necessary

# 3) fit linear models for each response
lm_φ   = lm(@formula(φ ~ 1 + B_term + conn + C_ratio + IS + epsilon + skew + scenario), df)
lm_Cmay= lm(@formula(C_may ~ 1 + B_term + conn + C_ratio + IS + epsilon + skew + scenario), df)
lm_depth  = lm(@formula(depth ~ 1 + B_term + conn + C_ratio + IS + epsilon + skew + scenario), df)
lm_unpred = lm(@formula(unpredict ~ 1 + B_term + conn + C_ratio + IS + epsilon + skew + scenario), df)
lm_niche  = lm(@formula(niche ~ 1 + B_term + conn + C_ratio + IS + epsilon + skew + scenario), df)

# 4) look at summaries
println("Spectral radius (ϕ) model:")
display(coeftable(lm_φ))

println("\nMay complexity (C_may) model:")
display(coeftable(lm_Cmay))

println("\nKnock-out depth model:")
display(coeftable(lm_depth))

println("\nTemporal unpredictability model:")
display(coeftable(lm_unpred))

println("\nBiotic niche contribution model:")
display(coeftable(lm_niche))

###############################################################################
###############################################################################
# analysis_pipeline.jl
df = df_collect_1
using CSV, DataFrames, GLM, MixedModels, CairoMakie, StatsBase

# 1. Load data
# df = CSV.read("collectivity_df.csv", DataFrame)
# categorical!(df, :scenario)

# metrics & steps
metrics = [:φ, :C_may, :depth, :unpredict, :niche]
steps   = sort(unique(df.step))

# 2. Makie plots in begin...end blocks

# 2a. Mean metric vs. ladder step, by scenario
begin
    for m in metrics
        fig = Figure(; size = (800, 600))
        ax  = Axis(
            fig[1,1];
            xlabel = "Ladder Step",
            ylabel = string(m),
            title  = "Mean $(m) vs. Step by Scenario"
        )
        for sc in levels(df.scenario)
            sub = df[df.scenario .== sc, :]
            means = [mean(sub[sub.step .== s, m]) for s in steps]
            lines!(ax, steps, means; label = string(sc))
        end
        axislegend(ax; position = :rt)
        fig
        display(fig)
    end
end

# 2b. Correlation with full (step == 1) for each metric
begin
    for m in metrics
        full_vals = df[df.step .== 1, m]
        corrs = [ cor(df[df.step .== s, m], full_vals) for s in steps ]

        fig = Figure(; size = (800, 600))
        ax  = Axis(
            fig[1,1];
            xlabel = "Ladder Step",
            ylabel = "corr($(m), full)",
            title  = "Correlation with Full for $(m)"
        )
        lines!(ax, steps, corrs;)
        fig
        display(fig)
    end
end

# 3. Step-wise linear models (console output)
for s in steps
    sub = df[df.step .== s, :]
    println("\n=== Step $s ===")
    for m in metrics
        # build a FormulaTerm dynamically:
        response   = Term(m)
        rhs_terms  = Term.(predictors)        # vector of Term(:B_term), Term(:conn), …
        formula    = response ~ sum(rhs_terms)  # Term(m) ~ Term(:B_term)+Term(:conn)+…
        lmres      = lm(formula, sub)
        println("Metric: ", m)
        display(coeftable(lmres))
    end
end

# 4. Mixed-effects models (console output)
predictors = [:B_term, :conn, :C_ratio, :IS, :epsilon, :skew]
# build a vector of Term objects for the predictors + scenario
all_terms = [ Term(p) for p in predictors ] 
push!(all_terms, Term(:scenario))

begin
    for m in metrics
        Δr2 = Float64[]
        lhs = Term(m)  # the left‐hand side of the formula

        # construct the full‐model RHS by summing all terms
        rhs_full = reduce(+, all_terms)
        fullf    = lhs ~ rhs_full

        for s in steps
            sub    = df[df.step .== s, :]
            fullr2 = r2(lm(fullf, sub))

            # now drop each predictor in turn
            for p in predictors
                # make a reduced RHS with p removed
                kept = filter(t -> t != Term(p), all_terms)
                rhs_red = reduce(+, kept)
                redf    = lhs ~ rhs_red

                push!(Δr2, fullr2 - r2(lm(redf, sub)))
            end
        end

        # reshape into (steps × predictors)
        mat = reshape(Δr2, length(predictors), length(steps))'

        # plot
        fig = Figure(; size=(900,400))
        ax  = Axis(fig[1,1];
                  xticks=(1:length(steps), string.(steps)),
                  yticks=(1:length(predictors), string.(predictors)),
                  xlabel="Step", ylabel="Predictor",
                  title="Partial R² for $(m)")
        MK.heatmap!(ax, mat; colormap=:viridis)
        display(fig)
    end
end
