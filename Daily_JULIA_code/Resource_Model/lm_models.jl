using DataFrames, GLM

# Join the species effects DataFrame with TrophInd.
# Here, we assume that `see_even` has a column :species_name and TrophInd has columns :Species and :TL.
df_join = innerjoin(see_even, TrophInd, on = :species_name => :Species)
avg_species_metrics
# Fit a linear model with average_effect as the response and TL as the predictor.
lm_model = lm(@formula(average_effect ~ TL), df_join)

# Print the model summary (coefficients, p-values, etc.)
println("Linear Model Summary:")
println(coeftable(lm_model))

using DataFrames, GLM, Statistics

# Join the two DataFrames on species name.
# Here, we assume that see_even has a column :species_name and :average_effect,
# and that avg_species_metrics has :species_name along with network metrics.
df_joined = innerjoin(see_even, avg_species_metrics, on = :species_name)

# Fit different models to explore relationships.
models = Dict()

# Model 1: Using only mean_indegree.
models["Model 1: average_effect ~ mean_indegree"] =
    lm(@formula(average_effect ~ mean_indegree), df_joined)

# Model 2: Using only mean_total_degree.
models["Model 2: average_effect ~ mean_total_degree"] =
    lm(@formula(average_effect ~ mean_total_degree), df_joined)

# Model 3: Using mean_betweenness and mean_closeness.
models["Model 3: average_effect ~ mean_betweenness + mean_closeness"] =
    lm(@formula(average_effect ~ mean_betweenness + mean_closeness), df_joined)

# Model 4: Using only mean_clustering.
models["Model 4: average_effect ~ mean_clustering"] =
    lm(@formula(average_effect ~ mean_clustering), df_joined)

# Model 5: Using all metrics additively.
models["Model 5: average_effect ~ mean_indegree + mean_total_degree + mean_betweenness + mean_closeness + mean_clustering"] =
    lm(@formula(average_effect ~ mean_indegree + mean_total_degree + mean_betweenness + mean_closeness + mean_clustering), df_joined)

# Model 6: Interaction model (all pairwise interactions among predictors).
models["Model 6: Interaction model"] =
    lm(@formula(average_effect ~ (mean_indegree + mean_total_degree + mean_betweenness + mean_closeness + mean_clustering)^2), df_joined)

# Print out summaries for each model.
for (name, model) in models
    println("---- ", name, " ----")
    println(coeftable(model))
    println("AIC: ", aic(model))
    println("")
end

