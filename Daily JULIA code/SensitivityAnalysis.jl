using CSV, DataFrames

# Load the results
results_df = CSV.read("drago_result_with_herbivory_new.csv", DataFrame)

# Ensure all relevant columns are present
required_cols = [:sigma, :epsilon, :alfa, :k_DA_name, :sigma_comp,
                 :avg_shannon, :avg_bbp, :richness_similarity,
                 :alive_predators, :mean_tl, :mean_n_of_species, :Threads, :NaNs]

missing_cols = setdiff(required_cols, names(results_df))
if !isempty(missing_cols)
    error("Missing columns in the data: ", join(missing_cols, ", "))
end

# Remove rows with NaNs
clean_df = filter(row -> row.NaNs == 0.0, results_df)

describe(clean_df)

using Plots

histogram(clean_df.sigma, bins=20, xlabel="Sigma", title="Distribution of Sigma")
# Repeat for other parameters and metrics

using Statistics

param_cols = [:sigma, :epsilon, :alfa, :sigma_comp]
metric_cols = [:avg_shannon, :avg_bbp, :richness_similarity,
               :alive_predators, :mean_tl, :mean_n_of_species]
param_cols_names = ["Sigma", "Epsilon", "Alfa", "Sigma_comp"]
metric_cols_names = ["Shannon", "BBP", "Richness", "Predators", "TL", "Species"]
analysis_df = clean_df[:, vcat(param_cols, metric_cols)]
corr_matrix = cor(Matrix(analysis_df))

# Convert to DataFrame for readability
corr_df = DataFrame(corr_matrix, :auto)

using GLM, StatsModels

for metric in metric_cols
    formula = Term(metric) ~ sum(Term.(param_cols))
    model = lm(formula, clean_df)
    println("Regression results for $(metric):")
    display(coeftable(model))
end

# Example with interactions
formula = @formula(avg_shannon ~ (sigma + epsilon + alfa + sigma_comp)^2)
model = lm(formula, clean_df)
coeftable(model)

for col in param_cols
    μ = mean(clean_df[!, col])
    σ = std(clean_df[!, col])
    clean_df[!, col] = (clean_df[!, col] .- μ) ./ σ
end

for param in param_cols
    for metric in metric_cols
        scatter(clean_df[!, param], clean_df[!, metric],
                xlabel=string(param), ylabel=string(metric),
                title="$(metric) vs $(param)")
    end
end

# Example: Effect of sigma and epsilon on avg_shannon
using StatsPlots
x_vals = unique(clean_df.sigma)
y_vals = unique(clean_df.epsilon)
x_vals = sort(x_vals)
y_vals = sort(y_vals)
using DataFrames

# Initialize an empty matrix
z_matrix = zeros(Float64, length(y_vals), length(x_vals))

# Populate the matrix
for i in 1:length(x_vals)
    for j in 1:length(y_vals)
        # Filter the DataFrame for the specific sigma and epsilon
        subset = clean_df[(clean_df.sigma .== x_vals[i]) .& (clean_df.epsilon .== y_vals[j]), :]
        if nrow(subset) > 0
            # Assuming there's only one avg_shannon value per combination
            z_matrix[j, i] = subset.avg_shannon[1]
        else
            # Assign NaN or a placeholder if no data is available
            z_matrix[j, i] = NaN
        end
    end
end

using Plots

heatmap(
    x_vals,         # x-axis values (sigma)
    y_vals,         # y-axis values (epsilon)
    z_matrix,       # z values (avg_shannon)
    xlabel = "Sigma",
    ylabel = "Epsilon",
    title = "Avg Shannon Index",
    colorbar_title = "Avg Shannon",
    aspect_ratio = :equal
)

