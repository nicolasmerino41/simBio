using Glob
using DataFrames
using JLD2
using Statistics
using LinearAlgebra
using CSV

# Define the folder containing the .jls files
folder_path = "D:/Results from Direct Sampling with herbivory drago 14-10-2024"

# List all .jls files in the folder
jls_files = sort(glob("$folder_path/*.jls"))

# Initialize an empty DataFrame to store the results
results_df = DataFrame(
    sigma = Float64[],
    epsilon = Float64[],
    alfa = Float64[],
    sigma_comp = Float64[],
    assymetry = Float64[],
    slope = Float64[]
)

# Define the bins for trophic levels
bin_size = 0.2
# Assuming TrophInd_vector is defined and contains trophic levels for each species
# TrophInd_vector = ... (Provide your vector here)

# Iterate over all .jls files
for file in jls_files
    # Load the data
    array_output = deserialize(file)
    
    # Check for NaNs in the data (adjust based on your data structure)
    if any(isnan, array_output)
        println("Skipping file with NaNs: ", file)
        continue
    end

    # Extract parameters from the filename
    # Filename format: "3.0-5.0-0.1-1.0-0.33.jls"
    filename = basename(file)
    filename_without_ext = splitext(filename)[1]
    params_str = split(filename_without_ext, "-")
    params = parse.(Float64, params_str)
    sigma, epsilon, alfa, sigma_comp, assymetry = params

    # Combine abundances for the entire region
    combined_abundances = array_output # Adjust if necessary
    biomass = zeros(256)  # Assuming 256 species
    for i in idx
        biomass += combined_abundances[i].a
    end

    # Define trophic levels if not already defined
    # TrophInd_vector = ... (Provide your vector here)
    
    # Bin the biomass by trophic level
    min_trophic_level = floor(minimum(TrophInd_vector) / bin_size) * bin_size
    max_trophic_level = ceil(maximum(TrophInd_vector) / bin_size) * bin_size
    bins = min_trophic_level:bin_size:max_trophic_level

    biomass_per_bin = zeros(length(bins) - 1)
    for i in 1:(length(bins) - 1)
        bin_indices = (TrophInd_vector .>= bins[i]) .& (TrophInd_vector .< bins[i + 1])
        biomass_per_bin[i] = sum(biomass[bin_indices])
    end

    # Calculate the steepness of the pyramid
    # Approach: Linear regression on log-transformed biomass vs. trophic level
    # Handle zeros to avoid log(0)
    nonzero_indices = biomass_per_bin .> 0
    if sum(nonzero_indices) < 2
        println("Not enough data points for regression in file: ", file)
        continue
    end
    log_biomass = log10.(biomass_per_bin[nonzero_indices])
    trophic_levels = (bins[1:end-1] .+ bin_size / 2)[nonzero_indices]  # Bin midpoints

    # Fit linear model
    X = [trophic_levels ones(length(trophic_levels))]
    coeffs = X \ log_biomass
    slope = coeffs[1]

    # Append the results to the DataFrame
    push!(results_df, (
        sigma = sigma,
        epsilon = epsilon,
        alfa = alfa,
        sigma_comp = sigma_comp,
        assymetry = assymetry,
        slope = slope
    ))
end

# Save the results to a CSV file
CSV.write("pyramid_slope_from_DirectSampling.csv", results_df)
