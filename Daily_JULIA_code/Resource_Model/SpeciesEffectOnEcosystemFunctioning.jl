# TO RUN THIS SCRIPT YOU'LL NEDD TO RUN Searching_keystone_species.jl FIRST SO
# YOU GET ALL RESULTS LIST, WICH CONTAINS THE INFO OF ALL REMOVALS

# --------------------------
# Initialize an empty DataFrame to store the results
species_ecosystem_effect = DataFrame(
    species_name = String[],
    average_effect = Float64[],
    average_effect_sd = Float64[],
    average_effect_standardized = Float64[],
    average_effect_standardized_sd = Float64[],
    clear_directional_efect = Bool[]
)

# Loop through all unique species in the dataset
all_species = unique(vcat([dd.sp_removed for dd in all_results_list]...))
filter!(x -> x != "none", all_species)

for species_name in all_species
    # Initialize storage for biomass impacts
    delta_total_biomass_values = Float64[]
    standardized_biomass_values = Float64[]

    # Loop through all results (cells) to check for this species
    for dd in all_results_list
        # Check if the species was removed in this cell (i.e., it was present initially)
        if species_name in dd.sp_removed
            # Extract the delta_total_biomass for the species in this cell
            delta_total_biomass = dd[dd.sp_removed .== species_name, :delta_total_biomass]

            # Add the biomass values to the list
            append!(delta_total_biomass_values, delta_total_biomass)
            # println(species_name)
            # If the species exists in the biomass table, standardize the value
            if species_name in birmmals_biomass_fixed.species
                species_biomass = birmmals_biomass_fixed[birmmals_biomass_fixed.species .== species_name, :biomass][1]
                append!(standardized_biomass_values, delta_total_biomass ./ species_biomass)
            else
                println("Warning: Species $species_name not found in birmmals_biomass_fixed. Skipping standardization.")
            end
        end
    end

    # Compute the mean delta biomass values if any data was found
    if !isempty(delta_total_biomass_values)
        average_effect = mean(delta_total_biomass_values)
        average_effect_sd = std(delta_total_biomass_values)
        average_effect_standardized = isempty(standardized_biomass_values) ? NaN : mean(standardized_biomass_values)
        average_effect_standardized_sd = isempty(standardized_biomass_values) ? NaN : std(standardized_biomass_values)
        if abs(average_effect_sd) < abs(average_effect)
            clear_directional_efect = true
        else
            clear_directional_efect = false
        end

        # Append results to the DataFrame
        push!(species_ecosystem_effect, (
            species_name,
            average_effect,
            average_effect_sd,
            average_effect_standardized,
            average_effect_standardized_sd,
            clear_directional_efect
        ))
    else
        println("Species $species_name was not found in any cell.")
    end
end

# Display the resulting DataFrame
println(species_ecosystem_effect)

