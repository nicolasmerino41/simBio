function SEEF(df_to_work_with)
    # Initialize an empty DataFrame to store the results.
    species_ecosystem_effect = DataFrame(
        species_name = String[],
        average_effect = Float64[],
        average_effect_sd = Float64[],
        average_effect_standardized = Float64[],
        average_effect_standardized_sd = Float64[],
        clear_directional_efect = Bool[]
    )

    # Get all unique species names across all cells in df_to_work_with.
    # Here, each element of df_to_work_with is expected to have a field/column named "sp_removed".
    all_species = unique(vcat([dd.sp_removed for dd in df_to_work_with if hasproperty(dd, :sp_removed)]...))
    filter!(x -> x != "none", all_species)

    # Loop over each species.
    for species_name in all_species
        # Initialize storage for biomass impact values.
        delta_total_biomass_values = Float64[]
        standardized_biomass_values = Float64[]

        # Loop through each cell's DataFrame in df_to_work_with.
        for dd in df_to_work_with
            if species_name in dd.sp_removed
                # Extract the delta_total_biomass values for rows where the species was removed.
                delta_total_biomass = dd[dd.sp_removed .== species_name, :delta_total_biomass]
                # Only add values that are within a reasonable range.
                if delta_total_biomass[1] < 50.0 && delta_total_biomass[1] > -50.0
                    append!(delta_total_biomass_values, delta_total_biomass)
                    # Standardize the value if the species exists in birmmals_biomass_fixed.
                    if species_name in birmmals_biomass_fixed.species
                        species_biomass = birmmals_biomass_fixed[birmmals_biomass_fixed.species .== species_name, :biomass][1]
                        append!(standardized_biomass_values, delta_total_biomass ./ species_biomass)
                    else
                        println("Warning: Species $species_name not found in birmmals_biomass_fixed. Skipping standardization.")
                    end
                end
            end
        end

        # If any biomass values were collected for this species, compute summary metrics.
        if !isempty(delta_total_biomass_values)
            average_effect = mean(delta_total_biomass_values)
            average_effect_sd = std(delta_total_biomass_values)
            average_effect_standardized = isempty(standardized_biomass_values) ? NaN : mean(standardized_biomass_values)
            average_effect_standardized_sd = isempty(standardized_biomass_values) ? NaN : std(standardized_biomass_values)
            clear_directional_efect = abs(average_effect_sd) < abs(average_effect)
            
            # Append the computed values as a new row in the results DataFrame.
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

    return species_ecosystem_effect
end
