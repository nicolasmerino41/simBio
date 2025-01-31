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
            # println("delta_total_biomass", delta_total_biomass)
            if delta_total_biomass[1] < 50.0 && delta_total_biomass[1] > -50.0
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

##### PLOTTING THE DATA ############
begin
    log = false
    avg_eff_or_avg_eff_stand = false  # true for average_effect, false for average_effect_standardized
    # 1) Sort descending by average_effect
    see = deepcopy(species_ecosystem_effect)
    sort!(see, avg_eff_or_avg_eff_stand ? :average_effect : :average_effect_standardized, rev=true)

    # 2) Extract sorted columns
    species = see.species_name
    avgs    = avg_eff_or_avg_eff_stand ? see.average_effect : see.average_effect_standardized
    sds     = avg_eff_or_avg_eff_stand ? see.average_effect_sd : see.average_effect_standardized_sd

    # 3) Define figure and axis
    fig = Figure(resolution=(1100, 700))
    ax = Axis(fig[1, 1],
        title = "Species Effects (± SD)",
        xlabel = "Species",
        ylabel = avg_eff_or_avg_eff_stand ? "Average Effect" : "Average Effect (Standardized)",
        xticks = (1:length(species), species),  # Label x-axis with species names
        xticklabelrotation = π/4,               # Rotate x-axis labels
        xticklabelalign = (:right, :center),     # Align labels to prevent overlap
        yscale = log ? log10 : identity
    )

    # 4) Plot error bars
    MK.errorbars!(ax, 1:length(species), avgs, sds, color=:black)

    # 5) Plot points
    MK.scatter!(ax, 1:length(species), avgs, color=:blue, markersize=10)

    # 6) Save and display
    # save("species_effect_plot.png", fig)
    display(fig)
end

AAAA = all_results_list[20]