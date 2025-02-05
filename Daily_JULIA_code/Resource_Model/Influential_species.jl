#### THIS CODE IS A COPY FROM A CHUNK OF SEARCHING_KEYSTON_SPECIES.JL
#### IT CAN BE USED AFTER RUNNIN Applying_results.jl and IT WILL GIVE 
#### INFORMATION ON THE MOST INFLUENTIAL SPECIES
df_to_work_with = all_results_list_even_pi
# SEEMS ALL_RESULTS_LIST IS CORRECT
comparison_df = DataFrame(big_p_sr = [], cell_from_big_p = [], cell_from_all = [], all_results_list_sr = [], diff = [])
for bp_row in eachrow(Big_P_results_maximised)
    cell_id = bp_row[:cell_id]
    all_results_row = findfirst(x -> x[1, :cell] == cell_id, df_to_work_with)
    
    if all_results_row !== nothing
        best_params_sr = bp_row[:survival_rate]
        all_results_list_sr = df_to_work_with[all_results_row][1, :survival_rate]
        diff = best_params_sr - all_results_list_sr
        push!(comparison_df, (best_params_sr, cell_id, df_to_work_with[all_results_row][1, :cell], all_results_list_sr, diff))
    end
end

#################################################################
#################################################################
function create_species_count_df(df_to_work_with)
    # FIND THE MOST INFLUENTIAL SPECIES
    # Initialize the DataFrame with proper column types
    species_count_df = DataFrame(
        species_name = String[],
        count = Int[],
        ind_ext_names_length = Int[],
        number_of_presences = Int[],
        cell_id = Vector{Int}[],  # Assuming cell IDs are integers
        ind_ext_names = Vector{String}[],
    )

    # Loop through all results
    for dd in df_to_work_with
        # Find the entry with the lowest survival_rate
        min_survival = minimum(dd.survival_rate)
        min_survived = minimum(dd.survived_herbs + dd.survived_preds)

        # Skip if the minimum survival condition matches the first row's survivors minus 1
        if min_survived == dd[1, :survived_herbs] + dd[1, :survived_preds] - 1
            continue
        end

        # Get all entries with the minimum survival_rate
        min_survival_entries = findall(x -> x == min_survival, dd.survival_rate)
    
        for min_survival_entry in min_survival_entries
            # Get the species name, individual extinction names, and cell ID
            species_name = dd[min_survival_entry, :sp_removed]
            ind_ext_names = dd[min_survival_entry, :ind_ext_name]
            cell_id = dd[min_survival_entry, :cell]

            # Check if species_name is already in the DataFrame
            if species_name in species_count_df.species_name
                # Find the index of the existing entry
                indexx = findfirst(x -> x == species_name, species_count_df.species_name)

                # Update the ind_ext_names and cell_id if they are not already present
                for ind_ext_name in ind_ext_names
                    if !(ind_ext_name in species_count_df[indexx, :ind_ext_names])
                        push!(species_count_df[indexx, :ind_ext_names], ind_ext_name)
                    end
                end

                if !(cell_id in species_count_df[indexx, :cell_id])
                    push!(species_count_df[indexx, :cell_id], cell_id)
                end

                # Increment the count
                species_count_df[indexx, :count] += 1
            else
                # Add a new entry for this species
                push!(species_count_df, (
                    species_name,                  # species_name
                    1,                             # count starts at 1
                    0,
                    count_number_of_presences(species_name; info = false),
                    [cell_id],                     # cell_id as a single-element array
                    ind_ext_names,                  # ind_ext_names as-is
                ))
            end
        end
    end

    for i in birmmals_names
        if !(i in species_count_df.species_name)
            push!(species_count_df, (
                i,
                0,              # :count
                0,              # :ind_ext_names_length
                count_number_of_presences(i; info = false),
                Int[],          # Empty Vector{Int} for :cell_id
                String[]        # Vector{String} for :ind_ext_names
            ))
        end
    end

    species_count_df[!, :ind_ext_names_length] = length.(species_count_df.ind_ext_names)
    insertcols!(species_count_df, 3,
        :stand_count => species_count_df[!, :count] ./ species_count_df[!, :number_of_presences]
    )

    return species_count_df
end
# species_count_df = create_species_count_df(all_results_list_even_pi)

#### PLOT MOST INFLUENTIAL SPECIES AND PI BIOMASS VALUES
begin
    ending = 205
    species_count_df = sort(species_count_df, :count, rev=true)
    fig = Figure(resolution = (1000, 800))
    ax = Axis(fig[1, 1], title="Frequency of most influencial species", xlabel="Species", ylabel="Frequency")
    MK.barplot!(ax, 1:length(species_count_df.species_name[1:ending]), species_count_df.count[1:ending])
    ax.xticks = (1:length(species_count_df.species_name[1:ending]), species_count_df.species_name[1:ending])
    ax.xticklabelrotation = π/2.5
    ax.xticklabelsize = 6

    # # species_count_df = sort(species_count_df, :stand_count, rev=true)
    # ax2 = Axis(fig[1, 2], title="Standardized frequency by prevalence", xlabel="Species", ylabel="Frequency")
    # MK.barplot!(ax2, 1:length(species_count_df.species_name), species_count_df.stand_count)
    # ax2.xticks = (1:length(species_count_df.species_name), species_count_df.species_name)
    # ax2.xticklabelrotation = π/2.5
    # ax2.xticklabelsize = 6
        
    sorted_birmmals_biomass_fixed = sort(birmmals_biomass_fixed, :biomass, rev=true)
    ax2 = Axis(fig[1, 2], title="H0_values Biomass", xlabel="Species", ylabel="Biomass")
    MK.barplot!(ax2, 1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.biomass)
    ax2.xticks = (1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.species)
    ax2.xticklabelrotation = π/2.5
    ax2.xticklabelsize = 6
    
    display(fig)

    # println(species_count_df.species_name)
end

begin
    ending = 205
    species_count_df = sort(species_count_df, :stand_count, rev=true)
    fig = Figure(resolution = (1000, 800))
    ax2 = Axis(fig[1, 1], title="Standardized frequency by prevalence", xlabel="Species", ylabel="Frequency")
    MK.barplot!(ax2, 1:length(species_count_df.species_name[1:ending]), species_count_df.stand_count[1:ending])
    ax2.xticks = (1:length(species_count_df.species_name[1:ending]), species_count_df.species_name[1:ending])
    ax2.xticklabelrotation = π/2.5
    ax2.xticklabelsize = 8

    display(fig)
end

#### PLOT PI BIOMASS VALUES (LOG SCALE or not)
begin
    dis = true
    save_the_plot = true
    log = true
    fig = Figure(resolution = (1000, 800))
    sorted_birmmals_biomass_fixed = sort(birmmals_biomass_fixed, :biomass, rev=true)
    ax2 = Axis(
        fig[1, 1], title="H0_values Biomass", 
        xlabel="Species", ylabel="Biomass",
        yscale=log ? log10 : identity
        )
    MK.barplot!(ax2, 1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.biomass)
    ax2.xticks = (1:nrow(birmmals_biomass_fixed), sorted_birmmals_biomass_fixed.species)
    ax2.xticklabelrotation = π/2.5
    ax2.xticklabelsize = 6
    
    if dis
        display(fig)
    end
    word = log ? "log" : "absolute"
    if save_the_plot
        save("Plots/H0_values_biomass_$word.png", fig) 
    end
end

#### PLOT CORRELATION BETWEEN PI BIOMASS VALUES AND MOST INFLUENTIAL SPECIES
begin
    dis = true
    log = false
    save_the_plot = true
    fig = Figure(resolution = (1000, 800))
    ax = Axis(
        fig[1, 1],
         title="Correlation between H0_values and frequency of most influencial species", 
         xlabel="Frequency", ylabel="H0_values",
         yscale=log ? log10 : identity,
        # xscale=log ? log10 : identity
    )
    matched_df = innerjoin(species_count_df, sorted_birmmals_biomass_fixed, on=:species_name => :species)
    MK.scatter!(ax, matched_df.count, matched_df.biomass)
    ax.xlabel = "Frequency of most influencial species"
    ax.ylabel = "Biomass of species"
    
    if dis
        display(fig)
    end
    word = log ? "logscale" : ""
    if save_the_plot
        save("Plots/H0_values_VS_most_influencial_$word.png", fig) 
    end
end

###### PLOT MOST INFLUENTIAL SPECIES, AND SOME OF THEIR NETWORK METRICS
function plot_species_metrics(species_count_df, new_all_results_list, selected_metric::Symbol)
    # Sort species by frequency of being most influential
    sorted_species = sort(species_count_df, :count, rev=true)
    
    # Extract species names & frequencies
    species_names = sorted_species.species_name
    species_frequencies = count_number_of_presences.(species_names; info = false)
    species_counts = sorted_species.count

    new_df = DataFrame(species_names = species_names, species_frequencies = species_frequencies)
    temp_df = DataFrame(
        species = species_names,
        temp_vect = zeros(length(species_names)),
        count = zeros(length(species_names))
    )

    for i in species_names
        for j in 1:length(new_all_results_list)
            # println(new_all_results_list[j].sp_removed)
            if i in new_all_results_list[j].sp_removed
                new_value = new_all_results_list[j][new_all_results_list[j].sp_removed .== i, selected_metric]
                # println("new_value: ", new_value)
                temp_df[temp_df.species .== i, :temp_vect] += new_value
                # println("hello", temp_df[temp_df.species .== i, :count])
                temp_df[temp_df.species .== i, :count] += Any[1.0]
            end
        end
    end

    mean_values = temp_df.temp_vect ./ temp_df.count
    new_df[!, selected_metric] = mean_values

    # Define figure
    fig = Figure(resolution = (1200, 600))

    # --- Frequency Plot (Left) ---
    ax1 = Axis(fig[1, 1], 
        title="Most Influential Species", 
        xlabel="Species", ylabel="Frequency",
        xticks=(1:length(species_names), species_names),
        xticklabelrotation=π/4, xticklabelsize=8)

    barplot!(ax1, 1:length(species_names), species_counts, color=:blue)

    # --- Species-Specific Metric Plot (Right) ---
    ax2 = Axis(fig[1, 2], 
        title="Average $selected_metric per Species", 
        xlabel="Species", ylabel="Value",
        xticks=(1:length(species_names), species_names),
        xticklabelrotation=π/4, xticklabelsize=8)
        # Regression line IT DOES NOT WORK FOR NOW
        # regression_data = DataFrame(x=1:length(species_names), y=mean_values)
        # model = lm(@formula(y ~ x), regression_data)
        # line_x = range(1, length(species_names), length=100)
        # line_y = coef(model)[1] .+ coef(model)[2] .* line_x
        # lines!(ax2, line_x, line_y, linewidth=2, color=:red)

    scatter!(ax2, 1:length(species_names), mean_values, color=:red)
    # errorbars!(ax2, 1:length(species_names), metric_values, [0.0 for _ in metric_values], color=:black)  # Adjust error bars if needed

    display(fig)
end

# Example usage:
selected_metric = :clustering # Change to :indegree, :outdegree, :total_degree, :closeness, or :clustering
plot_species_metrics(species_count_df, new_all_results_list, selected_metric)
