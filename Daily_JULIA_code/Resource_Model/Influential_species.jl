#### THIS CODE IS A COPY FROM A CHUNK OF SEARCHING_KEYSTONE_SPECIES.JL
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
using DataFrames

function create_species_count_df(df_to_work_with; only_most_influential::Bool=true)
    # Initialize the DataFrame with proper column types.
    species_count_df = DataFrame(
        species_name = String[],
        count = Int[],             # Number of cells where this species appears (for removal)
        total_ind_ext = Int[],       # Sum of indirect extinctions (across cells)
        number_of_presences = Int[], # Overall number of presences (from external info)
        cell_id = Vector{Int}[],     # List of cell IDs where this species appeared
    )

    if only_most_influential
        # For each cell, record only the entry (or entries) corresponding to the minimum survival rate.
        for dd in df_to_work_with
            # Find the minimum survival rate in this cell.
            min_survival = minimum(dd.survival_rate)
            # (The following check has been removed so that even if the removal causes 0 indirect extinctions, we record it.)
            # Get all entries with the minimum survival_rate.
            min_survival_entries = findall(x -> x == min_survival, dd.survival_rate)
            for entry in min_survival_entries
                species_name = dd[entry, :sp_removed]
                ind_ext_names = dd[entry, :ind_ext_name]  # This is assumed to be a vector (possibly empty)
                cell_id = dd[entry, :cell]
                ind_ext_count = length(ind_ext_names)  # This may be 0, but we record it.
                if species_name in species_count_df.species_name
                    indexx = findfirst(x -> x == species_name, species_count_df.species_name)
                    push!(species_count_df[indexx, :cell_id], cell_id)
                    species_count_df[indexx, :count] += 1
                    species_count_df[indexx, :total_ind_ext] += ind_ext_count
                else
                    push!(species_count_df, (
                        species_name,
                        1,                # count starts at 1
                        ind_ext_count,    # total_ind_ext is the number of indirect extinctions
                        count_number_of_presences(species_name; info=false),
                        [cell_id]         # cell_id as a single-element array
                    ))
                end
            end
        end
    else
        # For each cell, record for every species (each row where sp_removed != "none").
        for dd in df_to_work_with
            for row in eachrow(dd)
                if row.sp_removed == "none"
                    continue
                end
                species_name = row.sp_removed
                ind_ext_names = row.ind_ext_name  # A vector, possibly empty.
                cell_id = row.cell
                ind_ext_count = length(ind_ext_names)  # May be 0.
                if species_name in species_count_df.species_name
                    indexx = findfirst(x -> x == species_name, species_count_df.species_name)
                    push!(species_count_df[indexx, :cell_id], cell_id)
                    species_count_df[indexx, :count] += 1
                    species_count_df[indexx, :total_ind_ext] += ind_ext_count
                else
                    push!(species_count_df, (
                        species_name,
                        1,
                        ind_ext_count,
                        count_number_of_presences(species_name; info=false),
                        [cell_id]
                    ))
                end
            end
        end
    end

    # Compute average indirect extinctions per cell for each species.
    species_count_df[!, :avg_ind_ext] = species_count_df.total_ind_ext ./ species_count_df.number_of_presences

    # Also, add a standardized count (count divided by number of presences).
    insertcols!(species_count_df, 3,
        :stand_count => species_count_df.count ./ species_count_df.number_of_presences
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
function plot_species_metrics(
    species_count_df,
    new_all_results_list,
    selected_metric::Symbol;
    only_most_influential = true,
    count_or_stand_count = true,
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    by_name_or_by_TL = true,
    palette = custom_palette
)
    # Sort species and choose left-plot values based on the flag.
    if only_most_influential
        if count_or_stand_count
            sorted_species = sort(species_count_df, :count, rev=true)
            left_values = sorted_species.count
            left_ylabel = "Frequency"
            left_title = "Most Influential Species Not Standardised by Prevalence counting only_most_influential"
        else
            sorted_species = sort(species_count_df, :stand_count, rev=true)
            left_values = sorted_species.stand_count
            left_ylabel = "Standardized Frequency"
            left_title = "Most Influential Species Standardised by Prevalence counting only_most_influential"
        end
    else
        sorted_species = sort(species_count_df, :avg_ind_ext, rev=true)
        left_values = sorted_species.avg_ind_ext
        left_ylabel = "Average Indirect Extinctions"
        left_title = "All species, sorted by Average Species Indirect Extinctions"
    end

    # Extract species names.
    species_names = sorted_species.species_name

    # Prepare a temporary DataFrame to aggregate the selected metric from new_all_results_list.
    temp_df = DataFrame(
        species = species_names,
        temp_vect = zeros(length(species_names)),
        count = zeros(length(species_names))
    )

    # Loop over each species in the sorted list and over all cell results.
    for sp in species_names
        for dd in new_all_results_list
            # If the species is recorded in the cell (in the sp_removed column),
            # extract the value for the selected metric.
            if sp in dd.sp_removed
                new_value = dd[dd.sp_removed .== sp, selected_metric]
                temp_df[temp_df.species .== sp, :temp_vect] .+= new_value
                temp_df[temp_df.species .== sp, :count] .+= 1.0
            end
        end
    end

    # Calculate the mean value for the selected metric.
    mean_values = temp_df.temp_vect ./ temp_df.count

    # Prepare a new DataFrame that will contain the metric for each species.
    new_df = DataFrame(species_names = species_names)
    new_df[!, selected_metric] = mean_values

    # Color structure
    # Assign colors: if by_name_or_by_TL, use discrete colors; otherwise, use trophic level (continuous).
    if by_name_or_by_TL
        colors = [species_names[i] in herbivore_names ? :blue :
                  (species_names[i] in predator_names ? :red : :gray)
                  for i in 1:length(species_names)]
    else
        # Use TrophInd: assume it has columns :Species and :TL.
        # For each species in the sorted list, get its trophic level.
        troph_levels = [TrophInd[TrophInd.Species .== species_names[i], :TL][1] for i in 1:length(species_names)]
        tl_min = minimum(troph_levels)
        tl_max = maximum(troph_levels)
        norm_tl = [(tl - tl_min) / (tl_max - tl_min) for tl in troph_levels]
        # Get a continuous colormap, e.g., viridis.
        colormap = cgrad(palette)
        # Instead of calling colormap(x), we index into it:
        colors = [colormap[x] for x in norm_tl]
    end

    # Define the figure.
    fig = Figure(resolution = (1200, 600))

    # --- Left Plot: Barplot using left_values ---
    ax1 = Axis(fig[1, 1],
        title = left_title,
        xlabel = "Species",
        ylabel = left_ylabel,
        xticks = (1:length(species_names), species_names),
        xticklabelrotation = π/4, 
        xticklabelsize = 8
    )
    barplot!(ax1, 1:length(species_names), left_values, color = colors)

    # --- Right Plot: Scatter Plot of the selected metric ---
    ax2 = Axis(fig[1, 2],
        title = "Average $(selected_metric) per Species",
        xlabel = "Species",
        ylabel = "Value",
        xticks = (1:length(species_names), species_names),
        xticklabelrotation = π/4, 
        xticklabelsize = 8
    )
    scatter!(ax2, 1:length(species_names), mean_values, color = colors)
    by_name_or_by_TL ? identity : Colorbar(fig[1, 3], limits = (tl_min, tl_max), colormap = cgrad(palette))
    display(fig)
    return fig
end

# Example usage:
selected_metric = :clustering # Change to :indegree, :outdegree, :total_degree, :closeness, or :clustering
plot_species_metrics(species_count_df, new_all_results_list, selected_metric)
