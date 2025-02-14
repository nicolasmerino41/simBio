# TO RUN THIS SCRIPT YOU'LL NEDD TO RUN applying_results.jl FIRST SO
# YOU GET ALL RESULTS LIST, WICH CONTAINS THE INFO OF ALL REMOVALS

# --------------------------
# The SEEF function is defined in Functions/SEEF_function.jl
for i in 1:length(all_results_list)
    if length(names(all_results_list[i])) > 2
       
    else
        println("DataFrame $i does not have column :spo_removed")
    end
end
for i in 1:length(all_results_list)
    names(all_results_list[i])[2]
end
see_not_even = SEEF(all_results_list)
see_even     = SEEF(all_results_list_even_pi)

see_to_plot  = see_not_even
########## PLOTTING THE DATA ############
function plot_species_effects(
    see_to_plot;
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    log = false,
    standardise_by_H0 = false,
    resolution = (1100, 500),
    by_name_or_by_TL = true,  # if true, assign colors by name; if false, use continuous TL from TrophInd
    palette = custom_palette
)
    # 1) Sort descending by average_effect
    see = deepcopy(see_to_plot)
    sort!(see, standardise_by_H0 ? :average_effect_standardized : :average_effect, rev=true)

    # 2) Extract sorted columns
    species = see.species_name
    avgs    = standardise_by_H0 ? see.average_effect_standardized : see.average_effect
    sds     = standardise_by_H0 ? see.average_effect_standardized_sd : see.average_effect_sd

    # 3) Assign colors: if by_name_or_by_TL, use discrete colors; otherwise, use trophic level (continuous).
    if by_name_or_by_TL
        colors = [species[i] in herbivore_names ? :blue :
                  (species[i] in predator_names ? :red : :gray)
                  for i in 1:length(species)]
    else
        # Use TrophInd: assume it has columns :Species and :TL.
        # For each species in the sorted list, get its trophic level.
        troph_levels = [TrophInd[TrophInd.Species .== species[i], :TL][1] for i in 1:length(species)]
        tl_min = minimum(troph_levels)
        tl_max = maximum(troph_levels)
        norm_tl = [(tl - tl_min) / (tl_max - tl_min) for tl in troph_levels]
        # Get a continuous colormap, e.g., viridis.
        colormap = cgrad(palette)
        # Instead of calling colormap(x), we index into it:
        colors = [colormap[x] for x in norm_tl]
    end

    # 4) Define figure and axis.
    fig = Figure(resolution = resolution)
    ax = Axis(fig[1, 1],
        title = "Species Effects (± SD)",
        xlabel = "Species",
        ylabel = standardise_by_H0 ? "Average Effect (Standardized)" : "Average Effect",
        xticks = (1:length(species), species),  # Label x-axis with species names.
        xticklabelrotation = π/4,
        xticklabelalign = (:right, :center),
        yscale = log ? log10 : identity,
        xticklabelsize = 6
    )

    # 5) Plot error bars and points.
    MK.errorbars!(ax, 1:length(species), avgs, sds, color = colors)
    MK.scatter!(ax, 1:length(species), avgs, color = colors, markersize = 10)

    # 6) Plot a reference line at zero.
    MK.lines!(ax, 1:length(species), fill(0.0, length(species)), color = :red)

    by_name_or_by_TL ? identity : Colorbar(fig[1, 2], limits = (tl_min, tl_max), colormap = cgrad(palette))

    display(fig)
    return fig
end

##### PLOTING AVERAGE EFFECT OF EACH SPECIES VS THEIR METRICS ########
function plot_average_effect_vs_metrics(
    see_to_plot = see_even;
    avg_species_metrics = avg_species_metrics,
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    by_name_or_by_TL = true,  # if true, color by species type; if false, color continuously by TL from TrophInd
    palette = custom_palette
)
    avg_in_y = true

    # Join the two DataFrames on species name.
    joined_df = innerjoin(see_to_plot, avg_species_metrics, on = :species_name)
    println("Minimum average effect is ", minimum(joined_df.average_effect))

    # Define the list of metrics and labels to plot.
    metrics = [:mean_indegree, :mean_outdegree, :mean_total_degree, :mean_betweenness, :mean_closeness, :mean_clustering]
    labels  = ["Mean Indegree", "Mean Outdegree", "Mean Total Degree", "Mean Betweenness", "Mean Closeness", "Mean Clustering"]

    # 3) Assign colors based on species type or trophic level.
    if by_name_or_by_TL
        colors = [joined_df.species_name[i] in herbivore_names ? :blue :
                  (joined_df.species_name[i] in predator_names ? :red : :gray)
                  for i in 1:nrow(joined_df)]
    else
        # Use TrophInd: assume it has columns :Species and :TL.
        # For each species in the sorted list, get its trophic level.
        troph_levels = [TrophInd[TrophInd.Species .== species[i], :TL][1] for i in 1:length(species)]
        tl_min = minimum(troph_levels)
        tl_max = maximum(troph_levels)
        norm_tl = [(tl - tl_min) / (tl_max - tl_min) for tl in troph_levels]
        # Get a continuous colormap, e.g., viridis.
        colormap = cgrad(palette)
        # Instead of calling colormap(x), we index into it:
        colors = [colormap[x] for x in norm_tl]
    end

    # Create a figure with subplots (e.g., 2 rows × 3 columns).
    fig = Figure(resolution = (1200, 800))
    n = length(metrics)
    for (i, metric) in enumerate(metrics)
        row = div(i-1, 3) + 1
        col = mod(i-1, 3) + 1
        ax = Axis(fig[row, col],
            title = avg_in_y ? labels[i] * " vs Average Effect" : "Average Effect vs " * labels[i],
            xlabel = avg_in_y ? labels[i] : "Average Effect", 
            ylabel = avg_in_y ? "Average Effect" : labels[i]
        )
        scatter!(ax,
            avg_in_y ? joined_df[!, metric] : joined_df[!, :average_effect],
            avg_in_y ? joined_df[!, :average_effect] : joined_df[!, metric],
            markersize = 8,
            color = colors
        )
    end
    display(fig)
    return fig
end

######## PLOT SPECIES EFFECT VS THEIR METRIC VALUES IN EACH CELL ########
function plot_species_effect_vs_cell_metrics(species_identifier::Union{String, Int})
    # Prepare empty vectors for each metric and for the delta_total_biomass.
    indegree_vals    = Float64[]
    outdegree_vals   = Float64[]
    total_degree_vals = Float64[]
    betweenness_vals = Float64[]
    closeness_vals   = Float64[]
    clustering_vals  = Float64[]
    delta_vals       = Float64[]
    
    if species_identifier isa Int
        species_identifier = birmmals_names[species_identifier]
    elseif species_identifier isa String
    else
        error("species_identifier must be a String or an Int.")
    end
    # Loop over each DataFrame (each cell) in new_all_results_list.
    for dd in new_all_results_list
        
        # Determine indices where this species is recorded.        
        idxs = findall(x -> x == species_identifier, dd.sp_removed)
        
        # For each matching row, extract the values.
        for k in idxs
            push!(delta_vals, dd[k, :delta_total_biomass])
            push!(indegree_vals, dd[k, :indegree])
            push!(outdegree_vals, dd[k, :outdegree])
            push!(total_degree_vals, dd[k, :total_degree])
            push!(betweenness_vals, dd[k, :betweenness])
            push!(closeness_vals, dd[k, :closeness])
            push!(clustering_vals, dd[k, :clustering])
        end
    end

    # Prepare a vector of metric names and a corresponding vector of x-data.
    metrics = ["Indegree", "Outdegree", "Total Degree", "Betweenness", "Closeness", "Clustering"]
    data_x = [indegree_vals, outdegree_vals, total_degree_vals, betweenness_vals, closeness_vals, clustering_vals]

    # Create a figure with 2 rows × 3 columns.
    n_metrics = length(metrics)
    ncols = 3
    nrows = ceil(Int, n_metrics / ncols)
    fig = Figure(resolution = (1200, 800))

    for i in 1:n_metrics
        row = div(i-1, ncols) + 1
        col = mod(i-1, ncols) + 1
        ax = Axis(fig[row, col],
            title = "$(metrics[i]) vs Δ Total Biomass",
            xlabel = metrics[i],
            ylabel = "Δ Total Biomass"
        )
        scatter!(ax, data_x[i], delta_vals, markersize = 8, color = species_identifier in herbivore_names ? :blue : :red)
    end

    display(fig)
    return fig
end

####### COMPARE THE AVERAGE EFFECT WHEN EVEN PI AND NOT EVEN PI #######
begin
    # Assume SEEF returns a DataFrame with columns:
    #   species_name, average_effect, ...
    # Make copies and rename the average_effect column accordingly.
    even_df = deepcopy(see_even[:, [:species_name, :average_effect]])
    noteven_df = deepcopy(see_not_even[:, [:species_name, :average_effect]])
    rename!(even_df, :average_effect => :average_effect_even)
    rename!(noteven_df, :average_effect => :average_effect_noteven)

    # Join the two DataFrames on species_name.
    joined_df = innerjoin(even_df, noteven_df, on = :species_name)

    # Create the scatter plot:
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        title = "Species Average Effect: Even vs. Not Even PI",
        xlabel = "Average Effect (pi even)",
        ylabel = "Average Effect (pi not even)"
    )

    scatter!(ax, joined_df.average_effect_even, joined_df.average_effect_noteven,
        color = :blue, markersize = 10)

    # Optionally, compute and display the Pearson correlation coefficient.
    corr_val = cor(joined_df.average_effect_even, joined_df.average_effect_noteven)
    text!(ax, "r = $(round(corr_val, digits=2))", position = (0.05, 0.95),
        align = (:left, :top), color = :black, fontsize = 12)

    display(fig)
end

####### COMPARE THE SPECIES RANK IN TERMS OF AVERAGE EFFECT WHEN EVEN PI AND NOT EVEN PI #######
begin
    
    # Assume these two functions have been run and return DataFrames:
    # see_even     = SEEF(all_results_list_even_pi)
    # see_not_even = SEEF(all_results_list)
    #
    # Both DataFrames are expected to have at least these columns:
    #   species_name and average_effect

    # Make deep copies to avoid modifying the originals.
    even_df = deepcopy(see_even)
    noteven_df = deepcopy(see_not_even)

    # Join the two DataFrames on species_name.
    # (Only species that appear in both conditions will be compared.)
    joined_df = innerjoin(even_df, noteven_df, on = :species_name, makeunique=true)
    # After the join, assume the columns are named:
    #   species_name, average_effect, average_effect_1, ...
    # Rename them for clarity.
    rename!(joined_df, Dict("average_effect" => :average_effect_even, "average_effect_1" => :average_effect_noteven))

    n = nrow(joined_df)
    # Preallocate vectors to hold rank values.
    ranks_even = Vector{Int}(undef, n)
    ranks_noteven = Vector{Int}(undef, n)

    # Compute ranks: higher average effect gets rank 1.
    sorted_idx_even = sortperm(joined_df.average_effect_even, rev=true)
    sorted_idx_noteven = sortperm(joined_df.average_effect_noteven, rev=true)
    for i in 1:n
        ranks_even[sorted_idx_even[i]] = i
        ranks_noteven[sorted_idx_noteven[i]] = i
    end

    joined_df[!, :rank_even] = ranks_even
    joined_df[!, :rank_noteven] = ranks_noteven
    # Optionally, compute the difference in rank.
    joined_df[!, :rank_diff] = joined_df.rank_even .- joined_df.rank_noteven

    # Now, create a slope graph (bump chart) to show rank changes.
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        xlabel = "Scenario",
        ylabel = "Rank (1 = highest effect)",
        xticks = ([1, 2], ["Even", "Not Even"]),
        yreversed = true  # so that rank 1 is at the top
    )

    # Plot one line per species.
    for i in 1:n
        lines!(ax, [1, 2], [joined_df.rank_even[i], joined_df.rank_noteven[i]], color=:blue, linewidth=2)
        # Optionally, label species on the left side.
        text!(ax, joined_df.species_name[i], position=(1, joined_df.rank_even[i]),
              align=(:left, :center), color=:black, fontsize=10)

    end

    display(fig)
end

#### HERE WE JUST PRINT THE TOP SPECIES FOR EACH SCENARIO ####
begin
    sorted_even = sort(see_even, :average_effect, rev=true)
    sorted_noteven = sort(see_not_even, :average_effect, rev=true)

    println("Top species for even PI:")
    for i in 1:10
        println(sorted_even.species_name[i], " with average effect: ", sorted_even.average_effect[i])
    end

    println("\nTop species for not even PI:")
    for i in 1:10
        println(sorted_noteven.species_name[i], " with average effect: ", sorted_noteven.average_effect[i])
    end
end