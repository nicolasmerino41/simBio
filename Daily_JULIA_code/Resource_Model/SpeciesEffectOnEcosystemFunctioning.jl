# TO RUN THIS SCRIPT YOU'LL NEDD TO RUN applying_results.jl FIRST SO
# YOU GET ALL RESULTS LIST, WICH CONTAINS THE INFO OF ALL REMOVALS

# --------------------------
# The SEEF function is defined in Functions/SEEF_function.jl
see_not_even = SEEF(all_results_list)
see_even     = SEEF(all_results_list_even_pi)

see_to_plot  = see_even
########## PLOTTING THE DATA ############
begin
    log = false
    avg_eff_or_avg_eff_stand = true  # true for average_effect, false for average_effect_standardized
    # 1) Sort descending by average_effect
    see = deepcopy(see_to_plot)
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
        yscale = log ? log10 : identity,
        xticklabelsize = 6 
    )

    # 4) Plot error bars
    MK.errorbars!(ax, 1:length(species), avgs, sds, color=:black)

    # 5) Plot points
    MK.scatter!(ax, 1:length(species), avgs, color=:blue, markersize=10)

    # 6) Save and display
    # save("species_effect_plot.png", fig)
    display(fig)
end

##### PLOTING AVERAGE EFFECT OF EACH SPECIES VS THEIR METRICS ########
begin
    avg_in_y = true
    # First, join the two DataFrames on species name.
    # species_ecosystem_effect has column :species_names and avg_species_metrics has column :species.
    joined_df = innerjoin(
        see_to_plot, avg_species_metrics,
        on = :species_name
    )
    println("minimum avg_effect is", minimum(joined_df.average_effect))
    # Define the list of metrics you want to plot from avg_species_metrics.
    metrics = [:mean_indegree, :mean_outdegree, :mean_total_degree, :mean_betweenness, :mean_closeness, :mean_clustering]
    labels  = ["Mean Indegree", "Mean Outdegree", "Mean Total Degree", "Mean Betweenness", "Mean Closeness", "Mean Clustering"]

    # Create a figure with multiple subplots (e.g., 2 rows × 3 columns).
    fig = Figure(resolution = (1200, 800))
    n = length(metrics)
    for (i, metric) in enumerate(metrics)
        row = div(i-1, 3) + 1
        col = mod(i-1, 3) + 1
        ax = Axis(fig[row, col],
            title = avg_in_y ? labels[i] * " vs Average Effect" : "Average Effect vs " * labels[i],
            xlabel = avg_in_y ? labels[i] : "Average Effect", 
            ylabel = invert_y_and_x ? "Average Effect" : labels[i],
        )
        
        scatter!(
            ax,
            invert_y_and_x ? joined_df[!, metric] : joined_df[!, :average_effect],
            invert_y_and_x ? joined_df[!, :average_effect] : joined_df[!, metric],
            markersize = 8,
            color = :blue
        )
    end

    display(fig)
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

# COMPARE AVERAGE EFFECT VS H0_VALUES OR H_END_VALUES
