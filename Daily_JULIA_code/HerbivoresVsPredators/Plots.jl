######################## PLOTS ############################
# Extract data for plotting
predators = collect(keys(predator_prey_count))
total_prey = [predator_prey_count[p][1] for p in predators]
herbivores_eaten = [predator_prey_count[p][2] for p in predators]
predators_eaten = [predator_prey_count[p][3] for p in predators]

################## Stacked Bar Chart ####################
begin
    # Sort the predators by total prey eaten
    sorted_indices = sortperm(total_prey, rev=true)
    sorted_predators = predators[sorted_indices]
    sorted_herbivores_eaten = herbivores_eaten[sorted_indices]
    sorted_predators_eaten = predators_eaten[sorted_indices]

    fig1 = Figure(resolution=(800, 400))
    ax1 = Axis(fig1[1, 1], title="Prey Composition by Predator", xlabel="Predators", ylabel="Number of Prey")

    # Stacked bars
    bar_positions = 1:length(sorted_predators)
    barplot!(ax1, bar_positions, sorted_herbivores_eaten, label="Herbivores Eaten", color=:blue)
    barplot!(ax1, bar_positions, sorted_predators_eaten, offset=sorted_herbivores_eaten, label="Predators Eaten", color=:orange)

    # ax1.xticks = (bar_positions, sorted_predators)
    Legend(fig1[1, 2], ax1)
    fig1
end

################## Scatterplot ###########################
# Proportion of herbivores in the diet
proportion_herbivores = [herbivores_eaten[i] / total_prey[i] for i in 1:length(predators)]
begin
    
    fig2 = Figure(resolution=(800, 400))
    ax2 = Axis(fig2[1, 1], title="Predator Diversity and Dominance",
               xlabel="Total Prey Eaten", ylabel="Proportion of Herbivores in Diet")

    # Scatterplot
    scatter!(ax2, total_prey, proportion_herbivores, 
             markersize=0.2 .* predators_eaten, color=:green, label="Predators")

    # # Add text labels
    # for (i, predator) in enumerate(predators)
    #     text!(ax2, predator, position=(total_prey[i], proportion_herbivores[i] + 0.05), align=:center, fontsize=10)
    # end

    fig2
end

################# Heatmap for Prey Availability ########################
# Identify species that are eaten
filtered_spain_names = spain_names
# Plot the filtered heatmap
cut = true
begin 
    fig3 = Figure(resolution=(800, 400))
    ax3 = Axis(fig3[1, 1], title="Prey Availability per Predator", xlabel="Prey Species", ylabel="Predators")
    names_herb_pred = vcat(herbivore_names, predator_names)
    ibe = iberian_interact_NA[names_herb_pred, names_herb_pred]
    if cut == true
        ibe = ibe[156:256, :]
    end
    # Create heatmap with filtered matrix
    heatmap!(ax3, ibe', colormap=:viridis, interpolate=false)
    
    # Update ticks with filtered species names
    # ax3.xticks = (1:length(filtered_spain_names), filtered_spain_names)
    # ax3.yticks = (1:length(predators), predators)
        
    fig3
end

println("Corvus corone is in position: ", findfirst(x -> x == "Corvus corone", spain_names))
################### Predator Overlap Heatmap #####################
# Compute overlap matrix
overlap_matrix = [
    length(intersect(predator_prey_dict[p1], predator_prey_dict[p2])) for p1 in predators, p2 in predators
]

# Plot predator overlap heatmap
begin
    fig4 = Figure(resolution=(800, 400))
    ax4 = Axis(fig4[1, 1], title="Predator Overlap (Shared Prey Species)", xlabel="Predators", ylabel="Predators")
    heatmap!(ax4, overlap_matrix, colormap=:blues, colorrange=(0, maximum(overlap_matrix)))

    # Update ticks with predator names
    ax4.xticks = (1:length(predators), predators)
    ax4.yticks = (1:length(predators), predators)

    fig4
end

# Compute relative overlap matrix
relative_overlap_matrix = [
    length(intersect(predator_prey_dict[p1], predator_prey_dict[p2])) /
    max(length(predator_prey_dict[p1]), length(predator_prey_dict[p2])) for p1 in predators, p2 in predators
]

# Plot relative predator overlap heatmap
begin
    fig4 = Figure(resolution=(800, 400))
    ax4 = Axis(fig4[1, 1], title="Relative Predator Overlap (Shared Prey Species)", 
               xlabel="Predators", ylabel="Predators")
    
    # Plot the relative overlap heatmap
    heatmap!(ax4, relative_overlap_matrix, colormap=:blues, colorrange=(0, 1))

    # Update ticks with predator names
    # ax4.xticks = (1:length(predators), predators)
    # ax4.yticks = (1:length(predators), predators)

    fig4
end
