function species_distribution(species::String; resolution = (800, 600))
    the_map = deepcopy(float.(DA_sum))
    the_map .= 0.0
    if isone(species in keys(species_dict_herbivores_in_birmmals))
        herbivore = true
    elseif isone(species in keys(species_dict_predators_in_birmmals))
        herbivore = false
    else
        error("Species $species not found in species_dict_herbivores_in_birmmals or species_dict_predators_in_birmmals")
    end

    for i in idx
        if herbivore 
            if !iszero(DA_birmmals_with_pi[i].a[species_dict_herbivores_in_birmmals[species]])
                the_map[i] = 1.0
            else
                the_map[i] = 0.5
            end
        elseif !herbivore
            if !iszero(DA_birmmals_with_pi[i].a[species_dict_predators_in_birmmals[species]])
                the_map[i] = 1.0
            else
                the_map[i] = 0.5
            end
        end
    end

    fig = Figure()
    ax = Axis(fig[1, 1], title = "Distribution of $species", xlabel = "x", ylabel = "y")
    heatmap!(ax, the_map, colormap = :inferno, colorrange = (0, 1))
    Colorbar(fig[1, 2], colormap = :inferno, vertical = true)
    ax.yreversed = true
    display(fig)
end

# Example
# species_distribution("Mustela putorius")

function count_number_of_presences(species::String; info = true)
    if isone(species in keys(species_dict_herbivores_in_birmmals))
        herbivore = true
    elseif isone(species in keys(species_dict_predators_in_birmmals))
        herbivore = false
    else
        error("Species $species not found in species_dict_herbivores_in_birmmals or species_dict_predators_in_birmmals")
    end
    countt = 0
    for i in idx
        if herbivore 
            if !iszero(DA_birmmals_with_pi[i].a[species_dict_herbivores_in_birmmals[species]])
                countt += 1
            end
        elseif !herbivore
            if !iszero(DA_birmmals_with_pi[i].a[species_dict_predators_in_birmmals[species]])
                countt += 1
            end
        end
    end
    if info
    @info "Number of presences of $species: $countt / $(length(idx))"
    end
    return countt
end
# Example
# count_number_of_presences("Mustela putorius")