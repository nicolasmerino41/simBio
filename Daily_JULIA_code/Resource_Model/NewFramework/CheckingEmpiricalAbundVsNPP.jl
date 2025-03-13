begin
    protected = true
    ind = protected ? protected_idx : idx
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1])
    
    npp_points = []
    for cell in ind
        push!(npp_points, npp_DA_relative_to_1000[cell[1], cell[2]])
    end
    pi_points = []
    for cell in ind
        value = 0
        for herb_name in herbivore_names
        ind = species_dict_herbivores_in_birmmals[herb_name]
        value += DA_birmmals_with_pi_corrected[cell[1], cell[2]].a[ind]
        end
        push!(pi_points, value)
    end
    
    scatter!(ax, npp_points, pi_points)
    if !protected
        ax2 = Axis(fig[1, 2])
        pi_points_all = []
        for cell in ind
            push!(pi_points_all, DA_birmmals_with_pi_corrected[cell[1], cell[2]].b)
        end
        scatter!(ax2, npp_points, pi_points_all; markersize = 10, color = :red)
        display(fig)
    else
        display(fig)
    end
end

maximum(npp_DA[.!isnan.(npp_DA)])
