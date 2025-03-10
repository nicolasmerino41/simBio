begin
    
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1])
    
    npp_points = []
    for cell in idx
        push!(npp_points, npp_DA_relative_to_1000[cell[1], cell[2]])
    end
    pi_points = []
    for cell in idx
        value = 0
        for herb_name in herbivore_names
        ind = species_dict_herbivores_in_birmmals[herb_name]
        value += DA_birmmals_with_pi_corrected[cell[1], cell[2]].a[ind]
        end
        push!(pi_points, value)
    end
    
    scatter!(ax, npp_points, pi_points)

    ax2 = Axis(fig[1, 2])
    pi_points_all = []
    for cell in idx
        push!(pi_points_all, DA_birmmals_with_pi_corrected[cell[1], cell[2]].b)
    end
    scatter!(ax2, npp_points, pi_points_all; markersize = 10, color = :red)
    display(fig)
end

maximum(npp_DA[.!isnan.(npp_DA)])
