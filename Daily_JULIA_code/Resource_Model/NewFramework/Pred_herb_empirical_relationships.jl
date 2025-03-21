begin
    pred_herb_ratios = Float64[]

    for cell in 1:length(idx)
        i, j = idx[cell][1], idx[cell][2]
        sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[i, j])
        herbs_biomass = 0.0
        for herb_name in herbivore_names
            if herb_name in sp_nm
                herb_index = species_dict_herbivores_in_birmmals[herb_name]
                herbs_biomass += DA_birmmals_with_pi_corrected[i, j].a[herb_index]
            end
        end
        preds_biomass = 0.0
        for pred_name in predator_names
            if pred_name in sp_nm
                pred_index = species_dict_predators_in_birmmals[pred_name]
                preds_biomass += DA_birmmals_with_pi_corrected[i, j].a[pred_index]
            end
        end
        ratio = preds_biomass / herbs_biomass
        push!(pred_herb_ratios, ratio)
    end

    df = DataFrame(
        cell = 1:length(idx),
        pred_herb_ratios = pred_herb_ratios
    )

    df = sort(df, :pred_herb_ratios, rev = true)

    fig = Figure(resolution = (700, 450))
    ax = Axis(
        fig[1, 1],
        title = "Sorted predator / herbivore biomass ratio across cells",
        xlabel = "Cell",
        ylabel = "P/H Biomass Ratio",
        yticks = 0:1.0:11
    )

    MK.scatter!(ax, 1:nrow(df), df.pred_herb_ratios)

    hlines!(ax, [1.0], color = :red, linewidth = 2)
    hlines!(ax, [0.1], color = :orange, linewidth = 2)

    display(fig)
end
