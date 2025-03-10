begin
    npp_points = []
    pi_points = []
    for cell in 1:length(idx)
        i, j = idx[cell][1], idx[cell][2]
        if !isnan(npp_DA[i, j])
            push!(npp_points, npp_DA[i, j])
            push!(pi_points, DA_birmmals_with_pi_corrected[i, j].b)
        end
    end
    
    fig = Figure(resolution = (600, 400))
    ax = Axis(
        fig[1, 1], title = "NPP vs. BIOMASS",
        xlabel = "NPP", ylabel = "BIOMASS"
        )
    scatter!(ax, npp_points, pi_points)
    display(fig)
end

begin

    npp_points = []
    pi_points = []
    pi_points_pred = []
    for cell in 1:length(idx)
        i, j = idx[cell][1], idx[cell][2]
        pi_points_for_a_cell = []
        pi_points_for_a_cell_pred = []
        if !isnan(npp_DA[i, j])
            push!(npp_points, npp_DA[i, j])
            a_element = DA_birmmals_with_pi_corrected[i, j].a
            for herb_name in herbivore_names
                herb_index = species_dict_herbivores_in_birmmals[herb_name]
                if herb_index in bad_indices 
                else
                    push!(pi_points_for_a_cell, a_element[herb_index])
                end
            end
            for pred_name in predator_names
                pred_index = species_dict_predators_in_birmmals[pred_name]
                push!(pi_points_for_a_cell_pred, a_element[pred_index])
            end
            push!(pi_points, sum(pi_points_for_a_cell))
            push!(pi_points_pred, sum(pi_points_for_a_cell_pred))
        end
    end
    
    fig = Figure(resolution = (600, 400))
    ax = Axis(
        fig[1, 1], title = "NPP vs. BIOMASS (HERBIVORES ONLY)",
        xlabel = "NPP", ylabel = "BIOMASS"
        )
    scatter!(ax, npp_points, pi_points)
    ax2 = Axis(
        fig[1, 2], title = "NPP vs. BIOMASS (PREDATORS ONLY)",
        xlabel = "NPP", ylabel = "BIOMASS"
    )
    scatter!(ax2, npp_points, pi_points_pred)
    display(fig)
end

herbivore_biomass_df = birmmals_biomass_fixed[map(s -> s in herbivore_names, birmmals_biomass_fixed[:, :species]), :]
A = sort(birmmals_biomass_fixed, :biomass, rev=true)
herbivore_biomass_df = sort(herbivore_biomass_df, :biomass, rev=true)
top_herbivore_names = herbivore_biomass_df.species[1:10]
bad_indices = []
for i in top_herbivore_names
    push!(bad_indices, species_dict_herbivores_in_birmmals[i])
end

"""
    plot_npp_biomass(cell_indices; plot_category=:all, cap_top=0)

Loop over the cells given by `cell_indices` and for each cell extract the NPP value and the 
total biomass of herbivores and predators. If `cap_top > 0`, then the top `cap_top` heaviest 
herbivore species (as determined from the global DataFrame `birmmals_biomass_fixed`) are ignored 
(i.e. their contribution is not included in the herbivore biomass sum).

Arguments:
- `cell_indices::Vector{Int}`: indices of cells to process (e.g. 1:length(idx))
- `plot_category::Symbol`: one of `:herbivores`, `:predators`, `:total`, or `:all` (default is `:all`)
- `cap_top::Int`: number of top-heaviest herbivore species to remove (default is 0)

Returns a DataFrame with columns:
- `cell_id`, `NPP`, `herb_biomass`, `pred_biomass`, and `total_biomass`.

Also, displays a Makie figure according to the chosen plot category.
"""
function plot_npp_biomass(idx; plot_category = :all, cap_top::Int = 0)
    # Determine which herbivore indices to remove if cap_top > 0.
    bad_indices = Int[]
    if cap_top > 0
        # Select rows in birmmals_biomass_fixed corresponding to herbivores.
        herbivore_biomass_df = birmmals_biomass_fixed[map(s -> s in herbivore_names, birmmals_biomass_fixed[:, :species]), :]
        # Sort in descending order by biomass.
        herbivore_biomass_df = sort(herbivore_biomass_df, :biomass, rev=true)
        ncap = min(cap_top, nrow(herbivore_biomass_df))
        top_herbivore_names = herbivore_biomass_df.species[1:ncap]
        for name in top_herbivore_names
            push!(bad_indices, species_dict_herbivores_in_birmmals[name])
        end
    end

    # Initialize arrays for results.
    npp_points = Float64[]
    herb_points = Float64[]
    pred_points = Float64[]
    total_points = Float64[]
    cell_ids = Int[]

    # Loop over each cell.
    for cell in 1:length(idx)
        i, j = idx[cell][1], idx[cell][2]
        # Initialize per-cell biomass vectors.
        herb_biomass_vec = Float64[]
        pred_biomass_vec = Float64[]
        # Process only if NPP data is valid.
        if !isnan(npp_DA[i, j])
            push!(npp_points, npp_DA[i, j])
            push!(cell_ids, cell)
            a_element = DA_birmmals_with_pi_corrected[i, j].a
            # Sum contributions from herbivores (skip those in bad_indices if cap_top > 0).
            for herb_name in herbivore_names
                herb_index = species_dict_herbivores_in_birmmals[herb_name]
                if cap_top > 0 && (herb_index in bad_indices)
                    # Skip this species.
                    continue
                else
                    push!(herb_biomass_vec, a_element[herb_index])
                end
            end
            # Sum contributions from predators.
            for pred_name in predator_names
                pred_index = species_dict_predators_in_birmmals[pred_name]
                push!(pred_biomass_vec, a_element[pred_index])
            end
            # Compute the cell totals.
            push!(herb_points, sum(herb_biomass_vec))
            push!(pred_points, sum(pred_biomass_vec))
            push!(total_points, sum(herb_biomass_vec) + sum(pred_biomass_vec))
        end
    end

    # Build DataFrame with the results.
    df = DataFrame(cell_id = cell_ids, NPP = npp_points,
                   herb_biomass = herb_points, pred_biomass = pred_points,
                   total_biomass = total_points)

    # Plot according to the chosen category.
    if plot_category == :herbivores
        fig = Figure(resolution = (600, 400))
        ax = Axis(fig[1, 1], title = "NPP vs. Herbivore Biomass",
                  xlabel = "NPP", ylabel = "Herbivore Biomass")
        scatter!(ax, npp_points, herb_points)
        display(fig)
    elseif plot_category == :predators
        fig = Figure(resolution = (600, 400))
        ax = Axis(fig[1, 1], title = "NPP vs. Predator Biomass",
                  xlabel = "NPP", ylabel = "Predator Biomass")
        scatter!(ax, npp_points, pred_points)
        display(fig)
    elseif plot_category == :total
        fig = Figure(resolution = (600, 400))
        ax = Axis(fig[1, 1], title = "NPP vs. Total Biomass",
                  xlabel = "NPP", ylabel = "Total Biomass")
        scatter!(ax, npp_points, total_points)
        display(fig)
    elseif plot_category == :all
        fig = Figure(resolution = (800, 600))
        ax1 = Axis(fig[1, 1], title = "NPP vs. Herbivore Biomass",
                   xlabel = "NPP", ylabel = "Herbivore Biomass")
        scatter!(ax1, npp_points, herb_points)
        ax2 = Axis(fig[1, 2], title = "NPP vs. Predator Biomass",
                   xlabel = "NPP", ylabel = "Predator Biomass")
        scatter!(ax2, npp_points, pred_points)
        ax3 = Axis(fig[2, 1], title = "NPP vs. Total Biomass",
                   xlabel = "NPP", ylabel = "Total Biomass")
        scatter!(ax3, npp_points, total_points)
        display(fig)
    else
        @warn "Unknown plot_category: $plot_category. No plot will be generated."
    end

    return df
end

plot_npp_biomass(idx; plot_category=:herbivores, cap_top=10)
