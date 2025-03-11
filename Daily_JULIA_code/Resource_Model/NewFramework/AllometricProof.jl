using Statistics, Optim, CSV, DataFrames

# Function to compute the error for a single cell given an alpha value.
function cell_error(cell::Int, alpha::Float64; 
    only_protected_cells::Bool = false, 
    cap_heaviest_species::Int = 0
)
    # Get the cell indices: if only_protected_cells is true, use protected_idx; otherwise, use idx.
    ind = only_protected_cells ? protected_idx : idx
    if cell > length(ind)
        error("You specified more cells than there are in the ind of choice.\n
        Non-protected idx has length 5950 while protected idx has length 145.")
    end
    i, j = ind[cell][1], ind[cell][2]

    # Retrieve NPP for this cell.
    cell_NPP = npp_DA_relative_to_1000[i, j]
    if isnan(cell_NPP)
    return missing  # Skip cell if NPP is not available.
    end

    # Extract species names present in the cell.
    sp_names = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[i, j])

    # Filter for herbivores only.
    herb_sp = [sp for sp in sp_names if sp in herbivore_names]
    if isempty(herb_sp)
    return missing
    end

    # If cap_heaviest_species > 0, remove the top n heaviest herbivores from the full birmmals_biomass_fixed.
    if cap_heaviest_species > 0
    herbivore_biomass_df = filter(row -> row.species in herbivore_names, birmmals_biomass_fixed)
    herbivore_biomass_df = sort(herbivore_biomass_df, :biomass, rev=true)
    # Get the names of the top heaviest species.
    top_herbivore_names = herbivore_biomass_df.species[1:cap_heaviest_species]
    # Remove them from the herbivore list for this cell.
    herb_sp = [sp for sp in herb_sp if !(sp in top_herbivore_names)]
    end

    # Retrieve observed biomass and body mass data for the remaining herbivores.
    df_herb = filter(row -> row.species in herb_sp, birmmals_biomass_fixed)
    # Reorder rows so that they match the order in herb_sp.
    df_herb_sorted = DataFrame()
    for sp in herb_sp
    row_sp = filter(row -> row.species == sp, df_herb)
    append!(df_herb_sorted, row_sp)
    end

    # Extract observed biomass (H_obs) and body mass (B).
    H_obs = df_herb_sorted.biomass
    B = df_herb_sorted.bodyMass
    H_tot = sum(H_obs)
    if H_tot == 0
    return missing
    end
    # Compute relative abundances.
    h_vec = H_obs ./ H_tot

    # Compute predicted growth rates based on allometric scaling:
    # g_pred_i = (B_i / mean(B))^(–α)
    B_norm = B ./ mean(B)
    g_pred = B_norm .^ (-alpha)
    avg_g_pred = sum(g_pred .* h_vec)

    # Compute the observed per-unit-biomass production.
    obs_avg = cell_NPP / H_tot

    # Return the absolute error between prediction and observation.
    return abs(avg_g_pred - obs_avg)
end

# Function to compute the average error across all cells for a given alpha.
function global_error(n_cells::Int, alpha::Float64; only_protected_cells = false, cap_heaviest_species = 0)
    # Preallocate an array to hold the error for each cell.
    local_errs = Vector{Union{Missing, Float64}}(undef, n_cells)
    
    Threads.@threads for cell in 1:n_cells
        local_errs[cell] = cell_error(cell, alpha; only_protected_cells = only_protected_cells, cap_heaviest_species = cap_heaviest_species)
    end

    # Filter out any missing errors.
    valid_errs = [err for err in local_errs if err !== missing]
    return isempty(valid_errs) ? NaN : mean(valid_errs)
end

# Define a grid of alpha values to explore.
alpha_range = 0.0:0.01:2.0
errors = [global_error(5950, alpha; only_protected_cells = false, cap_heaviest_species = 5) for alpha in alpha_range]

# Find the alpha that minimizes the average error.
min_error, idx_min = findmin(errors)
optimal_alpha = alpha_range[idx_min]
println("Optimal alpha = ", optimal_alpha, " with average error = ", min_error)

# (Optional) Save results to CSV for further inspection.
df_results = DataFrame(alpha = alpha_range, error = errors)
CSV.write("optimal_alpha_results.csv", df_results)

# Initialize arrays to store NPP and total herbivore biomass for each cell.
npp_values = Float64[]
Htot_values = Float64[]

# Loop over all cells.
only_protected_cells = true
ind = only_protected_cells ? protected_idx : idx
for cell in 1:length(ind)
    i, j = ind[cell][1], ind[cell][2]
    
    # Get the NPP value for the cell.
    cell_npp = npp_DA_relative_to_1000[i, j]
    if isnan(cell_npp)
        continue  # Skip cells with NaN NPP.
    end
    
    # Get the abundance vector for the cell.
    cell_abundances = DA_birmmals_with_pi_corrected[i, j].a
    
    # Get indices corresponding to herbivores.
    herb_indices = [species_dict_herbivores_in_birmmals[sp] for sp in herbivore_names if haskey(species_dict_herbivores_in_birmmals, sp)]
    
    # Sum the herbivore abundances.
    H_tot = sum(cell_abundances[herb_indices])
    
    # Only consider cells where the total herbivore biomass is positive.
    if H_tot > 0
        push!(npp_values, cell_npp)
        push!(Htot_values, H_tot)
    end
end

# Compute the per-cell ratio NPP/H_tot and then average over cells.
ratios = npp_values ./ Htot_values
# sub = filter!(x->x>100.0, copy(ratios))
avg_ratio = mean(ratios)

# Plot the distribution of ratios.
fig = Figure(resolution = (600, 400))
ax = Axis(fig[1, 1], title = "Distribution of NPP/H_tot ratios", xlabel = "NPP/H_tot", ylabel = "Frequency")
hist!(ax, ratios, bins = 50, normalization = :pdf)
display(fig)

println("The average NPP/H_tot across cells is ", avg_ratio)

########################
function plot_predicted_vs_observed(n_cells::Int, alpha::Float64; only_protected_cells = false)
    predicted = Float64[]
    observed  = Float64[]

    ind = only_protected_cells ? protected_idx : idx
    if n_cells > length(ind)
        error("You specified more cells than there are in the ind of choice.\n
        Non-protected idx has length 5950 while protected idx has length 145.")
    end
    for cell in 1:n_cells
        # Retrieve cell indices.
        i, j = ind[cell][1], ind[cell][2]
        # Get the NPP value for the cell.
        cell_NPP = npp_DA_relative_to_1000[i, j]
        if isnan(cell_NPP)
            continue
        end

        # Extract species present in the cell.
        sp_names = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[i, j])
        # Filter to herbivores.
        herb_sp = [sp for sp in sp_names if sp in herbivore_names]
        if isempty(herb_sp)
            continue
        end

        # Get observed biomass and body mass data from birmmals_biomass_fixed.
        df_herb = filter(row -> row.species in herb_sp, birmmals_biomass_fixed)
        # Reorder rows to match the order in herb_sp.
        df_herb_sorted = DataFrame()
        for sp in herb_sp
            row_sp = filter(row -> row.species == sp, df_herb)
            append!(df_herb_sorted, row_sp)
        end

        H_obs = df_herb_sorted.biomass
        B = df_herb_sorted.bodyMass
        H_tot = sum(H_obs)
        if H_tot == 0
            continue
        end
        # Compute relative abundances.
        h_vec = H_obs ./ H_tot

        # Compute predicted growth rates using the scaling law:
        # g_pred = (B_i/mean(B))^(–alpha)
        B_norm = B ./ mean(B)
        g_pred = B_norm .^ (-alpha)
        avg_pred = sum(g_pred .* h_vec)

        push!(predicted, avg_pred)
        push!(observed, cell_NPP / H_tot)
    end

    # Create a scatter plot with Makie.
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1],
              title = "Predicted vs. Observed Average Growth",
              xlabel = "Predicted ⟨g⟩",
              ylabel = "Observed (NPP/H_tot)")
    scatter!(ax, predicted, observed, color = :blue, markersize = 8)

    # Add the identity line.
    all_vals = vcat(predicted, observed)
    min_val = minimum(all_vals)
    max_val = maximum(all_vals)
    lines!(ax, [min_val, max_val], [min_val, max_val], color = :red, linewidth = 2)

    display(fig)
end

# Example usage:
# Assume 'optimal_alpha' was computed previously.
plot_predicted_vs_observed(145, optimal_alpha; only_protected_cells = true)

