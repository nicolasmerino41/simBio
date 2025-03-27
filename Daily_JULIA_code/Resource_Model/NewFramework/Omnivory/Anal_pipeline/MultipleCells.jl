using DifferentialEquations, ForwardDiff, LinearAlgebra, StaticArrays
using NaNMath          # for safe log

# --- Helper: extract Float64 from Dual numbers ---
to_float(x) = ForwardDiff.value(x)

# --- PART 1. Find a stable configuration for each cell and record elasticities ---
# Define a list of cell identifiers. (Assume cells are indexed by integers or similar.)
num_cells = 2
cells = [i for i in 1:num_cells]  # adjust as needed

# We will store, for each cell, a tuple containing:
#  • the cell identifier,
#  • the chosen parameter (here, the conversion efficiency ε),
#  • the elasticity matrix (J_elasticity; dimensions: (# species)×3),
#  • and the species names vector (as extracted from the cell).
stable_results = []

# For each cell, we search over ε (third parameter) until a stable configuration is found.
for cell in cells
    stable_found = false
    chosen_eps = nothing
    result = nothing
    for eps_val in 0.0:0.01:1.0
        if stable_found; break; end
        for mu in 0.0:0.1:1.0
            if stable_found; break; end
            for mean_m_alpha in 0.0:0.1:1.0
                res = analytical_equilibrium(
                    cell,
                    mu, eps_val, mean_m_alpha;
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = false, pi_size = 1.0,
                    H_init = nothing,
                    P_init = nothing,
                    nu_omni_proportion = 1.0,
                    nu_b_proportion = 1.0,
                    r_omni_proportion = 1.0,
                    callbacks = false
                )
                # Test stability using the Jacobian (local stability if all real parts < 0).
                eigvals = eigen(res.Jacobian).values
                if all(real.(eigvals) .< 0)
                    chosen_mu = mu
                    chosen_eps = eps_val
                    chosen_mean_m_alpha = mean_m_alpha
                    result = res
                    stable_found = true
                    break
                end
            end
        end
    end
    if stable_found
        # Compute the elasticity matrix via com_log_eq.
        # We call com_log_eq with the same parameter values: [0.5, chosen_eps, 0.1].
        J_elast = ForwardDiff.jacobian(p -> com_log_eq(p, cell), [chosen_mu, chosen_eps, chosen_mean_m_alpha])
        # Extract species names from the cell using your extraction function.
        sp_names = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[idx[cell][1], idx[cell][2]])
        push!(stable_results, (cell = cell, eps = chosen_eps, elasticity = J_elast, species = sp_names))
        println("Cell $cell: stable configuration found at mu = $chosen_mu, eps = $chosen_eps, m_alpha = $chosen_mean_m_alpha with $(length(sp_names)) species.")
    else
        println("No stable configuration found for cell $cell")
    end
end

# --- PART 2. Combine elasticity data across cells ---
# We'll build a dictionary mapping each species name to a vector of elasticity row vectors (one per cell where it appears).
# Each elasticity row is a 3-element SVector (for the three parameters).
species_data = Dict{String, Vector{SVector{3,Float64}}}()

for res in stable_results  # iterate over all stable results
    elast = res.elasticity      # (# species) × 3 matrix
    sp_names = res.species      # species names for that cell
    n = min(length(sp_names), size(elast, 1))
    for i in 1:n
        sp = sp_names[i]
        row = SVector{3}(elast[i, :])
        if haskey(species_data, sp)
            push!(species_data[sp], row)
        else
            species_data[sp] = [row]
        end
    end
end

# For each species, compute the mean and standard deviation (elementwise) of its elasticity vectors.
species_avg = Dict{String, SVector{3,Float64}}()
species_std = Dict{String, SVector{3,Float64}}()

for (sp, vecs) in species_data
    n = length(vecs)
    mean_vec = reduce(+, vecs) / n
    # Compute standard deviation elementwise.
    sq_diffs = [ (v .- mean_vec).^2 for v in vecs ]
    std_vec = sqrt.(reduce(+, sq_diffs) / n)
    species_avg[sp] = mean_vec
    species_std[sp] = std_vec
end

# --- PART 3. Plot average elasticity (with error bars) for each parameter ---
# We want one bar plot per parameter.
# The x-axis will show species (sorted by average elasticity for that parameter)
# and each bar is colored according to species type (herbivore, omnivore, or predator).
# Define mapping from species to type:
species_types = Dict{String, String}()
for sp in keys(species_avg)
    species_types[sp] = sp in omnivore_names ? "omnivore" :
                          sp in herbivore_names ? "herbivore" :
                          sp in carnivore_names ? "predator" : "unknown"
end

# Define colors for each type.
type_colors = Dict("herbivore" => :blue, "omnivore" => :green, "predator" => :red, "unknown" => :gray)

# Parameter names (for display on plots).
param_names = ["μ", "ε", "m_α"]

# Get list of all species (from the union of cells).
all_species = collect(keys(species_avg))

for param in 1:3
    # For each species, extract the mean elasticity and std for this parameter.
    mean_vals = [species_avg[sp][param] for sp in all_species]
    std_vals  = [species_std[sp][param] for sp in all_species]
    
    # Sort species by mean elasticity (descending).
    sorted_idx = sortperm(mean_vals, rev=true)
    sorted_species = all_species[sorted_idx]
    sorted_means   = mean_vals[sorted_idx]
    sorted_stds    = std_vals[sorted_idx]
    
    # Determine colors based on species type.
    sorted_colors = [type_colors[species_types[sp]] for sp in sorted_species]
    println("For parameter $(param_names[param]) the elasticity values are: $sorted_means")
    # Create a figure and axis.
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1,1],
        title = "Average Sensitivity to Parameter $(param_names[param])",
        xlabel = "Species (sorted)",
        ylabel = "Elasticity",
        xticklabelrotation = π/4,
        xticklabelalign = (:right, :center)
    )

    # Plot error bars and points.
    errorbars!(ax, 1:length(sorted_species), sorted_means, sorted_stds, color = sorted_colors)
    scatter!(ax, 1:length(sorted_species), sorted_means, color = sorted_colors, markersize = 10)

    # Plot a horizontal reference line at zero.
    lines!(ax, 1:length(sorted_species), fill(0.0, length(sorted_species)), color = :red)

    # Optionally, set the x-tick labels to the species names.
    ax.xticks = (1:length(sorted_species), sorted_species)

    display(fig)
end

