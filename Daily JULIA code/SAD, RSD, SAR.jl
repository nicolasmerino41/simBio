############# Species Abundance Distribution ################
function SAD_plotting(array_output, position, type = nothing; num_samples=10, modified = false, caca = false, log_scale = false)
        
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
        if type == "abundance"
            # Ensure element-wise division by body_mass_vector for each cell
            for i in idx
                combined_abundances[i] = MyStructs256(SVector{256, Float64}(combined_abundances[i].a ./ body_mass_vector))
            end
        end
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    end

    # Sample 10 random cells
    random_cells = sample(idx, num_samples, replace=false)

    # Prepare the figure to hold all subplots; adjust the figure size based on the number of valid cells
    fig = Figure(resolution=(800, max(800, length(random_cells) * 200)))

    # Process each cell
    for (i, cell) in enumerate(random_cells)
        if !any(isnan, combined_abundances[cell].a)
            # Order the species by abundance in descending order
            if log_scale
                ordered_abundances = sort(log.(combined_abundances[cell].a), rev=true)
            else
                ordered_abundances = sort(combined_abundances[cell].a, rev=true)
            end
            
            # Generate species ranks (x-values)
            species_ranks = 1:length(ordered_abundances)

            # Create a subplot for this cell
            ax = Axis(fig[i, 1], title="Cell $cell", xlabel="Species Rank", ylabel="Abundance")
            barplot!(ax, species_ranks, ordered_abundances)

            # Dynamically set the x-axis limits based on the number of species
            xlims!(ax, 0, length(ordered_abundances) + 1)  # +1 to ensure the last bar is fully visible
            ylims!(ax, 0, maximum(ordered_abundances) * 1.1)
        else
            # Create a placeholder axis if data is not valid
            ax = Axis(fig[i, 1], title="Cell $cell - Data Not Available")
        end
    end
    
    # Display the figure
    fig
end

SAD_plotting(a, 1, "abundance"; modified = true, caca = false, log_scale = false)

# Helper function to calculate AIC
function aic_loglikelihood(data, distribution)
    log_likelihood = sum(logpdf(distribution, data))
    k = length(fieldnames(typeof(distribution)))  # Number of parameters in the distribution
    aic = 2 * k - 2 * log_likelihood
    return aic
end

# Main analysis function with Log-Normal, Gamma, and Normal distributions
function analyze_SAD(array_output, position, type = nothing; num_samples=10, modified = false, caca = false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
        if type == "abundance"
            # Ensure element-wise division by body_mass_vector for each cell
            for i in idx
                combined_abundances[i] = MyStructs256(SVector{256, Float64}(combined_abundances[i].a ./ body_mass_vector))
            end
        end
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    end

    # Sample 10 random cells
    sampled_cells = sample(idx, num_samples, replace=false)

    # Aggregate species abundances across sampled cells, considering body mass threshold
    species_counts = Dict{Int, Float64}()
    for cell in sampled_cells
        cell_data = combined_abundances[cell].a
        for (species_id, abundance) in enumerate(cell_data)
            if !isnan(abundance) && abundance > body_mass_vector[species_id]
                species_counts[species_id] = get(species_counts, species_id, 0.0) + abundance
            end
        end
    end

    # Filter out species with zero or one total abundance
    filtered_counts = Dict{Int, Float64}(
        species_id => abundance for (species_id, abundance) in species_counts if abundance > 1
    )

    # Convert the filtered dictionary to a DataFrame for easier manipulation
    if isempty(filtered_counts)
        println("No species meet the criteria for analysis.")
        return
    end
    
    aggregated_data = DataFrame(species_id=collect(keys(filtered_counts)), total_abundance=collect(values(filtered_counts)))

    # Check for variance in the data
    if std(aggregated_data.total_abundance) == 0
        println("Standard deviation is zero, which is not suitable for fitting.")
        return  # Exit the function if data is not suitable for fitting
    end

    # Initialize model fitting results
    aic_values = Dict()

    # Fit Log-Normal distribution to log-transformed abundances
    log_abundances = log.(aggregated_data.total_abundance)
    lognorm_params = nothing
    try
        lognorm_params = fit_mle(LogNormal, log_abundances)
        aic_values["LogNormal"] = aic_loglikelihood(log_abundances, LogNormal(lognorm_params.μ, lognorm_params.σ))
        println("Fitted LogNormal parameters: μ = $(lognorm_params.μ), σ = $(lognorm_params.σ)")
    catch e
        println("Error fitting LogNormal: ", e)
    end

    # Fit Gamma distribution
    gamma_params = nothing
    try
        gamma_params = fit_mle(Gamma, aggregated_data.total_abundance)
        aic_values["Gamma"] = aic_loglikelihood(aggregated_data.total_abundance, Gamma(gamma_params.α, gamma_params.θ))
        println("Fitted Gamma parameters: α = $(gamma_params.α), θ = $(gamma_params.θ)")
    catch e
        println("Error fitting Gamma: ", e)
    end

    # Fit Normal distribution
    normal_params = nothing
    try
        normal_params = fit_mle(Normal, aggregated_data.total_abundance)
        aic_values["Normal"] = aic_loglikelihood(aggregated_data.total_abundance, Normal(normal_params.μ, normal_params.σ))
        println("Fitted Normal parameters: μ = $(normal_params.μ), σ = $(normal_params.σ)")
    catch e
        println("Error fitting Normal: ", e)
    end

    # Find the best model based on AIC
    if isempty(aic_values)
        println("No models were successfully fitted.")
        return
    end

    best_model = argmin(aic_values)
    println("Best model: $best_model with AIC = $(aic_values[best_model])")

    # Visualize with Makie
    fig = Figure(size = (800, 400))
    ax1 = Axis(fig[1, 1], title="Aggregated SAD with Best Fit", xlabel="Total Abundance", ylabel="Frequency")

    # Plot the histogram using Makie
    # histogram_data = aggregated_data.total_abundance
    histogram_data = best_model == "LogNormal" ? log.(aggregated_data.total_abundance) : aggregated_data.total_abundance
    hist_plot = hist!(ax1, histogram_data, bins=50, color=:blue)  # Increased bins to 100 for more granularity

    # Extract the total number of frequencies (weights) manually from the histogram
    hist_frequencies = hist_plot[1][]  # Extract the frequencies (weights) observable
    
    if best_model == "LogNormal" && !isnothing(lognorm_params)
        x_vals = range(minimum(log.(aggregated_data.total_abundance)), stop=maximum(log.(aggregated_data.total_abundance)), length=100)
        y_vals = pdf(LogNormal(lognorm_params.μ, lognorm_params.σ), x_vals)
        # Overlay the best-fitting model on top of the histogram
        
    elseif best_model == "Gamma" && !isnothing(gamma_params)
        y_vals = pdf(Gamma(gamma_params.α, gamma_params.θ), x_vals)
        # Overlay the best-fitting model on top of the histogram
        x_vals = range(minimum(aggregated_data.total_abundance), stop=maximum(aggregated_data.total_abundance), length=100)
    elseif best_model == "Normal" && !isnothing(normal_params)
        # Ensure the x-axis range covers the full data spread for Normal
        x_vals = range(minimum(histogram_data) - 5 * normal_params.σ, stop=maximum(histogram_data) + 5 * normal_params.σ, length=100)
        y_vals = pdf(Normal(normal_params.μ, normal_params.σ), x_vals)
    else
        println("No valid model selected for plotting.")
        return
    end

    # Ensure correct scaling by normalizing the curve to match the histogram's area
    hist_area = sum(hist_frequencies) * mean(diff(histogram_data))  # Area under the histogram
    y_vals_scaled = y_vals * hist_area  # Scale the PDF to match histogram area
    
    lines!(ax1, x_vals, y_vals, color=:red, linewidth=2)

    fig
end

# Call the function to analyze SAD
analyze_SAD(a, 1, "abundance"; num_samples = 100, modified = true, caca = false)
############# Range Size Distribution ################
function calculate_RSD(array_output, position; modified = false, caca = false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    end
    
    # Number of species and cells in the grid
    num_species = length(combined_abundances[1].a)
    num_cells = length(idx)

    # Initialize an array to store the range size for each species
    range_sizes = zeros(Int, num_species)
    
    # Loop through all cells
    for cell in idx
        cell_abundances = combined_abundances[cell].a
        
        # For each species, check if its abundance exceeds its threshold from body_mass_vector
        for species in 1:num_species
            if cell_abundances[species] > body_mass_vector[species]
                range_sizes[species] += 1  # Increment the range size if species is present in this cell
            end
        end
    end

    # Return the range size distribution as a dictionary (frequency of each range size)
    range_size_distribution = Dict{Int, Int}()
    
    for range_size in range_sizes
        range_size_distribution[range_size] = get(range_size_distribution, range_size, 0) + 1
    end
    
    return range_size_distribution
end
rsd = calculate_RSD(b, 1; modified = true, caca = false)
println("Range Size Distribution: ", rsd)

function plot_RSD(rsd)
    fig = Figure(resolution=(800, 400))
    ax = Axis(fig[1, 1], title="Range Size Distribution", xlabel="Range Size (Number of Cells)", ylabel="Frequency")
    
    # Convert dictionary to vectors for plotting
    range_sizes = collect(keys(rsd))
    frequencies = collect(values(rsd))
    
    MK.barplot!(ax, range_sizes, frequencies, color=:blue, width=10.0)
    
    fig
end

plot_RSD(rsd)
############# Species Area Relationship ################
####### RANDOM SAMPLING ########
function random_SAR(array_output, position; max_cells=100, step_size=5, modified=false, caca=false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    end

    valid_cells = idx  # Or use `valid_cells = findall(DA_sum .== 1)` for valid cells from DA_sum
    num_valid_cells = length(valid_cells)
    max_cells = min(max_cells, num_valid_cells)

    # Store species richness for each area (i.e., number of cells sampled)
    species_area = Dict{Int, Float64}()

    # Loop over increasing numbers of sampled cells
    for area in 1:step_size:max_cells
        sampled_cells = sample(valid_cells, area, replace=false)

        # Set to track which species are present across the sampled cells
        present_species = Set{Int}()

        # Loop over the sampled cells and check species presence
        for cell in sampled_cells
            cell_abundances = combined_abundances[cell].a
            for species in 1:size(cell_abundances, 1)
                if cell_abundances[species] > body_mass_vector[species]
                    push!(present_species, species)
                end
            end
        end

        # Store the species richness (number of unique species) for the sampled area
        species_area[area] = length(present_species)
    end

    return species_area
end

sar = calculate_SAR(a, 1; max_cells=100, step_size=5, threshold=0.0, modified=true, caca=false)
println("Species Area Relationship: ", sar)
######## CONTIGUOUS SAMPLING ########
# Function for contiguous sampling SAR
function contiguous_sar(array_output, position; max_cells=100, threshold=0.0, grid_size=(10, 10), num_origins=1, modified=false, caca=false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    end

    valid_cells = idx  # Assuming `idx` contains valid CartesianIndex{2} cells.
    num_valid_cells = length(valid_cells)
    max_cells = min(max_cells, num_valid_cells)

    # For multiple origin points, store species richness at each area in parallel sequences
    all_richness = []

    # Define possible neighbor offsets (up, down, left, right, diagonals)
    neighbor_offsets = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)]

    # Repeat for multiple origins
    for _ in 1:num_origins
        # Start from a random valid CartesianIndex{2}
        origin_cell = sample(valid_cells, 1)[1]
        sampled_cells = Set([origin_cell])
        present_species = Set{Int}()

        # Track species richness at each step for this origin
        species_area = []

        # Add species from the origin cell
        origin_abundances = combined_abundances[origin_cell].a
        for species in 1:length(origin_abundances)
            if origin_abundances[species] > body_mass_vector[species]
                push!(present_species, species)
            end
        end

        # Record richness after the first cell
        push!(species_area, length(present_species))

        # Progressively add adjacent cells one by one
        while length(sampled_cells) < max_cells
            # Pick a random already sampled cell to expand from
            current_cell = sample(collect(sampled_cells), 1)[1]
            
            # Pick a random direction (up, down, left, right, diagonal)
            dr, dc = sample(neighbor_offsets, 1)[1]
            neighbor_cell = CartesianIndex(Tuple(current_cell) .+ (dr, dc))

            # Only add neighbor if it's valid and hasn't been sampled yet
            if neighbor_cell in valid_cells && neighbor_cell ∉ sampled_cells
                push!(sampled_cells, neighbor_cell)

                # Add species from the new cell
                cell_abundances = combined_abundances[neighbor_cell].a
                for species in 1:length(cell_abundances)
                    if cell_abundances[species] > body_mass_vector[species]
                        push!(present_species, species)
                    end
                end

                # Record the species richness after adding the new cell
                push!(species_area, length(present_species))
            end
        end

        # Add the species richness sequence for this origin to the overall list
        push!(all_richness, species_area)
    end

    return all_richness
end

sar_contiguous = contiguous_sar(c, 1; modified = true, caca = false, num_origins = 10, max_cells = 1000)

function plot_SAR(all_richness)
    fig = Figure(resolution=(800, 400))
    ax = Axis(fig[1, 1], title="Species-Area Relationship", xlabel="Number of Cells", ylabel="Species Richness")

    # Ensure all values in all_richness are numeric (convert to Float64 if needed)
    for richness in all_richness
        x_vals = 1:length(richness)
        
        # Ensure that richness is a vector of numeric types
        numeric_richness = Float64.(richness)  # Convert all elements to Float64

        # Plot each origin's species richness progression as a separate line
        lines!(ax, x_vals, numeric_richness, linewidth=2)
    end

    fig
end

# Assuming sar_contiguous is defined and contains numeric data
plot_SAR(sar_contiguous)

############# CALCULATE SAR ################
using Statistics

function calculate_sar(all_richness)
    slopes = []
    p_values = []
    constants = []

    # Loop over each origin's SAR progression
    for richness in all_richness
        # Number of cells sampled (x-values)
        x_vals = 1:length(richness)

        # Log-transform both species richness and number of cells (to fit the power law model)
        log_x_vals = log.(x_vals)
        log_richness = log.(richness)

        # Perform linear regression on the log-transformed values using least squares
        n = length(log_x_vals)
        x_mean = mean(log_x_vals)
        y_mean = mean(log_richness)
        
        # Calculate slope (z) using covariance/variance
        numerator = sum((log_x_vals .- x_mean) .* (log_richness .- y_mean))
        denominator = sum((log_x_vals .- x_mean) .^ 2)
        z = numerator / denominator  # This is the exponent in the power law (slope)

        # Calculate intercept (log(c))
        log_c = y_mean - z * x_mean
        
        # Recover the constant c from log(c)
        c = exp(log_c)
        
        # Estimate the residuals
        residuals = log_richness .- (log_c .+ z .* log_x_vals)
        
        # Calculate the standard error of the slope
        s_squared = sum(residuals.^2) / (n - 2)
        standard_error = sqrt(s_squared / denominator)

        # Calculate t-statistic for the slope
        t_statistic = z / standard_error

        # Calculate p-value from the t-statistic
        p_value = 2 * (1 - cdf(TDist(n - 2), abs(t_statistic)))

        # Append slope, p-value, and constant to the results
        push!(slopes, z)
        push!(p_values, p_value)
        push!(constants, c)
    end

    # Calculate the average slope (z), p-value, and constant (c) across all origins
    avg_z = mean(slopes)
    avg_p_value = mean(p_values)
    avg_c = mean(constants)

    return avg_z, avg_p_value, avg_c
end

# Calculate the average slope and p-value for the SAR
avg_slope, avg_p_value, avg_c = calculate_sar(sar_contiguous)

# Print the results
println("Average Slope: ", avg_slope)
println("Average p-value: ", avg_p_value)

########################################################
# Create a figure
fig = Figure(resolution = (800, 400))

# Log-Normal distribution with μ = 0 and σ = 0.2 (Continuous decline)
ax = Axis(fig[1, 1], title = "Log-Normal Distribution (μ = 0, σ = 0.2)", xlabel = "x", ylabel = "Density")
x_vals = range(0.1, 20, length = 100)  # Small range to show the decline
lognormal_dist = LogNormal(2.69, 0.24)
lines!(ax, x_vals, pdf(lognormal_dist, x_vals), color=:red)

fig

#########################################################
######### BIOMASS PYRAMID PLOTTING ######################
#########################################################
#########################################################
# List all .jls files in the 'outputs' folder
using Glob
jls_files = sort(glob("outputs/*.jls"))
jls_filess = sample(jls_files, 20, replace = false)
function biomass_distribution_plotting(
    array_output,
    position,
    type = "region";
    cell = nothing,
    bin_size = 0.2,
    caca = false,
    ax = nothing,
    logscale = false
)
    if type == "region"
        # Combine abundances for the entire region
        if caca
            combined_abundances = array_output .* lambda_DA[position]
        else
            combined_abundances = array_output[end].state .* lambda_DA[position]
        end
        biomass = MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
        for i in idx
            biomass = biomass + combined_abundances[i]
        end
        biomass = biomass.a
    elseif type == "cell"
        if isnothing(cell)
            error("Please provide a valid position for cell mode.")
        end
        if caca
            combined_abundances = array_output[cell[1], cell[2]] * lambda_DA[position][cell[1], cell[2]]
        else
            combined_abundances = deepcopy(array_output[end].state[cell[1], cell[2]]) * lambda_DA[position][cell[1], cell[2]]
        end
        biomass = combined_abundances.a
    else
        error("Invalid type provided. Use 'region' for the whole region or 'cell' for a single cell.")
    end

    # Define the range of trophic levels
    min_trophic_level = floor(minimum(TrophInd_vector) / bin_size) * bin_size
    max_trophic_level = ceil(maximum(TrophInd_vector) / bin_size) * bin_size

    # Create bins for trophic levels
    bins = min_trophic_level:bin_size:max_trophic_level

    # Calculate total biomass per bin
    biomass_per_bin = zeros(length(bins) - 1)
    for i in 1:(length(bins) - 1)
        bin_indices = (TrophInd_vector .>= bins[i]) .& (TrophInd_vector .< bins[i + 1])
        biomass_per_bin[i] = sum(biomass[bin_indices])
    end

    # Apply log scale to biomass if requested
    if logscale
        # Add a small epsilon to avoid log(0)
        epsilon = 1e-70
        biomass_per_bin = log10.(biomass_per_bin .+ epsilon)
    end

    # If no axis is provided, create one
    if ax === nothing
        # Create the figure and axis
        fig = Figure(resolution = (600, 400))
        ax = Axis(
            fig[1, 1],
            title = "Biomass Distribution by Trophic Level",
            xlabel = logscale ? "Log10(Total Biomass)" : "Total Biomass",
            ylabel = "Trophic Level"
        )
        # Plot horizontal bars
        barplot!(ax, bins[1:end-1], biomass_per_bin; color = :blue, direction = :x)
        # Display the figure
        fig
    else
        # Plot on the provided axis
        # Adjust axis labels if necessary
        ax.xlabel = logscale ? "Log10(Total Biomass)" : "Total Biomass"
        ax.ylabel = "Trophic Level"
        # Plot horizontal bars
        barplot!(ax, bins[1:end-1], biomass_per_bin; color = :blue, direction = :x)
        # Adjust axis limits if needed
        # xlims!(ax, 0, maximum(biomass_per_bin) * 1.1)
        # ylims!(ax, minimum(bins), maximum(bins))
    end
end

# Number of rows and columns
num_rows = 4
num_columns = 5
# Create a figure with a 4x5 grid layout
fig = Figure(resolution = (1000, 800))
saxes = [Axis(fig[i, j]) for i in 1:num_rows, j in 1:num_columns]
for (idx_file, file) in enumerate(jls_filess)
    # Load the data
    array_output = deserialize(file)
    
    # Extract parameters from the filename
    # Filename format: "3.0-5.0-0.1-1.0-0.33.jls"
    filename = basename(file)
    filename_without_ext = splitext(filename)[1]
    params_str = split(filename_without_ext, "-")
    params = parse.(Float64, params_str)
    sigma, epsilon, alfa, sigma_comp, assymetry = params
    
    # Determine the grid position
    row = div((idx_file - 1), num_columns) + 1
    col = mod(idx_file - 1, num_columns) + 1

    # Check that row and col are within the grid dimensions
    if row > num_rows
        println("Warning: More files than grid positions. Skipping file: ", file)
        continue
    end
    
    ax = saxes[row, col]
    
    # Plot on the axis using your function
    biomass_distribution_plotting(array_output, 1, "cell"; cell = (35, 35), caca = true, ax = ax, logscale = true)
    
    # Set the axis title with parameter values
    ax.title = "σ=$(sigma), ε=$(epsilon)\nα=$(alfa), σc=$(sigma_comp), asy=$(assymetry)"
    xlims!(ax, 0, nothing)
    # Adjust axis labels (optional, since labels might overlap in a grid)
    if row == num_rows
        ax.xlabel = "Total Biomass"
    else
        hidespines!(ax, :b)
        ax.xlabelvisible = false
    end
    if col == 1
        ax.ylabel = "Trophic Level"
    else
        hidespines!(ax, :l)
        ax.ylabelvisible = false
    end
end

######## PLOTTING RICHNESS/BIOMASS MAP FOR n OUTPUTS ########
######################################################
######################################################
function richness_biomass_map_plotting()
type = "richness"
jls_filess = sample(jls_files, 20, replace = false)
# Number of rows and columns
num_rows = 4
num_columns = 5
# Create a figure with a 4x5 grid layout
fig = Figure(resolution = (1000, 800))
saxes = [Axis(fig[i, j]) for i in 1:num_rows, j in 1:num_columns]
for (idx_file, file) in enumerate(jls_filess)
    # Load the data
    array_output = deserialize(file)
    
    # Extract parameters from the filename
    # Filename format: "3.0-5.0-0.1-1.0-0.33.jls"
    filename = basename(file)
    filename_without_ext = splitext(filename)[1]
    params_str = split(filename_without_ext, "-")
    params = parse.(Float64, params_str)
    sigma, epsilon, alfa, sigma_comp, assymetry = params
    
    # Determine the grid position
    row = div((idx_file - 1), num_columns) + 1
    col = mod(idx_file - 1, num_columns) + 1

    # Check that row and col are within the grid dimensions
    if row > num_rows
        println("Warning: More files than grid positions. Skipping file: ", file)
        continue
    end
    
    ax = saxes[row, col]
    
    # Plot on the axis using your function
    if type == "richness"
        MK.image!(ax, array_output; colormap = custom_palette, colorrange = (0, 256))
    elseif type == "biomass"
        Makie.heatmap!(ax, array_output; interpolate=false, colormap=custom_palette, colorrange = (0, m))
    end
    
    # Set the axis title with parameter values
    ax.title = "σ=$(sigma), ε=$(epsilon)\nα=$(alfa), σc=$(sigma_comp), asy=$(assymetry)"
    
    ax.yreversed[] = true
    return fig
end
end
richness_biomass_map_plotting()
########## LOG-BIOMASS VS TROPHIC LEVEL LINEAR REGRESSION ##########
