############# Species Abundance Distribution ################
function SAD_plotting(array_output, position; num_samples=10, modified = false, caca = false, log_scale = false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
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
SAD_plotting(a, 1; modified = true, caca = false, log_scale = true)

# Helper function to calculate AIC
function aic_loglikelihood(data, distribution)
    log_likelihood = sum(logpdf(distribution, data))
    k = length(fieldnames(typeof(distribution)))  # Number of parameters in the distribution
    aic = 2 * k - 2 * log_likelihood
    return aic
end

# Main analysis function with Log-Normal, Gamma, and Normal distributions
function analyze_SAD(array_output, position; num_samples=10, modified = false, caca = false)
    # Initialize combined_abundances depending on the conditions provided
    if !modified && !caca
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
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
analyze_SAD(a, 1; num_samples = 100, modified = true, caca = false)
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
rsd = calculate_RSD(a, 1; modified = true, caca = false)
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

# Create a figure
fig = Figure(resolution = (800, 400))

# Log-Normal distribution with μ = 0 and σ = 0.2 (Continuous decline)
ax = Axis(fig[1, 1], title = "Log-Normal Distribution (μ = 0, σ = 0.2)", xlabel = "x", ylabel = "Density")
x_vals = range(0.1, 20, length = 100)  # Small range to show the decline
lognormal_dist = LogNormal(2.69, 0.24)
lines!(ax, x_vals, pdf(lognormal_dist, x_vals), color=:red)

fig