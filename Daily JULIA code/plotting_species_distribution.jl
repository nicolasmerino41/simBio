const bwp = [RGB(0.0, 0.0, 0.0), RGB(1.0, 1.0, 1.0)]  # Black for 0, White for 1
############## REAL species distribution ##############
# Define the species distribution plotting function
function plot_species_distribution(species_name::String)
    # Create a matrix to store species presence/absence
    species_matrix = fill(NaN, size(utmraster_DA))
    
    # Loop through the DataFrame and populate the species matrix
    for row in eachrow(species_df)
        value = row[:Value]
        utm_index = findfirst(x -> x == value, utmraster_DA)
        if utm_index !== nothing
            species_matrix[Tuple(utm_index)...] = row[Symbol(species_name)]
        end
    end
    
    # Plot the species distribution using map_plot
    map_plot(species_matrix; type = "heatmap", palette = :greys, title = "Distribution of $species_name")
end

# Example usage
plot_species_distribution("Alytes cisternasii");

# Define the function to plot species distribution from simulation result
function plot_simulation_species_distribution(array_output, species_name::String)
    species_index = species_to_index[species_name]
    combined_abundances = (deepcopy(array_output[end].birmmals) + deepcopy(array_output[end].herps)).*lambda_DA
    
    # Create a matrix to store species presence/absence
    species_matrix = fill(NaN, size(combined_abundances))
    
    # Loop through the specified indices and populate the species matrix
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances[species_index] > body_mass_vector[species_index]
        species_matrix[cell] = presence ? 1.0 : 0.0
        presence ? println("In cell $cell is present") : println("In cell $cell is absence")
    end
    
    # Force contrast by setting the bottom-left corner
    unique_values = unique(species_matrix[.!isnan.(species_matrix)])  # Get unique non-NaN values
    
    if length(unique_values) == 1
        if unique_values[1] == 0.0
            species_matrix[1, 1] = 1.0  # Set bottom-left corner to 1 if all are 0
        elseif unique_values[1] == 1.0
            species_matrix[1, 1] = 0.0  # Set bottom-left corner to 0 if all are 1
        end
    end
    
    # Plot the species distribution using map_plot with a fixed discrete color palette
    map_plot(species_matrix; type = "heatmap", palette = bw, 
             color=:bwp, title = "Simulation Distribution of $species_name")
end

# Example usage
plot_simulation_species_distribution(p, "Canis lupus")

############ COMBINED species distribution ###########
# Define the function to plot both real and simulated maps side by side
# Function to plot real and simulated species distributions side by side
function combined_species_distribution_plot(array_output, species_name)
    # Create a matrix to store real species presence/absence
    real_species_matrix = fill(NaN, size(utmraster_DA))
    
    # Loop through the DataFrame and populate the real species matrix
    for row in eachrow(species_df)
        value = row[:Value]
        utm_index = findfirst(x -> x == value, utmraster_DA)
        if utm_index !== nothing
            real_species_matrix[utm_index[1], utm_index[2]] = row[species_name]
        end
    end

    # Create a matrix to store simulated species presence/absence
    species_index = species_to_index[species_name]
    # combined_abundances = (deepcopy(array_output[end].birmmals) + deepcopy(array_output[end].herps)).*lambda_DA
    combined_abundances = deepcopy(array_output[end].state).*lambda_DA
    simulated_species_matrix = fill(NaN, size(utmraster_DA))
    
    for cell in idx
        abundances = combined_abundances[cell].a
        presence = abundances .> body_mass_vector[species_index]
        simulated_species_matrix[cell] = presence[species_index] ? 1.0 : 0.0
    end
    
    # Plotting both maps side by side
    fig = Figure(resolution = (800, 400))
    ax1 = Axis(fig[1, 1], title = "Real Distribution of $species_name", xlabel = "Longitude", ylabel = "Latitude")
    heatmap!(ax1, real_species_matrix, colormap = :inferno)
    ax1.yreversed = true
    
    ax2 = Axis(fig[1, 2], title = "Simulated Distribution of $species_name", xlabel = "Longitude", ylabel = "Latitude")
    heatmap!(ax2, simulated_species_matrix, colormap = :inferno)
    ax2.yreversed = true

    return fig
end

combined_species_distribution_plot(p, "Sus scrofa")

using CairoMakie
using PDFmerger
# Function to save species plots to a PDF
function save_species_plots_to_pdf(array_output, species_names, filename)
    # Create a temporary directory to store the individual PDFs
    temp_dir = mktempdir()

    # Save each plot as a temporary PDF
    for species_name in species_names
        fig = combined_species_distribution_plot(array_output, species_name)
        temp_pdf_path = joinpath(temp_dir, "$species_name.pdf")
        save(temp_pdf_path, fig)
        append_pdf!(filename, temp_pdf_path, cleanup=false)
    end

    # Cleanup temporary directory
    rm(temp_dir, recursive=true)
end

# Example usage
save_species_plots_to_pdf(s, spain_names, "species_distribution.pdf")

###################### PLOTTING SIMULATION RICHNESS OUTPUTS ####################
using CairoMakie
using Serialization
using PDFmerger

function load_and_plot_results(sigmas, epsilons, alphas, output_dir)
    # Define the path for the final merged PDF
    output_pdf = joinpath(output_dir, "result_2000ts_var_epsi_var_sig_var_alpha.pdf")

    # Create a temporary directory to store the individual PDFs
    temp_dir = mktempdir()

    # Counter for missing files
    missing_files_count = 0

    # Loop over each sigma, epsilon, and alpha value
    for sigma in sigmas
        for epsilon in epsilons
            for alpha in alphas
                # Construct the file path for the serialized data
                file_path = joinpath(output_dir, "2000ts_s$(sigma)_e$(epsilon)_a$(alpha).jls")

                # Debugging: Print the path to verify it’s correct
                println("Checking file: $file_path")
                println("Absolute path: $(abspath(file_path))")

                # Check if the file exists
                if isfile(file_path)
                    # Load the serialized data
                    state = deserialize(file_path)

                    # Create the plot using your custom map_plot function with sigma, epsilon, and alpha as the title
                    fig = map_plot(state; type = "image", palette=custom_palette, title="Sigma = $(sigma), Epsilon = $(epsilon), Alpha = $(alpha)", colorrange = (0, 256))

                    # Save each plot as a temporary PDF
                    temp_pdf_path = joinpath(temp_dir, "plot_$(sigma)_$(epsilon)_$(alpha).pdf")
                    save(temp_pdf_path, fig)

                    # Append the individual PDF to the final output PDF
                    append_pdf!(output_pdf, temp_pdf_path, cleanup=false)
                else
                    # Increment the missing files counter
                    missing_files_count += 1
                    println("File not found: $(file_path). Skipping...")
                end
            end
        end
    end

    # Cleanup the temporary directory
    rm(temp_dir, recursive = true)

    # Print a message confirming where the merged PDF was saved
    println("Merged PDF saved as: $output_pdf")
    println("Number of files not found: $missing_files_count")
end


# Example usage
sigmas = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9]
epsilons = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
alpha = [0.1, 0.3, 0.5, 0.7, 0.9]
output_directory = "results"  # Directory containing the serialized files

# Call the function to generate the PDF with all the plots
load_and_plot_results(sigmas, epsilons, alpha, output_directory)

########## PLOTTING theoretical OUTPUTS ########
using CairoMakie
using Serialization
using PDFmerger

file_path = joinpath(output_dir, "matrix_hp_0.9_a0.7_c0.9_s0.8.jls")

function load_and_plot_results(sigmas, alpha_values, connectances, herbivore_proportions, output_dir, sample_size=500)
    # Define the path for the final merged PDF
    output_pdf = joinpath(output_dir, "result_plots_subsample.pdf")

    # Create a temporary directory to store the individual PDFs
    temp_dir = mktempdir()

    # Generate all combinations
    combinations = [(sigma, alpha, connectance, hp) for sigma in sigmas for alpha in alpha_values for connectance in connectances for hp in herbivore_proportions]

    # Randomly select a subsample of 500 combinations
    sampled_combinations = sample(combinations, min(sample_size, length(combinations)))

    # Counter for missing or corrupted files
    skipped_files_count = 0

    # Loop over the sampled combinations
    for (sigma, alpha, connectance, herbivore_proportion) in sampled_combinations
        # Construct the file path for the serialized data
        file_path = joinpath(output_dir, "matrix_hp_$(herbivore_proportion)_a$(alpha)_c$(connectance)_s$(sigma).jls")

        # Check if the file exists
        if isfile(file_path)
            try
                # Attempt to load the serialized data
                state_matrix = deserialize(file_path)

                # Create the plot using a heatmap
                fig = Figure()
                ax = Axis(fig[1, 1], title = "Sigma = $(sigma), Alpha = $(alpha), Connectance = $(connectance), HP = $(herbivore_proportion)")
                heatmap!(ax, state_matrix)

                # Save each plot as a temporary PDF
                temp_pdf_path = joinpath(temp_dir, "plot_s$(sigma)_a$(alpha)_c$(connectance)_hp$(herbivore_proportion).pdf")
                save(temp_pdf_path, fig)

                # Append the individual PDF to the final output PDF
                append_pdf!(output_pdf, temp_pdf_path, cleanup=false)
            catch e
                # If there is an issue with the file, print a warning and skip it
                println("Error loading or processing file $(file_path): $(e). Skipping...")
                skipped_files_count += 1
            end
        else
            skipped_files_count += 1
            println("File not found: $(file_path). Skipping...")
        end
    end

    # Cleanup the temporary directory
    rm(temp_dir, recursive = true)

    # Print a message confirming where the merged PDF was saved
    println("Merged PDF saved as: $output_pdf")
    println("Number of files skipped (missing/corrupted): $skipped_files_count")
end

# Define your variables
sigmas = [0.0001, 0.001, 0.005, 0.008, 0.01, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
alpha_values = [0.01, 0.1, 0.3, 0.7, 0.9, 1.2, 1.5, 1.8]
connectances = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
herbivore_proportions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# Specify the output directory
output_dir = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio\\theoretical outputs"

# Run the function with a subsample of 500 combinations
load_and_plot_results(sigmas, alpha_values, connectances, herbivore_proportions, output_dir, 500)



