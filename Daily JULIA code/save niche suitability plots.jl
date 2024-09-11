using CairoMakie
using PDFmerger  # To merge individual PDFs into a single PDF

# Function to save suitability plots to a single PDF
function save_suitability_plots_to_pdf(k_DA, filename)
    # Create a temporary directory to store individual PDFs
    temp_dir = mktempdir()

    # List of plot types and corresponding titles
    suitability_types = [
        ("Multiplicative Suitability", Matrix(k_DA.DA_multiplicative)),
        ("Additive Suitability", Matrix(k_DA.DA_additive)),
        ("Geometric Suitability", Matrix(k_DA.DA_geometric)),
        ("Minimum Suitability", Matrix(k_DA.DA_min)),
        ("Harmonic Suitability", Matrix(k_DA.DA_harmonic))
    ]

    # First, calculate the maximum value across all matrices for consistent color scaling
    total_max = maximum([maximum(Matrix(k_DA.DA_multiplicative)), maximum(Matrix(k_DA.DA_additive)),
                         maximum(Matrix(k_DA.DA_geometric)), maximum(Matrix(k_DA.DA_min)),
                         maximum(Matrix(k_DA.DA_harmonic))]).b

    # Save each plot as a temporary PDF
    for (title, matrix) in suitability_types
        fig = Figure(resolution = (600, 400))  # Create a new figure for each plot
        ax = Axis(fig[1, 1], title = title)
        MK.heatmap!(matrix; colormap = custom_palette, title = title, colorrange = (0, total_max))
        ax.yreversed[] = true
        temp_pdf_path = joinpath(temp_dir, "$title.pdf")
        save(temp_pdf_path, fig)
        append_pdf!(filename, temp_pdf_path, cleanup=false)  # Append to the main PDF
    end

    # Cleanup temporary directory
    rm(temp_dir, recursive=true)
end

# Example usage
save_suitability_plots_to_pdf(k_DA, "suitability_plots.pdf")