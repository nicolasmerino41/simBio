include("../Daily JULIA code/One-click code.jl")

include("../Daily JULIA code/human_footprint.jl")

m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

##### RUN SIMULATION #####
using CSV, DataFrames, CairoMakie, PDFmerger  # Added PDFmerger for appending PDFs

# Precompute DA_with_abundances once
DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)

# Function to save parameters, grid type (k_DA name), and metrics to CSV and append plots to the final PDFs
function run_simulation(sigma, epsilon, alpha, position, output_pdf_heatmap, output_pdf_image)

    # Assign the appropriate grid name based on the position
    if position == 1
        k_DA_name = "multiplicative"
    elseif position == 2
        k_DA_name = "additive"
    elseif position == 3
        k_DA_name = "min"
    elseif position == 4
        k_DA_name = "harmonic"
    elseif position == 5      
        k_DA_name = "geometric"  
    end

    # Get the corresponding k_DA grid from the list
    k_DA = k_DA_list[position]

    # Set up the state for the simulation
    pepe_state = (
        state = Matrix(DA_with_abundances),
        k_DA = Matrix(k_DA),
        npp_DA = Matrix(npp_DA)
    )
    
    outdisp = OutwardsDispersal{:state, :state}(;
        formulation=CustomKernel(alpha),
        distancemethod=AreaToArea(30),
        maskbehavior = Dispersal.CheckMaskEdges()
    )
    
    caca = deepcopy(iberian_interact_NA)
    full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon))

    println("sigma  = ", sigma, " epsilon = ", epsilon, " alpha = ", alpha, " k_DA = ", k_DA_name)

    # Run the simulation
    array_output = ResultOutput(
        pepe_state; tspan = 1:rand([3,4,5,6]),
        mask = Matrix(DA_sum)
    )
    p = sim!(array_output, Ruleset(biotic_rule_k, outdisp; boundary = Reflect()))

    # Step 1: Compute metrics from the last timestep of the simulation (p[end])
    avg_shannon = average_shannon_index(p; modified = true)
    avg_bbp = average_bbp(p; modified = true)
    richness_sim = richness_similarity(p; modified = true)
    # alive_preds = alive_predators(p; modified = true)
    # mean_tl = calculate_mean_tl(p; modified = true)

    # Step 2: Save the parameters, grid type, and metrics in a CSV
    results_row = DataFrame(
        sigma = sigma,
        epsilon = epsilon,
        alpha = alpha,
        k_DA_name = k_DA_name,
        avg_shannon = round(avg_shannon, digits = 2),
        avg_bbp = avg_bbp,
        richness_similarity = richness_sim,
        # alive_predators = alive_preds,
        # mean_tl = mean_tl  # Scalar value
    )
    
    # Append or create the CSV file
    csv_filename = "simulation_results.csv"
    if isfile(csv_filename)
        CSV.write(csv_filename, results_row, append = true)
    else
        CSV.write(csv_filename, results_row)
    end

    # Step 3: Append the plots to the final PDFs

    # Create a temporary directory to store the individual plots
    temp_dir = mktempdir()

    # Final state of the simulation
    final_state = p[end].state

    # Heatmap plot
    fig_heatmap = Figure(resolution = (800, 600))
    ax_heatmap = Axis(fig_heatmap[1, 1], title = "sigma = $(sigma), epsilon = $(epsilon), alpha = $(alpha), k_DA = $(k_DA_name)")
    heatmap!(ax_heatmap, final_state; colormap = custom_palette)
    ax_heatmap.yreversed[] = true
    temp_pdf_path_heatmap = joinpath(temp_dir, "temp_heatmap_$(k_DA_name).pdf")
    save(temp_pdf_path_heatmap, fig_heatmap)
    append_pdf!(output_pdf_heatmap, temp_pdf_path_heatmap, cleanup = false)

    # Image plot
    fig_image = Figure(resolution = (800, 600))
    ax_image = Axis(fig_image[1, 1], title = "sigma = $(sigma), epsilon = $(epsilon), alpha = $(alpha), k_DA = $(k_DA_name)")
    image!(ax_image, final_state; colormap = custom_palette)
    ax_image.yreversed[] = true
    temp_pdf_path_image = joinpath(temp_dir, "temp_image_$(k_DA_name).pdf")
    save(temp_pdf_path_image, fig_image)
    append_pdf!(output_pdf_image, temp_pdf_path_image, cleanup = false)

    # Cleanup temporary directory after appending the PDFs
    rm(temp_dir, recursive = true)
end

# Simulation parameters
sigmas = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.5, 2.0, 2.5, 3.0]
epsilons = [0.1, 0.5, 1.0, 2.0, 3.0]
alpha_values = [0.01, 0.05, 0.1, 0.3, 0.6, 0.9]
k_DA_list = [k_DA.DA_multiplicative, k_DA.DA_additive, k_DA.DA_min, k_DA.DA_harmonic, k_DA.DA_geometric]
positions = [1, 2, 3, 4, 5]

# Create the PDF files beforehand
output_pdf_heatmap = "final_results_heatmap.pdf"
output_pdf_image = "final_results_image.pdf"

# Use Threads.@threads to parallelize the loop
Threads.@threads for sigma in sigmas
    for epsilon in epsilons
        for alpha in alpha_values
            for position in positions  
                run_simulation(sigma, epsilon, alpha, position, output_pdf_heatmap, output_pdf_image)
            end
        end
    end
end