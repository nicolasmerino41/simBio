include("One-click code.jl")

include("human_footprint.jl")

m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

##### RUN SIMULATION #####
# using CairoMakie, PDFmerger  # Added PDFmerger for appending PDFs
# caca = deepcopy(iberian_interact_NA)
# sigma = 1.0
# epsilon = 1.0
# full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon))
# Precompute DA_with_abundances once
DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)

# Function to save parameters, grid type (k_DA name), and metrics to CSV and append plots to the final PDFs
function run_simulation(sigma, epsilon, alfa, position)

    k_DA_name = k_DA_names[position]
    
    # Get the corresponding k_DA grid from the list
    k_DA = k_DA_list[position]

    # Set up the state for the simulation
    pepe_state = (
        state = Matrix(DA_with_abundances),
        k_DA = Matrix(k_DA),
        npp_DA = Matrix(k_DA)
    )
    
    outdisp = OutwardsDispersal{:state, :state}(;
        formulation=CustomKernel(alfa),
        distancemethod=AreaToArea(30),
        maskbehavior = Dispersal.CheckMaskEdges()
    )
    
    caca = deepcopy(iberian_interact_NA)
    full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon))

    function trophic_optimized(abundances, full_IM)
        # Calculate the weighted interaction directly
        interaction = full_IM * abundances.a
        return MyStructs256(SVector{256, Float64}(interaction .* abundances.a))
    end
    biotic_rule_k = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
        # if any(isinf, state.a) || any(isnan, state.a)
        #     @error "state has NA values"
        #     println(I)
        # end
        merged_state = state + 
            int_Gr_for_biotic_k(state, self_regulation, k_DA)  +
            trophic_optimized(state, full_IM)
        return MyStructs256(SVector{256, Float64}(max.(0.0, merged_state.a)))
    end

    println("sigma  = ", sigma, " epsilon = ", epsilon, " alfa = ", alfa, " k_DA = ", k_DA_name)

    # Run the simulation
    array_output = ResultOutput(
        pepe_state; tspan = 1:1000,
        mask = Matrix(DA_sum)
    )
    # println("output done")
    p = sim!(array_output, Ruleset(biotic_rule_k, outdisp; boundary = Reflect()))
    # println("simulation done")
    # Step 1: Compute metrics from the last timestep of the simulation (p[end])
    avg_shannon = average_shannon_index(p, position; modified = true)
    # println("shannon done")
    avg_bbp = average_bbp(p, position; modified = true)
    # println("bbp done")
    richness_sim = richness_similarity(p, position; modified = true)
    # println("richness done")
    alive_preds = alive_predators(p, position; modified = true)
    # println("predators done")
    mean_tl = calculate_mean_tl(p, position; modified = true)
    # println("meantl done")
    
    final_state = p[end].state
    NaNs = any(i -> any(isnan, final_state[idx[i][1], idx[i][2]].a), 1:length(idx)) ? 1.0 : 0.0

    # println(NaNs)
    # Step 2: Save the parameters, grid type, and metrics in a CSV
    results_row = DataFrame(
        sigma = sigma,
        epsilon = epsilon,
        alfa = alfa,
        k_DA_name = k_DA_name,
        avg_shannon = round(avg_shannon, digits = 2),
        avg_bbp = round(avg_bbp, digits = 2),
        richness_similarity = round(richness_sim, digits = 2),
        alive_predators = round(alive_preds, digits = 2),
        mean_tl = round(mean_tl, digits = 2),
        NaNs = NaNs
    )

    # Append or create the CSV file
    csv_filename = "simulation_results_biode_new.csv"
    if isfile(csv_filename)
        CSV.write(csv_filename, results_row, append = true)
    else
        CSV.write(csv_filename, results_row)
    end
end

# Simulation parameters
sigmas = [0.001, 0.005, 0.008, 0.01, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.5]
epsilons = [0.1, 0.5, 1.0, 2.0, 3.0]
alfa_values = [0.01, 0.05, 0.1, 0.3, 0.6, 0.9]
k_DA_list = [k_DA.DA_multiplicative, k_DA.DA_additive, k_DA.DA_min, k_DA.DA_harmonic, k_DA.DA_geometric]
k_DA_names = ["multiplicative", "additive", "min", "harmonic", "geometric"]
positions = [1, 2, 3, 4, 5]

# Use Threads.@threads to parallelize the loop
Threads.@threads for sigma in sigmas
    for epsilon in epsilons
        for alfa in alfa_values
            for position in positions  
                run_simulation(sigma, epsilon, alfa, position)
            end
        end
    end
end

################## GLV THREADING PART #############################
sigmas = [0.0001, 0.001, 0.005, 0.008, 0.01, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
alfa_values = [0.01, 0.1, 0.3, 0.7, 0.9, 1.2, 1.5, 1.8]
connectances = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
herbivore_proportions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# sigmas = [0.001]
# alfa_values = [0.01, 0.1]
# connectances = [0.01, 0.1]
# herbivore_proportions = [0.1, 0.5]

include("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\Controlled GLV experiment.jl")
Threads.@threads for sigma in sigmas
    for herbivore_proportion in herbivore_proportions
        for alfa in alfa_values
            for connectance in connectances  
                A_matrix = Adjacency_Matrix_Maker(num_species, connectance, herbivore_proportion)[1]
                I_matrix = Interaction_Matrix_Maker(A_matrix, sigma)
                I_matrix = fill_diagonal!(I_matrix, self_regulation)
                function GLV(state::MyStructs256, k_DA::MyStructs256)
                     return MyStructs256(
                    SVector{num_species, Float64}(
                   state.a + (state.a .* (k_DA.a - state.a) + ((I_matrix * state.a) .* state.a)) 
                   )
                 )
                end

biotic = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return MyStructs256(SVector{num_species, Float64}(max.(0.0, GLV(state, k_DA).a)))
end

disp = OutwardsDispersal{:state, :state}(
    formulation = CustomKernel(alfa),
    distancemethod = AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges(),
)

pepe = (
    state = Matrix(raster_with_abundances),
    k_DA = Matrix(raster_k_DA),
)

array_output = ResultOutput(
    pepe, tspan = 1:100;
    mask = Matrix(raster_sum),
)

@time a = sim!(array_output, Ruleset(biotic, disp; boundary = Reflect(), proc = ThreadedCPU()))
                # Lock for thread safety when modifying shared resources
                serialize("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio\\theoretical outputs\\matrix_hp_$(herbivore_proportion)_a$(alfa)_c$(connectance)_s$(sigma).jls", a[end].state)
            end
        end
    end
end
