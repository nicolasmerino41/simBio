using DifferentialEquations
using Plots
########### ONE HERBIVORE ###########
begin
    # Parameter values
    NPP = 1.0        # Net Primary Productivity
    l = 0.1          # Leaching rate
    a = 0.02         # Attack rate
    epsilon = 0.5    # Assimilation efficiency
    d = 0.1          # Herbivore mortality rate

    # Define the ODE function
    function resource_herbivore!(du, u, p, t)
        R, H = u
        du[1] = NPP - l * R - a * R * H                # dR/dt
        du[2] = epsilon * a * R * H - d * H^2          # dH/dt
    end

    # Initial conditions
    R0 = 10.0    # Initial resource biomass
    H0 = 5.0     # Initial herbivore biomass
    u0 = [R0, H0]

    # Time span for the simulation
    tspan = (0.0, 200.0)

    # Define the ODE problem
    prob = ODEProblem(resource_herbivore!, u0, tspan)

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Plot the results
    PL.plot(sol.t, sol[1, :], label="Resource (R)", xlabel="Time", ylabel="Biomass", linewidth=2)
    PL.plot!(sol.t, sol[2, :], label="Herbivore (H)", linewidth=2)
end

########### MULTIPLE HERBIVORES ###########
begin 
     # For generating random numbers

    # Seed for reproducibility
    Random.seed!(1234)

    # Number of herbivore species
    num_herbivores = 3

    # Parameter values
    NPP = 1.0                                      # Net Primary Productivity
    l = 0.1                                        # Leaching rate
    a = rand(num_herbivores) * 0.05 .+ 0.01        # Attack rates (between 0.01 and 0.06)
    epsilon = rand(num_herbivores) * 0.5 .+ 0.25   # Assimilation efficiencies (between 0.25 and 0.75)
    d = rand(num_herbivores) * 0.01 .+ 0.005       # Mortality rates (between 0.005 and 0.015)

    # Print parameter values for reference
    println("Attack rates (a): ", a)
    println("Assimilation efficiencies (epsilon): ", epsilon)
    println("Mortality rates (d): ", d)

    # Define the ODE function
    function resource_multiple_herbivores!(du, u, p, t)
        R = u[1]                           # Resource biomass
        H = u[2:end]                       # Herbivore biomasses (array)
    
        # Resource equation
        du[1] = NPP - l * R - sum(a .* R .* H)
    
        # Herbivore equations
        for i in 1:num_herbivores
            du[i+1] = epsilon[i] * a[i] * R * H[i] - d[i] * H[i]^2
        end
    end

    # Initial conditions
    R0 = 10.0                               # Initial resource biomass
    H0 = rand(num_herbivores) * 5.0 .+ 1.0  # Initial herbivore biomasses (between 1.0 and 6.0)
    u0 = [R0; H0]                           # Combine into a single array

    # Time span for the simulation
    tspan = (0.0, 200.0)

    # Define the ODE problem
    prob = ODEProblem(resource_multiple_herbivores!, u0, tspan)

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Prepare data for plotting
    t = sol.t
    R = sol[1, :]
    H = sol[2:end, :]  # Each row corresponds to a herbivore species

    # Create a plot object
    plt = PL.plot(title="Resource and Herbivore Dynamics", xlabel="Time", ylabel="Biomass")

    # Add resource biomass to the plot
    PL.plot!(plt, t, R, label="Resource (R)", linewidth=2)

    # Add herbivore biomasses to the plot
    for i in 1:num_herbivores
        PL.plot!(plt, t, H[i, :], label="Herbivore H$i", linewidth=2)
    end

    # Customize plot attributes
    PL.xlims!(plt, 0, 200)              # Set x-axis limits
    PL.ylims!(plt, 0, maximum(R) * 1.1) # Set y-axis limits (adjust as needed)
    # PL.legend!(plt, :right)             # Position the legend
    # PL.grid!(plt, true)                 # Enable grid lines

    # Display the plot
    display(plt)
end

########### MULTIPLE HERBIVORES & ONE PREDATOR ###########
begin
    # Seed for reproducibility
    Random.seed!(1234)

    # Number of herbivore species
    num_herbivores = 3

    # Parameter values
    NPP = 1.0                                      # Net Primary Productivity
    l = 0.1                                        # Leaching rate

    # Herbivore parameters
    a = rand(num_herbivores) * 0.05 .+ 0.01        # Attack rates (between 0.01 and 0.06)
    epsilon = rand(num_herbivores) * 0.5 .+ 0.25   # Assimilation efficiencies (between 0.25 and 0.75)
    d = rand(num_herbivores) * 0.01 .+ 0.005       # Mortality rates (between 0.005 and 0.015)

    # Make herbivore one the most efficient
    epsilon[1] = 0.9                               # Set to a high value
    epsilon[2:end] .= epsilon[2:end] .* 0.5        # Reduce efficiency of others

    # Predator parameters
    b = 0.02          # Attack rate of predator on herbivore one
    epsilon_p = 0.5   # Assimilation efficiency of the predator
    m = 0.01          # Mortality rate of the predator

    # Print parameter values for reference
    println("Attack rates (a): ", a)
    println("Assimilation efficiencies (epsilon): ", epsilon)
    println("Mortality rates (d): ", d)
    println("Predator parameters - b: ", b, ", epsilon_p: ", epsilon_p, ", m: ", m)

    # Define the ODE function
    function resource_predator_herbivores!(du, u, p, t)
        R = u[1]                             # Resource biomass
        H = u[2:num_herbivores+1]            # Herbivore biomasses (array)
        P = u[end]                           # Predator biomass

        # Resource equation
        du[1] = NPP - l * R - sum(a .* R .* H)

        # Herbivore equations
        for i in 1:num_herbivores
            if i == 1
                # Herbivore one includes predation
                du[i+1] = epsilon[i] * a[i] * R * H[i] - d[i] * H[i]^2 - b * H[i] * P
            else
                # Other herbivores
                du[i+1] = epsilon[i] * a[i] * R * H[i] - d[i] * H[i]^2
            end
        end

        # Predator equation
        du[end] = epsilon_p * b * H[1] * P - m * P
    end

    # Initial conditions
    R0 = 10.0                               # Initial resource biomass
    H0 = rand(num_herbivores) * 5.0 .+ 1.0  # Initial herbivore biomasses (between 1.0 and 6.0)
    H0[1] = 5.0                             # Set initial biomass of herbivore one
    u0 = [R0; H0; 1.0]                      # Include predator initial biomass

    # Time span for the simulation
    tspan = (0.0, 200.0)

    # Define the ODE problem
    prob = ODEProblem(resource_predator_herbivores!, u0, tspan)

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Prepare data for plotting
    t = sol.t
    R = sol[1, :]
    H = sol[2:num_herbivores+1, :]  # Each row corresponds to a herbivore species
    P = sol[end, :]                 # Predator biomass

    # Calculate total herbivore biomass
    total_H = sum(H, dims=1)  # Sum across herbivore species for each time point

    # Create a plot object
    plt = PL.plot(title="Resource, Herbivores, and Predator Dynamics", xlabel="Time", ylabel="Biomass")

    # Add resource biomass to the plot
    PL.plot!(plt, t, R, label="Resource (R)", linewidth=2)

    # Add herbivore biomasses to the plot
    for i in 1:num_herbivores
        PL.plot!(plt, t, H[i, :], label="Herbivore H$i", linewidth=2)
    end

    # Add total herbivore biomass in red
    PL.plot!(plt, t, total_H[:], label="Total Herbivores", linewidth=2, color=:red)

    # Add predator biomass to the plot
    PL.plot!(plt, t, P, label="Predator (P)", linewidth=2, linestyle=:dashdot)

    # Customize plot attributes
    PL.xlims!(plt, 0, 200)                            # Set x-axis limits
    PL.ylims!(plt, 0, maximum([maximum(R), maximum(total_H), maximum(P)]) * 1.1)  # Set y-axis limits
    # PL.legend!(plt, :right)                           # Position the legend
    # PL.grid!(plt, true)                               # Enable grid lines

    # Display the plot
    display(plt)
end

########### FLEXIBLE FRAMEWORK ###########
begin
    # Number of herbivore species
    num_herbivores = 3

    # Number of predator species
    num_predators = 2  # You can increase this number as needed

    # Resource parameters
    params_resource = Dict(
        :NPP => 1.0,  # Net Primary Productivity
        :l => 0.1     # Leaching rate
    )

    # Herbivore parameters (arrays)
    params_herbivores = []

    for i in 1:num_herbivores
        push!(params_herbivores, Dict(
            :a => rand() * 0.05 + 0.01,       # Attack rate
            :epsilon => rand() * 0.5 + 0.25,  # Assimilation efficiency
            :d => rand() * 0.01 + 0.005       # Mortality rate
        ))
    end

    # Make herbivore one the most efficient
    params_herbivores[1][:epsilon] = 0.9  # Set to a high value
    for i in 2:num_herbivores
        params_herbivores[i][:epsilon] *= 0.5  # Reduce efficiency of others
    end
    # Predator parameters (array)
    params_predators = []

    for i in 1:num_predators
        push!(params_predators, Dict(
            :b => 0.02,          # Attack rate on herbivore one
            :epsilon_p => 0.5,   # Assimilation efficiency
            :m => 0.01           # Mortality rate
        ))
    end

    # Print parameter values for reference
    println("Resource parameters: ", params_resource)
    println("Herbivore parameters:")
    for i in 1:num_herbivores
        println("  H$i: ", params_herbivores[i])
    end
    println("Predator parameters:")
    for j in 1:num_predators
        println("  P$j: ", params_predators[j])
    end

    # Initial resource biomass
    R0 = 10.0

    # Initial herbivore biomasses (array)
    H0 = [rand() * 5.0 + 1.0 for _ in 1:num_herbivores]
    H0[1] = 5.0  # Set initial biomass of herbivore one

    # Initial predator biomasses (array)
    P0 = [1.0 for _ in 1:num_predators]

    # Combine initial conditions
    u0 = [R0; H0; P0]

    # Time span for the simulation
    tspan = (0.0, 200.0)

    # Combine all parameters into a single dictionary
    params = Dict(
        :num_herbivores => num_herbivores,
        :num_predators => num_predators,
        :params_resource => params_resource,
        :params_herbivores => params_herbivores,
        :params_predators => params_predators
    )

    # Define the ODE function with parameters captured via closure
    function create_ecosystem!(params)

        num_herbivores = params[:num_herbivores]
        num_predators = params[:num_predators]
        params_resource = params[:params_resource]
        params_herbivores = params[:params_herbivores]
        params_predators = params[:params_predators]

        function ecosystem!(du, u, t)
            # Unpack state variables
            R = u[1]                                         # Resource biomass
            H = u[2:1+num_herbivores]                        # Herbivore biomasses
            P = u[2+num_herbivores:end]                      # Predator biomasses

            # Initialize derivatives
            du .= 0.0

            # Resource equation
            du[1] = params_resource[:NPP] - params_resource[:l] * R - sum([params_herbivores[i][:a] * R * H[i] for i in 1:num_herbivores])

            # Herbivore equations
            for i in 1:num_herbivores
                predation = 0.0
                # If there are predators, calculate predation on herbivore one
                if num_predators > 0 && i == 1
                    for j in 1:num_predators
                        predation += params_predators[j][:b] * H[i] * P[j]
                    end
                end
                du[1 + i] = params_herbivores[i][:epsilon] * params_herbivores[i][:a] * R * H[i] - params_herbivores[i][:d] * H[i]^2 - predation
            end

            # Predator equations
            for j in 1:num_predators
                # Predator preys only on herbivore one
                du[1 + num_herbivores + j] = params_predators[j][:epsilon_p] * params_predators[j][:b] * H[1] * P[j] - params_predators[j][:m] * P[j]
            end
        end

        return ecosystem!
    end

    # Create the ODE function with parameters captured via closure
    ecosystem!! = create_ecosystem!(params)

    # Define the ODE problem without passing parameters
    prob = ODEProblem(ecosystem!!, u0, tspan)

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Prepare data for plotting
    t = sol.t
    R = sol[1, :]
    H = sol[2:1+num_herbivores, :]   # Each row corresponds to a herbivore species
    P = sol[2+num_herbivores:end, :] # Each row corresponds to a predator species

    # Calculate total herbivore biomass
    total_H = sum(H, dims=1)  # Sum across herbivore species for each time point

    # Create a plot object
    plt = plot(title="Ecosystem Dynamics", xlabel="Time", ylabel="Biomass")

    # Add resource biomass to the plot
    plot!(plt, t, R, label="Resource (R)", linewidth=2)

    # Add herbivore biomasses to the plot
    for i in 1:num_herbivores
        plot!(plt, t, H[i, :], label="Herbivore H$i", linewidth=2)
    end

    # Add total herbivore biomass in red
    plot!(plt, t, total_H[:], label="Total Herbivores", linewidth=2, color=:red)

    # Add predator biomasses to the plot
    for j in 1:num_predators
        plot!(plt, t, P[j, :], label="Predator P$j", linewidth=2, linestyle=:dashdot)
    end

    # Customize plot attributes
    xlims!(plt, 0, tspan[2])  # Set x-axis limits
    ylims!(plt, 0, maximum([maximum(R), maximum(total_H), maximum(P)]) * 1.1)  # Set y-axis limits
    legend!(plt, :right)      # Position the legend
    grid!(plt, true)          # Enable grid lines

    # Display the plot
    display(plt)
end

########### NEW RESOURCE EXPLICIT, NPP-BASED MODEL ###########
begin
    # Seed for reproducibility
    Random.seed!(rand(1:100, 1)[1])

    # Number of species
    num_producers = 10   # Number of producer species
    num_herbivores = 5  # Number of herbivore species

    # Resource parameters
    NPP = 1.0    # Net Primary Productivity
    l = 0.1      # Leaching rate

    # Producer parameters
    params_producers = []
    for alpha in 1:num_producers
        push!(params_producers, Dict(
            :c => rand() * 0.05 + 0.01,     # Uptake rate
            :d => rand() * 0.01 + 0.005     # Mortality rate
        ))
    end

    # Herbivore parameters
    params_herbivores = []
    for i in 1:num_herbivores
        push!(params_herbivores, Dict(
            :d => rand() * 0.01 + 0.005     # Mortality rate
        ))
    end

    # Mortality rates for producers and herbivores
    d_p = [params_producers[alpha][:d] for alpha in 1:num_producers]
    d_h = [params_herbivores[i][:d] for i in 1:num_herbivores]

    # Interaction parameters
    epsilon = 0.5  # Assimilation efficiency (assumed same for all)

    # Attack rates (a_i_alpha)
    a = [rand() * 0.05 + 0.01 for i in 1:num_herbivores, alpha in 1:num_producers]

    # Initial resource biomass
    R0 = 10.0

    # Initial producer biomasses
    P0 = [rand() * 5.0 + 1.0 for _ in 1:num_producers]

    # Initial herbivore biomasses
    H0 = [rand() * 2.0 + 1.0 for _ in 1:num_herbivores]

    # Combine all initial conditions
    u0 = [R0; P0; H0]

    # Time span for the simulation
    tspan = (0.0, 200.0)

    # Define the ODE function
    function ecosystem!(du, u, p, t)
        # Parameters are captured via closure
        # Unpack state variables
        R = u[1]
        P = u[2:1+num_producers]
        H = u[2+num_producers:end]

        # Initialize derivatives
        du .= 0.0

        # Resource dynamics
        consumption = sum([params_producers[alpha][:c] * P[alpha] * R for alpha in 1:num_producers])
        du[1] = NPP - l * R - consumption

        # Producer dynamics
        for alpha in 1:num_producers
            # Total herbivory on producer alpha
            herbivory = sum([a[i, alpha] * H[i] for i in 1:num_herbivores])
            # Corrected producer dynamics with density-dependent mortality
            du[1 + alpha] = epsilon * params_producers[alpha][:c] * R * P[alpha] - d_p[alpha] * P[alpha]^2 - P[alpha] * herbivory
        end

        # Herbivore dynamics
        for i in 1:num_herbivores
            # Total consumption of producers by herbivore i
            consumption_i = sum([a[i, alpha] * P[alpha] for alpha in 1:num_producers])
            # Corrected herbivore dynamics with density-dependent mortality
            du[1 + num_producers + i] = epsilon * H[i] * consumption_i - d_h[i] * H[i]^2
        end
    end

    # Define the ODE problem without passing parameters
    prob = ODEProblem(ecosystem!, u0, tspan)

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Prepare data for plotting
    t = sol.t
    R = sol[1, :]
    P = sol[2:1+num_producers, :]
    H = sol[2+num_producers:end, :]

    # Create a plot object
    plt = Plots.plot(title="Ecosystem Dynamics with Density-Dependent Mortality", xlabel="Time", ylabel="Biomass")

    # Plot resource biomass
    Plots.plot!(plt, t, R, label="Resource (R)", linewidth=2)

    # Plot producer biomasses
    for alpha in 1:num_producers
        Plots.plot!(plt, t, P[alpha, :], label="Producer P$alpha", linewidth=2)
    end

    # Plot herbivore biomasses
    for i in 1:num_herbivores
        Plots.plot!(plt, t, H[i, :], label="Herbivore H$i", linewidth=2, linestyle=:dash)
    end  

    # Customize plot attributes
    Plots.xlims!(plt, 0, tspan[2])
    Plots.ylims!(plt, 0, maximum([maximum(R), maximum(P), maximum(H)]) * 1.1)

    # Display the plot
    display(plt)
end

########### NEW RESOURCE EXPLICIT, NPP-BASED MODEL WITH PREDATORS ###########
begin
    # Seed for reproducibility
    Random.seed!(1234)

    # Verbose flag
    verbose = true  # Set to `false` to disable printing extinction messages

    # Number of species
    num_producers = 200    # Number of producer species
    num_herbivores = 50   # Number of herbivore species
    num_predators = 20    # Number of predator species

    # Resource parameters
    NPP = 1.0    # Net Primary Productivity
    l = 0.1      # Leaching rate

    # Producer parameters
    params_producers = []
    for alpha in 1:num_producers
        push!(params_producers, Dict(
            :c => rand() * 0.05 + 0.01,     # Uptake rate
            :d => rand() * 0.01 + 0.005     # Mortality rate
        ))
    end

    # Herbivore parameters
    params_herbivores = []
    for i in 1:num_herbivores
        push!(params_herbivores, Dict(
            :d => rand() * 0.01 + 0.005     # Mortality rate
        ))
    end

    # Predator parameters
    params_predators = []
    for k in 1:num_predators
        push!(params_predators, Dict(
            :d => rand() * 0.01 + 0.005     # Mortality rate
        ))
    end

    # Mortality rates
    d_p = [params_producers[alpha][:d] for alpha in 1:num_producers]
    d_h = [params_herbivores[i][:d] for i in 1:num_herbivores]
    d_c = [params_predators[k][:d] for k in 1:num_predators]

    # Assimilation efficiencies
    epsilon_P = 0.5  # Producers
    epsilon_H = 0.5  # Herbivores
    epsilon_C = 0.5  # Predators

    # Attack rates (trophic matrices)
    a = [rand() * 0.05 + 0.01 for i in 1:num_herbivores, alpha in 1:num_producers]   # Herbivore on producer
    b = [rand() * 0.05 + 0.01 for k in 1:num_predators, i in 1:num_herbivores]       # Predator on herbivore

    # Initial conditions
    R0 = 10.0
    P0 = [rand() * 5.0 + 1.0 for _ in 1:num_producers]
    H0 = [rand() * 2.0 + 1.0 for _ in 1:num_herbivores]
    C0 = [rand() * 1.0 + 0.5 for _ in 1:num_predators]
    u0 = [R0; P0; H0; C0]

    # Time span
    tspan = (0.0, 200.0)

    # Species names for printing (optional)
    producer_names = ["Producer P$α" for α in 1:num_producers]
    herbivore_names = ["Herbivore H$i" for i in 1:num_herbivores]
    predator_names = ["Predator C$k" for k in 1:num_predators]
    species_names = vcat(producer_names, herbivore_names, predator_names)  # Concatenate species names

    # Define the ODE function
    function ecosystem!(du, u, p, t)
        # Unpack state variables
        R = u[1]
        P = u[2:1+num_producers]
        H = u[2+num_producers : 1+num_producers+num_herbivores]
        C = u[2+num_producers+num_herbivores : end]

        # Initialize derivatives
        du .= 0.0

        ## Resource dynamics
        resource_consumption = sum([params_producers[alpha][:c] * P[alpha] * R for alpha in 1:num_producers])
        du[1] = NPP - l * R - resource_consumption

        ## Producer dynamics
        for alpha in 1:num_producers
            # Total herbivory on producer alpha
            herbivory = sum([a[i, alpha] * H[i] for i in 1:num_herbivores])
            # Producer dynamics
            du[1 + alpha] = epsilon_P * params_producers[alpha][:c] * R * P[alpha] - d_p[alpha] * P[alpha]^2 - P[alpha] * herbivory
        end

        ## Herbivore dynamics
        for i in 1:num_herbivores
            # Total consumption of producers by herbivore i
            consumption_i = sum([a[i, alpha] * P[alpha] for alpha in 1:num_producers])
            # Predation on herbivore i by predators
            predation = sum([b[k, i] * C[k] for k in 1:num_predators])
            # Herbivore dynamics
            du[1 + num_producers + i] = epsilon_H * H[i] * consumption_i - d_h[i] * H[i]^2 - H[i] * predation
        end

        ## Predator dynamics
        for k in 1:num_predators
            # Total consumption of herbivores by predator k
            consumption_k = sum([b[k, i] * H[i] for i in 1:num_herbivores])
            # Predator dynamics
            du[1 + num_producers + num_herbivores + k] = epsilon_C * C[k] * consumption_k - d_c[k] * C[k]^2
        end
    end

    # Define the condition function for the callback
    function condition_extinction(u, t, integrator)
        min_biomass = minimum(u[2:end])  # Exclude resource
        return min_biomass
    end

    # Define the affect function for the callback with verbose option
    function affect_extinction!(integrator)
        u = integrator.u
        t = integrator.t
        # Check and set negative biomasses to zero
        for idx in 2:length(u)  # Exclude resource
            if u[idx] < 0.0
                u[idx] = 0.0  # Set biomass to zero
            end
        end
        # Print extinction messages if verbose is true
        if verbose
            for idx in 2:length(u)  # Exclude resource
                if integrator.u[idx] == 0.0 && integrator.uprev[idx] > 0.0
                    species_idx = idx - 1  # Adjust index since u[1] is resource
                    if species_idx ≤ num_producers
                        species_name = producer_names[species_idx]
                    elseif species_idx ≤ num_producers + num_herbivores
                        species_name = herbivore_names[species_idx - num_producers]
                    else
                        species_name = predator_names[species_idx - num_producers - num_herbivores]
                    end
                    println("Species $(species_name) went extinct at time t = $(t).")
                end
            end
        end
    end

    # Create the ContinuousCallback
    extinction_callback = ContinuousCallback(condition_extinction, affect_extinction!; save_positions=(true, true))

    # Define the ODE problem
    prob = ODEProblem(ecosystem!, u0, tspan)

    # Solve the ODEs with the callback
    sol = solve(prob, Tsit5(); callback=extinction_callback)

    # Prepare data for plotting
    t = sol.t
    R = sol[1, :]
    P = sol[2:1+num_producers, :]
    H = sol[2+num_producers : 1+num_producers+num_herbivores, :]
    C = sol[2+num_producers+num_herbivores : end, :]

    # Create a plot object
    plt = Plots.plot(title="Ecosystem Dynamics with Predators and Extinction Events", xlabel="Time", ylabel="Biomass")

    # Plot resource biomass
    Plots.plot!(plt, t, R, label=false, linewidth=2)

    # Plot producer biomasses
    for alpha in 1:num_producers
        Plots.plot!(plt, t, P[alpha, :], label=false, linewidth=2)
    end

    # Plot herbivore biomasses
    for i in 1:num_herbivores
        Plots.plot!(plt, t, H[i, :], label=false, linewidth=2, linestyle=:dash)
    end

    # Plot predator biomasses
    for k in 1:num_predators
        Plots.plot!(plt, t, C[k, :], label=false, linewidth=2, linestyle=:dot)
    end

    # Customize plot attributes
    Plots.xlims!(plt, 0, tspan[2])
    Plots.ylims!(plt, 0, maximum([maximum(R), maximum(P), maximum(H), maximum(C)]) * 1.1)

    # Display the plot
    display(plt)
end

