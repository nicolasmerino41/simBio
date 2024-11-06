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
prob = ODEProblem(ecosystem!, u0, tspan)

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
