using Evolutionary

# Define the objective function
function objective(x)
    return sum((x .- 2).^2)  # Minimize distance to [2,2,2,...]
end

# Configure the Genetic Algorithm (GA) parameters
ga_algorithm = GA(
    populationSize = 100,      # Number of individuals in the population
    selection = tournament(3), # Tournament selection with size 3
    mutationRate = 0.2,        # Mutation probability
    crossoverRate = 0.9        # Crossover probability
)

# Set the optimization options
options = Evolutionary.Options(
    iterations = 100,  # Number of generations
    show_trace = true  # Display progress
)

# Perform the optimization
result = Evolutionary.optimize(
    objective,  # The function to minimize
    rand(5),    # Initial solution (5-dimensional vector)
    ga_algorithm,
    options
)

# Display the results
println("Best Solution Found: ", Evolutionary.minimizer(result))
println("Best Fitness Value: ", Evolutionary.minimum(result))
