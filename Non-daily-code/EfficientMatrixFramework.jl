# Load a DataFrame from a serialized file ('.jls' format).
iberian_interact_df = deserialize("iberian_interact_df.jls")
# Convert the DataFrame to a matrix for easier manipulation.
iberian_interact_matrix = iberian_interact_df |> Matrix
# Convert the modified matrix back to a DataFrame, preserving the original column names.
iberian_interact_df = DataFrame(iberian_interact_matrix, names(iberian_interact_df))
# Create a NamedArray from the matrix, using the DataFrame's column names for both dimensions.
iberian_interact_NA = NamedArray(
    iberian_interact_matrix, 
    (names(iberian_interact_df), names(iberian_interact_df)),
    ("Species", "Species")
)

random_iberian_interact = deepcopy(iberian_interact_NA)
for row in axes(iberian_interact_NA, 1), col in axes(iberian_interact_NA, 2)
    if iberian_interact_NA[row, col] != 0.0
        random_iberian_interact[col, row] = -1.0
    end
end
# new = turn_adj_into_inter(iberian_interact_NA, 100, 0.5)

# Initialize an empty OrderedDict to hold the resulting matrices
results = OrderedDict{Float64, OrderedDict{Float64, Matrix}}()

# Iterate over epsilon values
@time for epsilon in [0.1, 0.5, 1.0, 1.5, 3.0]
    # Initialize an OrderedDict for the current epsilon
    epsilon_results = OrderedDict{Float64, Matrix}()
    
    # Iterate over sigma values from 1.0 to 0.01, decrementing by 0.001 each time
    for sigma in 0.001
        caca = deepcopy(iberian_interact_NA)
        
        # Call the turn_adj_into_inter function with the current sigma and epsilon values
        result_matrix = turn_adj_into_inter(caca, sigma, epsilon)
        
        # Append the result to the epsilon_results OrderedDict
        epsilon_results[sigma] = result_matrix
    end
    
    # Store the epsilon_results OrderedDict in the main results OrderedDict
    results[epsilon] = epsilon_results
end

# for i in 1:length(results)
# mean_v = []
# for row in axes(results[i], 1), col in axes(results[i], 2)
#     if random_iberian_interact[row, col] != 0.0
#         push!(mean_v, abs(results[i][row, col]))
#     end
# end
# println(mean(mean_v))
# end

# # Loop for everything
# sim_results = OrderedDict{Float64, OrderedDict{Float64, Matrix}}()
# @time for epsilon in [0.1, 0.5, 1.0, 1.5, 3.0]
#     # Initialize an OrderedDict for the current epsilon
#     epsilon_results = OrderedDict{Float64, ArrayOutput}()

#     for sigma in 1.0:-0.01:0.1
#         full_IM = results[epsilon][sigma]
        
#         array_output = ArrayOutput(
#             DA_with_abundances_p; tspan = 1:100,
#             mask = DA_sum_p
#         )
        
#         r = sim!(array_output, ruleset)

#         # Append the result to the epsilon_results OrderedDict
#         epsilon_results[sigma] = r
#         println("epsilon", epsilon, " sigma: ", sigma)
#     end
#     # Store the epsilon_results OrderedDict in the main results OrderedDict
#     sim_results[epsilon] = epsilon_results
# end

# for i in 0.1:-0.001:0.01
#     println(i)
    
# end

# r = turn_adj_into_inter(iberian_interact_NA, 0.001, 0.5)