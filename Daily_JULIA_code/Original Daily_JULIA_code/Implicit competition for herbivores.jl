function turn_comp_into_inter(adjacencyy, sigma, assymetry, self_regulation = 0.0)
    adjacency = deepcopy(adjacencyy)
    u = adjacency  # Initialize the interaction matrix
    
    # Loop through the matrix
    for i in names(adjacency, 1)
        for j in names(adjacency, 2)
            # Ensure we only process species pairs where there is competition (non-zero entries)
            if i != j && adjacency[i, j] != 0.0 && adjacency[j, i] != 0.0
                # Skip if the reverse interaction is already calculated
                if u[j, i] != 0.0 && u[j, i] != 1.0
                    continue
                end

                # Draw interaction strength from normal distribution
                normal_dist = Normal(0, sigma)
                x = round(abs(rand(normal_dist)), digits = 20)  # Interaction for i, j

                # For j, i: adjust based on assymetry (correlation factor)
                y = assymetry * x + (1 - assymetry) * round(abs(rand(normal_dist)), digits = 10)

                # Assign interactions with asymmetry
                u[i, j] = - x  # Species i competes with j
                u[j, i] = - y  # Species j competes with i (asymmetrically based on assymetry)

            elseif i == j
                u[i, j] = -self_regulation  # Self-regulation on diagonal
            end
        end
    end
    return u
end