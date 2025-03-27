# Assume A_j is the Jacobian and A_p is the parameter tuple containing S and R,
# and that sp_nm is the ordered species list: first S entries for herbivores (and omnivores), then predators.
# Also assume that omnivore_names and predator_names are defined.

S = A_p.S
R = A_p.R
n = S + R

# Compute eigenvalues and eigenvectors
ev = eigen(A_j)
vals = ev.values
vecs = ev.vectors

# Create a vector to store guild labels for each eigenvalue
guild_labels = Vector{String}(undef, n)
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[idx[2][1], idx[2][2]])
# For each eigenvector, find the index of the largest magnitude component.
for i in 1:n
    v = vecs[:, i]
    idx_max = argmax(abs.(v))
    if idx_max <= S
        # For indices in the herbivore block, check if the species is an omnivore.
        if sp_nm[idx_max] in omnivore_names
            guild_labels[i] = "Omnivore"
        else
            guild_labels[i] = "Herbivore"
        end
    else
        guild_labels[i] = "Predator"
    end
end

# Map guild labels to colors.
colors = map(g -> g == "Herbivore" ? :blue : g == "Omnivore" ? :green : :red, guild_labels)

# Extract real and imaginary parts of eigenvalues.
re_vals = real.(vals)
println(sort(re_vals))
im_vals = imag.(vals)

# Plot eigenvalues in the complex plane.
Plots.scatter(re_vals, im_vals, xlabel="Real Part", ylabel="Imaginary Part",
        title="Eigenvalue Plot by Guild", color=colors, marker=:circle)
hline!([0], linestyle=:dash, color=:black)
vline!([0], linestyle=:dash, color=:black)

J = A_j
# Suppose J is the Jacobian matrix computed at equilibrium (for example, from ForwardDiff.jacobian)
# and t_val is the time at which you want to compute the sensitivity matrix.
t_val = 10.0  # you can choose any time point of interest

# Compute the sensitivity matrix S(t) as the matrix exponential of J*t.
S_t = exp(J * t_val)

println("Sensitivity matrix S(t) = exp(J*t):")
println(S_t)

# Optionally, if you want to visualize the sensitivity matrix as a heatmap:
Plots.heatmap(S_t, title = "Sensitivity Matrix at t = $t_val", xlabel = "Initial State Index", ylabel = "State Index at t")

