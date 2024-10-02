# Function to get combined abundances
function get_combined_abundances(array_output, position; modified = false, caca = false)
    if !modified && !caca
        # Merge birmmals and herps
        combined_abundances = (deepcopy(array_output[end].herps) + deepcopy(array_output[end].birmmals)) .* lambda_DA[position]
    elseif modified && !caca
        combined_abundances = array_output[end].state .* lambda_DA[position]
    elseif caca
        combined_abundances = array_output .* lambda_DA[position]
    else
        error("Invalid combination of parameters")
    end
    return combined_abundances
end

combined_abundances = get_combined_abundances(array_output, position; modified = false, caca = false)

# Create a list of species presence/absence vectors
species_presence = []

for cell in idx
    abundances = combined_abundances[cell].a
    if !any(isnan, abundances)
        presence = abundances .> body_mass_vector
        push!(species_presence, presence)
    else
        push!(species_presence, nothing)
    end
end

# Filter out cells with missing data
valid_indices = findall(x -> x != nothing, species_presence)
valid_cells = [idx[i] for i in valid_indices]
species_presence = species_presence[valid_indices]

