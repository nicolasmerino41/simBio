# Suppose iberian_interact_NA is a NamedMatrix, dimension (num_species x num_species).
# The row/column names might be accessible via:
row_species = names(iberian_interact_NA, 1)  # row names as a Vector{String}
col_species = names(iberian_interact_NA, 2)  # col names as a Vector{String}
spain_names = row_species
# herbivore_csv = CSV.read("DFs/herbivore_names.csv", DataFrame)
# herbivore_names = herbivore_csv.value

# In a well-formed NamedMatrix, row_species == col_species.

# We'll define a dictionary that maps each species name to its integer index.
# The index is for both row and column, since it's presumably symmetrical in names.

species_dict = Dict{String, Int}()

for (idx, sp) in enumerate(row_species)
    species_dict[sp] = idx
end

# SAME but only to identify the herbivores positions in a MyBirmmals struct
birmmals_names = spain_names[50:256]
species_dict_herbivores_in_birmmals = Dict{String, Int}()

for (idx, sp) in enumerate(birmmals_names)
    if sp in herbivore_names
        species_dict_herbivores_in_birmmals[sp] = idx
    end
end
