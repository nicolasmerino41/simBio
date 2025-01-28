# You need to have k_DA loaded for this step and transform it so all species have a decent K, 
# even predators and even if they're not going to use it
# Preferably if you already have the coordinates of the cell (i, j)

# Initialize the empty DA arrays for each suitability method (same shape as k_DA)
H0_DA = DimArray(reshape([MyBirmmals(SVector{205, Float64}(fill(0.0, 205))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

# Loop through the axes of the DA arrays
for row in axes(H0_DA, 1), col in axes(H0_DA, 2)
    if isone(DA_sum[row, col])
        # Calculate the scaled deviations for bio5, bio6, and bio12
        S_bio5 = 1 ./ (1 .+ abs.(bio5_DA[row, col] .- lax_species_niches.mean_bio5[50:end]) ./ lax_species_niches.sd_bio5[50:end])
        S_bio6 = 1 ./ (1 .+ abs.(bio6_DA[row, col] .- lax_species_niches.mean_bio6[50:end]) ./ lax_species_niches.sd_bio6[50:end])
        S_bio12 = 1 ./ (1 .+ abs.(bio12_DA[row, col] .- lax_species_niches.mean_bio12[50:end]) ./ lax_species_niches.sd_bio12[50:end])
        # println("S_bio5: ", S_bio5, "\n")
        # println("S_bio6: ", S_bio6, "\n")
        # println("S_bio12: ", S_bio12, "\n")

        # # 1. Multiplicative Approach (Original)
        multiplicative_suitability = S_bio5 .* S_bio6 .* S_bio12
        H0_DA[row, col] = MyBirmmals(SVector{205, Float64}(multiplicative_suitability))
        # println("multiplicative_suitability: ", multiplicative_suitability, "\n")
    end
end

# DA_BIRMMALS_WITH_PI CORRECTED FOR PREDATORS WITHOUT PREY
DA_birmmals_with_pi_corrected = deepcopy(DA_birmmals_with_pi)

for i in 1:length(idx)
    
    local_i, local_j = idx[i][1], idx[i][2]
    
    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey = check_predator_has_prey(sp_nm)

    if !predator_has_prey[1]
        filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        # @info("In cell $i, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
        
        vectoro = Vector(DA_birmmals_with_pi[local_i, local_j].a)
        
        for j in predator_has_prey[3]
            index = species_dict_predators_in_birmmals[j]
            vectoro[index] = 0.0
        end

        DA_birmmals_with_pi_corrected[local_i, local_j] = MyBirmmals(SVector{205, Float64}(vectoro))
    end
end