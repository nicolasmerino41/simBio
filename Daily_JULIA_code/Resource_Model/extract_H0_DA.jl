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
