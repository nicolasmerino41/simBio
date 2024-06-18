herbivore_names = CSV.File(joinpath(meta_path, "herbivore_names.csv")) |> DataFrame
herbivore_names = convert(Vector{String}, herbivore_names[:, 2])

competition_matrix = deepcopy(iberian_interact_NA)

for i in names(competition_matrix, 1), j in names(competition_matrix, 2)
    competition_matrix[i, j] = 0.0
end

for i in names(competition_matrix, 1), j in names(competition_matrix, 2)
    # println(i, " ", j, i in herbivore_names, j in herbivore_names)
    if i in herbivore_names && i != j && j in herbivore_names
        competition_matrix[i, j] = 1.0
    end 
end

sigma = 0.1
for i in axes(competition_matrix, 1), j in axes(competition_matrix, 2)
    if competition_matrix[i, j] == 1.0
        normal_dist = Normal(0, sigma)
        x = round(abs(rand(normal_dist)), digits = 20)
        competition_matrix[i, j] = -x 
    end 
end

