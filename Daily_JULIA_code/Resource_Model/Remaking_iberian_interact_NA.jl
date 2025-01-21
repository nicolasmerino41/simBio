temp_names = names(iberian_interact_NA, 1)
filtered_names = filter(name -> name in birmmals_names, temp_names)

println(predator_prey_dict["Coracias garrulus"])

zero_rows = []
for i in axes(iberian_interact_NA, 1)
    if all(x -> x == 0, iberian_interact_NA[i, :])
        push!(zero_rows, names(iberian_interact_NA, 1)[i])
    end
end

setdiff(spain_names, zero_rows)
println()