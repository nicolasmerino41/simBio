matrix = [
    0 0 0
    0 0 0
    0 0 0
]

links = OrderedDict{}()
value = 0
for i in size(matrix, 1), j in size(matrix, 2)
    if !iszero(matrix[i, j])
        value =+ 1
        links = append!(links[i => j])
    end
end

using EcologicalNetworksDynamics

food_web = FoodWeb(links)

