# Helper function to get the neighbors
function get_neighbors(matrix, row, col)
    neighbors = []
    rows, cols = size(matrix)
    
    for r in max(row-1, 1):min(row+1, rows)
        for c in max(col-1, 1):min(col+1, cols)
            if (r != row || c != col) && !isnan(matrix[r, c])
                push!(neighbors, matrix[r, c])
            end
        end
    end
    return neighbors
end

matrix_abundances = Matrix(DA_with_abundances)
x = 0
y = 0
for index in CartesianIndices(matrix_abundances)
    if isnan(prec_DA[index]) && !iszero(matrix_abundances[index])
        x += 1
        neighbors = get_neighbors(prec_DA, index[1], index[2])
        if !isempty(neighbors)
            prec_DA[index] = mean(neighbors)
        end
    end

    if isnan(temp_DA[index]) && !iszero(matrix_abundances[index])
        y += 1
        neighbors = get_neighbors(temp_DA, index[1], index[2])
        if !isempty(neighbors)
            temp_DA[index] = mean(neighbors)
        end
    end

end

println("Number of changes made: $x")
println("Number of changes made: $y")