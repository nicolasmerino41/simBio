#### We needed to loaded this way cause there was column problems but it works fine
ga_cell_results = CSV.read("C:/Users/MM-1/Downloads/ga_cell_results.csv", DataFrame, header=false)

function plotting_ga_results(ga_cell_results)
    grid = deepcopy(float(DA_sum))
    
    for i in axes(grid, 1), j in axes(grid, 2)
        if grid[i, j] == 0
            grid[i, j] = NaN
        elseif grid[i, j] == 1.0
            grid[i, j] = NaN
        end
    end

    for i in 1:nrow(ga_cell_results)
        coord = idx[ga_cell_results[i, 1]]
        local_i, local_j = coord[1], coord[2]
        grid[local_i, local_j] = ga_cell_results[i, 3]
        
    end

    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1])
    MK.heatmap!(
        ax, grid, interpolate = false, colormap = custom_palette,
        colorrange = (0.2, 0.8)
    )
    ax.yreversed[] = true
    display(fig)
    
end

plotting_ga_results(ga_cell_results)