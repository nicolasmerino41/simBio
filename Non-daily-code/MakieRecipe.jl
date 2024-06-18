using Pkg
# Desktop PC
# Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

# meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
################## MYSTRUCTS256 MAKIE RECIPE ######################
###################################################################
###################################################################
@MK.recipe(MyName, abundance) do scene
    Attributes(
        abundance_color = :inferno
    )
end

function Makie.plot!(plot::MyName)
    # Extract the abundance matrix from the plot object
    abundance = plot[1][]  # Use `[]` to get the value from the Observable
    
    println("Type of abundance: ", typeof(abundance))
    # Extract the `b` elements from the DimArray
    b_matrix = [mystruct.b for mystruct in abundance]
        
    # Reshape the `b_matrix` to match the original dimensions
    reshaped_b_matrix = reshape(b_matrix, size(abundance))

    # Create the heatmap plot
    MK.heatmap!(plot, reshaped_b_matrix, colormap = plot[:abundance_color])
    
    println("Type of plot: ", typeof(Observable(plot)))
    # Return the updated plot
    return Observable(plot)
    
end



# example_matrix = [MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0)) for _ in 1:3, _ in 1:3]
# b_matrix = [mystruct.b for mystruct in example_matrix]
# MK.heatmap(b_matrix, colormap = :viridis)

# Create an example matrix of `MyStructs256`
example_matrix = [MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0)) for _ in 1:3, _ in 1:3]
example_matrix[2, 2] = MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0))
observable_matrix = Observable(example_matrix)
# example_matrix = DimArray(example_matrix, (Dim{:a}(1:3), Dim{:b}(1:3)))
# Plot the matrix using the custom recipe
fig, ax, pl = myname(example_matrix, axis = (aspect = DataAspect(),))
#= NOTE FOR NEXT TIME YOU COME BACK: It works up to here but then it looks like myname!()
it's not able to update the observable of the frame properly, I have no idea on how to fix it=#
##################################MakieOutput###################################
init_simple = (a = example_matrix,
        b = example_matrix,
        c = example_matrix,
        d = example_matrix
)
makie_output = MakieOutput(init_simple, tspan = 1:100; fps = 50, ruleset = Ruleset(merge_rule)) do (; layout, frame)
        ax1 = Axis(layout[1, 1])
        ax2 = Axis(layout[1, 2])
        ax3 = Axis(layout[2, 1])
        ax4 = Axis(layout[2, 2])
        # hideydecorations!(ax)
        # hidexdecorations!(ax)
        # Plot the frame data
        println(typeof(frame))
        myname!(ax1, frame.a; interpolate=false, colormap=:inferno)
        myname!(ax2, frame.b; interpolate=false)
        myname!(ax3, frame.c; interpolate=false)
        myname!(ax4, frame.d; interpolate=false)
end

cell_m = Cell{:a, :a}() do data, a, I
    return a *rand(10:100)
end

cell_2 = Cell{:b, :b}() do data, a, I
    return a*1000
end

output = MakieOutput(init_simple; tspan=1:10, ruleset=Ruleset(rule), fps=5) do (; layout, frame)

    # Define the global color limits
    color_limits = (10.0, 2000.0)
    
    # Setup the keys and titles for each plot
    plot_keys = [:a, :b, :c, :d]  # Assuming you have these keys in your NamedTuple `pepe_for_makie`
    titles = ["Prey", "Top_predator", "Meso-predator", "Cell carrying capacity"]  # Custom titles for each plot

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]
    
    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        myname!(ax, frame[key]; interpolate=false)
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
    end
    
    # Optional: Add a colorbar to one of the axes to show the color scale
    # colorbar!(layout[1:2, end], axes[1], vertical=false)
end





