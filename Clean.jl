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
using NamedArrays, StaticArrays
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
###################################FUNCTIONS###########################
#######################################################################
############### GROWTH ########################
function growth(abundance::AbstractFloat, self_regulation::AbstractFloat, K::AbstractFloat)
    # Assuming growth is a simple function of abundance, self-regulation, and carrying capacity (K)
    return self_regulation * K * abundance * (1 - abundance / K)
end
################### SIMBIO_EFFICIENT_MAP ######################
###############################################################
function simbio_efficient_map!(K, self_regulation, abundances, A_matrix, num_steps)
    num_species = length(abundances)
    
    # Check if lengths match
    if num_species != length(self_regulation) || num_species != length(K)
        throw(ArgumentError("The number of species must be equal to the length of the self_regulation and K vectors"))
    end

    for t in 1:num_steps
        # Calculate growth rates
        growth_rates = growth.(abundances, self_regulation, K)
        # println(growth_rates)
        # Calculate total change in abundances
        delta_abundances = growth_rates .+ sum(A_matrix .* self_regulation .* adjoint(abundances) .*abundances, dims=2)[:, 1]
        # println(delta_abundances[])
        # Update abundances in-place
        @. abundances += delta_abundances
        # println(abundances)
    end
    
    return abundances
end

# aaaa = NamedArray(aaa, (species_names, species_names), ("Species", "Species"))
# function lv_function(abundance, K, A_matrix, species::AbstractString)
#     return abundance + growth(abundance, self_regulation, K) .+ sum(parent(A_matrix[species, :]) .* self_regulation .* adjoint(abundances) .*abundances, dims=2)
# end
# parent(aaaa["Anthus campestris", :])
########################### CODE ############################
#############################################################
self_regulation = 0.001
bioclim_5 = Raster(joinpath(meta_path, "Rasters", "iberian_temperature.tif"))
spain = bioclim_5[X(-10 .. 4), Y(36 .. 45)]
# Smaller Raster for layout #####################
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling"
bioclim_paths = RasterDataSources.getraster(WorldClim{BioClim}, (5,7,8,12))
bioclim_stack = RasterStack(WorldClim{BioClim}, (5, 7, 8, 12), res="10m")
bioclim_stack = RasterStack(bioclim_paths)
bioclim_5 = bioclim_stack[:bio5]
spain = bioclim_5[X(-10 .. 4), Y(36 .. 45)]
MK.plot(spain);

################### NPP #########################
################ Random NPP ##################
npp_array = deepcopy(spain)
npp_array = replace_missing(npp_array, 0)
for row in axes(npp_array, 1)
    for col in axes(npp_array, 2)
        if npp_array[row, col] != 0 #&& (row in 90:100) && (col in 90:100)
            npp_array[row, col] = rand(1000.0:10000.0)
        end
    end
end
# sizi = (100,100)
# npp_array = rand(100.0:1000.0, sizi)
MK.plot(npp_array);

############## Real NPP ####################
bio_npp = RS.Raster(joinpath("npp_cropped.tif"))
MK.plot(bio_npp);
bio_npp_downsized = RS.resample(bio_npp, to = spain)
MK.plot(bio_npp_downsized);
for row in axes(bio_npp_downsized, 1)
    for col in axes(bio_npp_downsized, 2)
        if spain[row, col] == -3.4f38 || bio_npp_downsized[row, col] == 0.0
            bio_npp_downsized[row, col] = NaN
        end
    end 
end
MK.plot(npp_array);
npp_array = replace_missing(bio_npp_downsized, 0)
# TODO, Remove Africa and Balear Islands
# spain_shape = Shapefile.Table("IP_shape.shp") 
# bio_npp_downsized = RS.crop(bio_npp_downsized, to = spain_shape)
# MK.plot(spain_shape.geometry)

################## INIT ################################
init = deepcopy(npp_array)
init_trimmed = Rasters.trim(spain) # This basically creates the smallest bounding box of non-Na values
init = replace_missing(init, 0) # We need to replace the missing for 0s or it'l spread everywhere
# init = rand(10.0:100.0, sizi)

function create_species_inits(init_raster, npp_raster, values)
    species_inits = []
    for value in values
        species_init = deepcopy(init_raster)
        for row in axes(species_init, 1), col in axes(species_init, 2)
            if !isnan(species_init[row, col]) && (row in 45:50) && (col in 45:50)
                species_init[row, col] = npp_raster[row, col] * (1 - value)
            else
                species_init[row, col] = 0.0
            end
        end
        push!(species_inits, species_init)
    end
    return species_inits
end

values = [0.2, 0.5, 0.8]
inits = create_species_inits(init, npp_array, values)

#### INIT TUPLE
# For GIF
transposed_init = permutedims(inits[1], (2, 1))
transposed_npp = permutedims(npp_array, (2, 1))
# For Makie
reversed_init = reverse(init, dims=2) 
reversed_npp = reverse(npp_array, dims=2)
reversed_inits = [reverse(inits[i], dims=2) for i in 1:length(inits)]
## SETTINGS
tspan = 1.0:100.0
masklayer_for_gif = boolmask(transposed_npp)
masklayer_for_makie = boolmask(reverse(npp_array, dims=2))

################# OUTPUT for multiple rasters #############
###########################################################
###########################################################
###########################################################
#### PEPE
pepe_for_gif = (
    a = transposed_init, 
    b = transposed_init,
    c = transposed_init,
    d = transposed_npp
) 
pepe_for_makie = (
    a = reversed_inits[1], 
    b = reversed_inits[2],
    c = reversed_inits[3],
    # d = reversed_inits[4],
    # e = reversed_inits[5],
    # f = reversed_inits[6],
    # g = reversed_inits[7],
    # h = reversed_inits[8],
    # i = reversed_inits[9],
    # j = reversed_inits[10],
    # k = reversed_inits[11],
    d = reversed_npp)
#### RULES
## CELL
cell_m = Cell{Tuple{:a,:b}, :a}() do data, (a, b), I
    return a + growth(a, self_regulation, b)
end
cell_m_sp1 = Cell{Tuple{:a,:b,:d}, :a}() do data, (a, b, d), I 
    # if I == (30,20)
    #     println("a: ", a, " b: ", b, " d: ", d)
    # end 
    return max(0.0, a  + growth(a, self_regulation, d) - 0.005 * a * b )
end

cell_m_sp2 = Cell{Tuple{:a,:b,:c,:d}, :b}() do data, (a, b, c, d), I
    return max(0.0, b + growth(b, self_regulation, d * 0.2) + 0.005 * b * a - 0.005 * b * c)
end

cell_m_sp3 = Cell{Tuple{:b,:c,:d}, :c}() do data, (b, c, d), I
    return max(0.0, c + growth(c, self_regulation, d * 0.05) + 0.005 * c * b)
end

## INWARDS DISPERSAL
# Define inwards dispersal for 'a', affecting 'a'
indisp_m = InwardsDispersal{Tuple{:a,:b,:c}, Tuple{:a,:b,:c}}(;
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
# Define inwards dispersal for species 'a', affecting only 'a'
indisp_m_sp1 = InwardsDispersal{:a, :a}(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
# Define inwards dispersal for species 'b', affecting only 'b'
indisp_m_sp2 = InwardsDispersal{:b, :b}(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
# Define inwards dispersal for species 'c', affecting only 'c'
indisp_m_sp3 = InwardsDispersal{:c, :c}(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)

ruleset_m = Ruleset(cell_m_sp1, cell_m_sp2, cell_m_sp3, indisp_m_sp1, indisp_m_sp2, indisp_m_sp3)
###################     # VISUALISATION ####################
## MAKIE
# Default way
# output = DG.MakieOutput(pepe_for_makie; tspan,
#         ruleset = ruleset_m, 
#         fps=1, 
#         mask = masklayer_for_makie,
#         maskcolor = RGB24(),
#         minval=0,
#         maxval= maximum(transposed_npp) + 50
# )

# # Setup for multi-grid simulation using NamedTuple of Observables
# pepe_for_makie = (a = Observable(reversed_init), b = Observable(reversed_npp))
# Personalised way
output = MakieOutput(pepe_for_makie; tspan=1:1000, ruleset=ruleset_m,
        mask=masklayer_for_makie, fps=5) do (; layout, frame)

    # Define the global color limits
    color_limits = (10.0, 2000.0)

    # Setup the keys and titles for each plot
    plot_keys = [:a, :b, :c, :d]  # Assuming you have these keys in your NamedTuple `pepe_for_makie`
    titles = ["Prey", "Top_predator", "Meso-predator", "Cell carrying capacity"]  # Custom titles for each plot

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]
    
    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        Makie.image!(ax, frame[key]; interpolate=false, colormap=:inferno, clims=color_limits)
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
    end
    
    # Optional: Add a colorbar to one of the axes to show the color scale
    # colorbar!(layout[1:2, end], axes[1], vertical=false)
end

## GIF
output = GifOutput(pepe_for_gif;  
    filename="MultipleRasters.gif", 
    tspan=tspan, 
    fps=1000, 
    mask = masklayer_for_gif,
    maskcolor = RGB24(0.0),
    minval=0, 
    maxval=maximum(transposed_npp), 
    scheme=ColorSchemes.inferno,
    zerocolor=RGB24(0.0)
);
sim!(output, ruleset_m);

################# OUTPUT for single raster ################
###########################################################
###########################################################
###########################################################
"""
This is fully operating, it is very important to ensure type_stability, or it will not work.
Notify any change from now on to keep the default version always avalaible.
"""
#### RULES
## NEIGHBOR
neighbor = DG.Neighbors(Moore(1)) do data, neighborhood, cell, I
    abundance = data[I...]
    K = npp_array[I...]
    
    binomial_event = rand(Binomial(1, 0.5))
    # leaving = cell * binomial_event * 0.01
    # entering = sum(neighborhood) * 0.005

    new_abundance = cell + growth(cell, self_regulation, K+0.0000001) 
    new_abundance = max(new_abundance, 0.0)
    data[I...] = new_abundance
end

## CELL
#= (it can work with reversed_npp[I...] now thanks to type-Float64 in growth()) but you
need a different ruleset for makie or gif due to I[...] coordinates being reversed =#
cell_s_makie = Cell() do data, cell, I
      # Generates a random integer between 100 and 1000
      # println("a: ", a, " npp: ", reversed_npp[I...], " growth: ", growth(a, self_regulation, reversed_npp[I...]))
      return a + growth(cell, self_regulation, reversed_npp[I...])
end
cell_s_gif = Cell() do data, a, I
      return a + growth(a, self_regulation, transposed_npp[I...])
end
MK.plot(reversed_init);

## INWARDS DISPERSAL
indisp_s = InwardsDispersal(
        formulation=ExponentialKernel(λ = 1.125),
        distancemethod = AreaToArea(30)
)

############### VISUALISATION ############################
tspan = 1.0:1000.0
## MAKIE
ruleset_s_makie = Ruleset(cell_s_makie, indisp_s) # You need a different ruleset for makie or gif
# Default way
# output = DG.MakieOutput(
#     reversed_init;
#     tspan = tspan, 
#     ruleset = ruleset_s_makie, 
#     fps = 10, 
#     mask = masklayer_for_makie
# );

# Personalised way
output = MakieOutput(
        reversed_init; 
        tspan = tspan,
        ruleset = ruleset_s_makie,
        mask = masklayer_for_makie) do (; layout, frame)
            ax = Axis(layout[1, 1])
            image!(ax, frame; interpolate=false, colormap=:inferno)
end

## GIF At least this one works
ruleset_s_gif = Ruleset(cell_s_gif, indisp_s)
output = GifOutput(transposed_init;  
    filename="SingleRaster.gif", 
    tspan=tspan, 
    fps=10, 
    mask = masklayer_for_gif,
    maskcolor = RGB24(0.0),
    minval=0, 
    maxval=maximum(transposed_npp), 
    scheme=ColorSchemes.inferno,
    zerocolor=RGB24(0.0)
);

sim!(output, ruleset_s_gif);
## LIFE
ruleset_life = Ruleset(Life())
# And the time-span for it to run
tspan = 1:100
# Create our own plots with Makie.jl
output = MakieOutput(rand(Bool, 200, 300); tspan, ruleset=ruleset_life) do (; layout, frame)
    image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
end


