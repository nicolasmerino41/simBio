using Pkg
# Desktop PC
Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
# Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling" # Desktop
# meta_path = "C:\\Users\\nicol\\OneDrive\\PhD\\Metaweb Modelling" # Laptop

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
#################################################################################################
###################################FUNCTIONS###########################
#######################################################################
############### GROWTH ########################
function growth(abundance, self_regulation, K)
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

########################### CODE ############################
#############################################################
self_regulation = 0.01
bioclim_5 = Raster(joinpath(meta_path, "Rasters", "iberian_temperature.tif"))
spain = bioclim_5[X(-10 .. 4), Y(36 .. 45)]
# Smaller Raster for layout #####################
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"
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
for row in 1:size(npp_array, 1)
    for col in 1:size(npp_array, 2)
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
for row in 1:size(bio_npp_downsized, 1)
    for col in 1:size(bio_npp_downsized, 2)
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
for row in 1:size(init, 1)
    for col in 1:size(init, 2)
        if !isnan(init[row, col]) && (row in 45) && (col in 45)
            init[row, col] = npp_array[row, col] - 10
        elseif isnan(init[row, col]) || (!in(row, 45) || !in(col, 45))
            init[row, col] = 0.0
        end
    end
end
#################### INIT TUPLE ############################
transposed_init = permutedims(init, (2, 1))
transposed_npp = permutedims(npp_array, (2, 1))
##################### RULES ################################
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
cell = Cell() do data, a, I
    b = rand(100:1000)  # Generates a random integer between 100 and 1000
    max(a, a + growth(a, self_regulation, b))
end
## INWARDS DISPERSAL
indisp = InwardsDispersal(;
    formulation=ExponentialKernel(λ=1.125),
    distancemethod=AreaToArea(30)
)
## CUSTOM KERNEL
# custdisp = InwardsDispersal{:a, :a}(;
#     formulation=CustomKernel(3),
#     distancemethod=AreaToArea(10)
# )
 
## SETTINGS
tspan = 1.0:1000.0
ruleset = Ruleset(cell, indisp)

masklayer_for_gif = boolmask(transposed_npp)
masklayer_for_makie = boolmask(reverse(npp_array, dims=2))

################# OUTPUT for multiple rasters #############
#### PEPE
pepe_for_gif = (a = transposed_init, 
        b = transposed_npp
) 
pepe_for_makie = (a = reverse(init, dims=2),
        b = reverse(npp_array, dims=2)
)
#### RULES
cell_m = Cell{Tuple{:a,:b}, :a}() do data, (a, b), I
    max(a, a + growth(a, self_regulation, b))
end
## INWARDS DISPERSAL
# Define inwards dispersal for 'a', affecting 'a'
indisp_m = InwardsDispersal{:a, :a}(;
    formulation=ExponentialKernel(λ=0.00125),
    distancemethod=AreaToArea(30)
)

ruleset_m = Ruleset(cell_m, indisp_m)

## MAKIE
output = DG.MakieOutput(pepe_for_makie; tspan,
        ruleset = ruleset_m, 
        fps=100, 
        mask = masklayer_for_makie,
        maskcolor = RGB24(),
        minval=0,
        maxval=maximum(npp_array)+50,
)

# output = DG.MakieOutput(pepe_for_makie; tspan = tspan, ruleset = ruleset_m) do layout, frame, t
## GIF
output = GifOutput(pepe_for_gif;  
    filename="MultipleRasters.gif", 
    tspan=tspan, 
    fps=100, 
    mask = masklayer_for_gif,
    maskcolor = RGB24(0.0),
    minval=0, 
    maxval=1000, 
    scheme=ColorSchemes.rainbow,
    zerocolor=RGB24(0.3)
);
sim!(output, ruleset_m);

################# OUTPUT for single raster ################
## MAKIE
output = DG.MakieOutput(
    reverse(init, dims=2);
    tspan = tspan, 
    ruleset = ruleset, 
    fps = 1000, 
    mask = boolmask(reverse(npp_array, dims=2)),

)

output = MakieOutput(reverse(init, dims=2); tspan = tspan, ruleset = ruleset) do layout, frame, t
    image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno) 
end

using DynamicGrids
  using GLMakie
  
  # Define a rule, basic game of life
  ruleset = Ruleset(Life())
  # And the time-span for it to run
  tspan = 1:100
  
 
  
  # Create our own plots with Makie.jl
  output = MakieOutput(rand(Bool, 200, 300); tspan, ruleset) do layout, frame, t
      image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
  end

## GIF
output = GifOutput(transposed_init;  
    filename="SingleRaster.gif", 
    tspan=1.0:100.0, 
    fps=100, 
    mask = masklayer_for_gif,
    maskcolor = RGB24(0.0),
    minval=0, 
    maxval=1000, 
    scheme=ColorSchemes.rainbow,
    zerocolor=RGB24(0.3)
);
sim!(output, ruleset);
