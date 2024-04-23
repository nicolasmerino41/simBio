using Pkg

Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio")
cd("C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\simBio") 

using DynamicGrids, Crayons, DynamicGridsGtk, Plots, ColorSchemes
using GeoData, NCDatasets, ArchGDAL, RasterDataSources
using Dispersal, Dates, Plots, GrowthMaps, Unitful, Statistics
using Unitful: °C, K, cal, mol, mm
using Rasters
using RasterDataSources
using GeoInterface
############# NOT NEEDED ############################
init = rand(Bool, 150, 200)
output = REPLOutput(init; tspan=1:200, fps=30, color=Crayon(foreground=:red, background=:black, bold=true))
sim!(output, Life())

# Or define it from scratch (yes this is actually the whole implementation!)
life = Neighbors(Moore(1)) do data, hood, state, I
    born_survive = (false, false, false, true, false, false, false, false, false), 
                   (false, false, true, true,  false, false, false, false, false)
    born_survive[state + 1][sum(hood) + 1]
end
sim!(output, life)

init = rand(Float32, 100, 100)
output = ArrayOutput(init; tspan=1:100)
sim!(output, ruleset; init=init)

########### Example Forest Fire ##############################
using DynamicGrids, ColorSchemes, Colors, BenchmarkTools, ImageMagick

DEAD, ALIVE, BURNING = 1, 2, 3

neighbors_rule = let prob_combustion=0.0001, prob_regrowth=0.01
    Neighbors(Moore(1)) do data, neighborhood, cell, I
        if cell == ALIVE
            if BURNING in neighborhood
                BURNING
            else
                rand() <= prob_combustion ? BURNING : ALIVE
            end
        elseif cell == BURNING
            DEAD
        else
            rand() <= prob_regrowth ? ALIVE : DEAD
        end
    end
end

# Set up the init array and output (using a Gtk window)
init = fill(ALIVE, 400, 400)
output = GifOutput(init; 
    filename="forestfire.gif", 
    tspan=1:200, 
    fps=10, 
    minval=DEAD, maxval=BURNING, 
    scheme=ColorSchemes.rainbow,
    zerocolor=RGB24(0.0)
)

# Run the simulation, which will save a gif when it completes
sim!(output, neighbors_rule)

################# MY SIMULATION ######################
# Define the growth function
function growth(abundance, self_regulation, K)
    return self_regulation * K * abundance * (1 - abundance / K)
end

# Generate random initial abundance values for each cell
num_cells = 5929
# Refactored code block to create a matrix of size 77x77
init_abundances_array = rand(10.0:100.0, round(Int, sqrt(num_cells)), round(Int, sqrt(num_cells)))

# Get the dimensions of the init_abundances_array
grid_size = size(init_abundances_array)

# Create constant vectors for carrying capacity (K) and self-regulation
K_values = fill(rand(100:1000.0), grid_size)
self_regulation_values = fill(0.01, grid_size)

# Define the rule for updating population abundance
abundance_rule = Neighbors(Moore(1)) do data, neighborhood, cell, I
        
    # Get the abundance, K, and self-regulation values for the current cell
    abundance = data[I...]
    K = K_values[I...]
    self_regulation = self_regulation_values[I...]
    
    # Update the abundance using the growth function
    new_abundance = abundance + growth(abundance, self_regulation, K)
    
    # Ensure abundance stays non-negative
    new_abundance = max(new_abundance, 0.0)
    
    return new_abundance
end

# Set up the output
output = ArrayOutput(init_abundances_array; tspan=1:100)

# Run the simulation
caca = sim!(output, abundance_rule)

output = GifOutput(init_abundances_array;
    filename="population_growth.gif",  # Renamed file to reflect content
    tspan=1:100,  # Adjusted to create a range of timesteps
    fps=1,
    scheme=ColorSchemes.rainbow,
    zerocolor=RGB24(0.0),
    minval=0.0,  # Ensures that the minimum value is zero
    maxval=1000.0,
    store=true  # Ensures that frames are stored
)
sim!(output, abundance_rule)

# output = REPLOutput(init_abundances_array; tspan=1:100, fps=30, color=Crayon(foreground=:red, background=:black, bold=true))

################## DISPERSAL #################################
# Generate random initial abundance values for each cell
num_cells = 5929
# Refactored code block to create a matrix of size 77x77
init_abundances_array = rand(10.0:2000.0, round(Int, sqrt(num_cells)), round(Int, sqrt(num_cells)))
# Get the dimensions of the init_abundances_array
grid_size = size(init_abundances_array)
# Create constant vectors for carrying capacity (K) and self-regulation
K_values = rand(100:1000.0,grid_size)
timestep = Month(1);
tspan = DateTime(2022, 1):timestep:DateTime(2030, 12)
# Define your Gaussian kernel formulation
gaussian_kernel = GaussianKernel(0.1)
localdisp = InwardsDispersal(;
    radius=1,
    formulation=ExponentialKernel(; λ = 0.0125),
    distancemethod=AreaToArea(30),
)

logGrowth = LogisticGrowth(; rate=0.002, carrycap=1000, timestep=Day(1))

# Use tspan consistent with timestep
ruleset = Ruleset(localdisp, logGrowth, timestep=timestep)

# Make sure tspan matches the timestep
output = ArrayOutput(init_abundances_array; tspan=tspan)

sim!(output, ruleset; init=init_abundances_array)

output = GifOutput(init_abundances_array;
    # Core keywords
    tspan=tspan,
    # Visualisation keywords
    scheme=ColorSchemes.rainbow, 
    fps=1,
    minval=0.0, maxval=2000.0, 
    filename="dispersal.gif"
);

sim!(output, ruleset; init=init_abundances_array)

#### Creating our own kernel
struct MyCustomKernelFormulation{Rho} <: KernelFormulation
    ro::Rho  # Rho is a type parameter for the ro value in your formula
end
# Implement the method for calculating the probability density based on distance
(k::MyCustomKernelFormulation)(distance) = exp(-distance^2 / (2 * k.ro^2))


# Create an instance of your custom kernel formulation with a specific ro value
custom_formulation = MyCustomKernelFormulation(0.2)  # Example value for ro
# Create a dispersal kernel using your custom formulation
custom_kernel = DispersalKernel(formulation=custom_formulation)
localdisp = InwardsDispersal(; radius=1, formulation=custom_kernel, distancemethod=AreaToArea(30))

DispersalKernel(formulation=custom_formulation)

################ Ralf Example ####################
struct Lon
    range::Tuple{Int, Int}
end
struct Lat
    range::Tuple{Int, Int}
end
function Between(a, b)
    return (a, b)
end

basedir = realpath(joinpath(dirname(Base.pathof(Dispersal)), "../docs"))
timestep = Month(1);
tspan = DateTime(2022, 1):timestep:DateTime(2030, 12)
aust = Lon(Between(113, 154)),
       Lat(Between(-45, -10));

localdisp = InwardsDispersal(;
    radius=1, formulation=ExponentialKernel(; λ = 0.0125),
    distancemethod=AreaToArea(30))

allee = AlleeExtinction(; minfounders=10.0);

humanpop_url = "https://github.com/cesaraustralia/DSuzukiiInvasionPaper/blob/master/data/population_density.tif?raw=true"
humanpop_filename = "population_density.tif"
humanpop_filepath = joinpath(basedir, humanpop_filename)
isfile(humanpop_filepath) || download(humanpop_url, humanpop_filepath);

humanpop = Raster(humanpop_filepath; mappedcrs=EPSG(4326))[At(Band{Int64}(1)), aust...] |>
    A -> replace_missing(A, missing) |>
    A -> permutedims(A, (Lat, Lon)) |> 
    A -> reorder(A, Lat=ReverseArray, Lon=ForwardArray)

    plot(humanpop)
savefig(joinpath(basedir, "build/assets/humanpop.png"));

cucu = Raster(humanpop_filepath; mappedcrs=EPSG(4326))[At(Band{Int64}(1)), aust...]
plot(cucu)

raster_data = Raster(humanpop_filepath; mappedcrs=EPSG(4326))
# Define the bounding box
top_left_x = 
top_left_y = ...
bottom_right_x = ...
bottom_right_y = ...
aust = Lon(Between(113, 154)),
       Lat(Between(-45, -10));
       
A = raster_data
aust = A[X(113 .. 155), Y(-40 .. -25), RS.Band(1)]

shp = Shapefile.Table("all_primary_countries/all_primary_countries.shp") |> DataFrame
plot(shp.geometry)

raster = Raster("Harmonized_DN_NTL_2020_simVIIRS.tif")
# plot(harmon)

raster[X(Near(77.1025)), Y(Near(28.7041)), Band(1)] # value of image near longitude = 77.1025 and latitude = 28.7041
raster[X(Near(77.1025)), Y(Near(28.7041)), Band(1)] |> Int 

shp = dropmissing(shp)
# plot(shp.geometry)

australia_index = findall(x -> x == "Australia", shp.name)
# Doesn't work| australia = Rasters.crop(raster; to = shp.geometry[australia_index])

bounds = X(Rasters.Between(113, 154)), Y(Rasters.Between(-45, -10))
australia = raster[bounds...] # |> plot
# Count the number of non-missing pixels in the raster
pixel_count = count(!ismissing, australia)

australia_masked = Rasters.mask(australia, with = shp.geometry[australia_index], missingvalue = 0)
sol_by_country = Rasters.zon(sum, raster; of=shp.geometry)

### Load and plot the file
# Set the environment variable for raster data sources to a default path
ENV["RASTERDATASOURCES_PATH"] = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"

### Mask
awap_masked = mask(awap; with=wc_mask)
b = plot(awap_masked; clims=(10, 45))

