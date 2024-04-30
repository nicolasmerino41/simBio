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
using Rasters, NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using Dates, DimensionalData, Distributions
using NCDatasets, RasterDataSources
using DynamicGrids
using Plots, ColorSchemes
# Dispersal
#################################################################################################

############# NOT NEEDED ############################
#####################################################
init = rand(Bool, 150, 200)
output = REPLOutput(init; tspan=1:200, fps=30)
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
##############################################################
##############################################################
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

################# MY SIMULATION ######################################
######################################################################
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
ruleset = Neighbors(Moore(1)) do data, neighborhood, cell, I
        
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

using DynamicGridsInteract, Blink
# Set up the output
output = DynamicGridsInteract.InteractOutput(init_abundances_array; 
    ruleset=ruleset,
    tspan=1:100, 
    store=false, 
)

# Run the simulation
caca = sim!(output, abundance_rule)

output = GifOutput(init_abundances_array;
    filename="population_growth.gif",  # Renamed file to reflect content
    tspan=1:100,  # Adjusted to create a range of timesteps
    fps=1,
    # scheme=ColorSchemes.rainbow,
    zerocolor=RGB24(0.0),
    minval=0.0,  # Ensures that the minimum value is zero
    maxval=1000.0,
    store=true  # Ensures that frames are stored
)
sim!(output, abundance_rule)

# output = REPLOutput(init_abundances_array; tspan=1:100, fps=30, color=Crayon(foreground=:red, background=:black, bold=true))

################## DISPERSAL #################################
##############################################################
##############################################################
# Generate random initial abundance values for each cell
num_cells = 5929
# Refactored code block to create a matrix of size 77x77
init_abundances_array = rand(10.0:2000.0, round(Int, sqrt(num_cells)), round(Int, sqrt(num_cells)))
# Get the dimensions of the init_abundances_array
grid_size = size(init_abundances_array)
# Create constant vectors for carrying capacity (K) and self-regulation
K_values = rand(100:1000.0,grid_size)

timestep = 1
tspan = 1:timestep:100
# Define your Gaussian kernel formulation
gaussian_kernel = GaussianKernel(0.1)
localdisp = InwardsDispersal(;
    radius=1,
    formulation=ExponentialKernel(; λ = 0.0125),
    distancemethod=AreaToArea(30),
)

logGrowth = LogisticGrowth(; rate=0.002, carrycap=1000, timestep=1)

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
##############################################################
init
################# Building the model again from scratch ########################
################################################################################
################################################################################
using DynamicGrids
using Dispersal
using GrowthMaps
# Creating the real npp array
npp_utm_df = CSV.read("npp_absolute_df.csv", DataFrame)
rename!(npp_utm_df, 3 => :npp, 4 => :X, 5 => :Y)
select(npp_utm_df, :X, :Y, :npp)
npp_values = npp_utm_df.npp
npp_values = push!(npp_values, 0.0)
# Create a 77x77 matrix filled with the values from npp_utm_df.npp
matrix_size = 77
npp_array = reshape(npp_values[1:end], matrix_size, matrix_size)
# Creating a simplified self-regulation array
self_regulation_array = fill(0.01, 77, 77)
# Creating abundance_array
init_abundances_array = rand(10.0:100.0, 77, 77)

# Random rule DOESN'T WORK
random = Cell{Tuple{:abun, :K, :self_reg}, :abun
    }() do abun, K, self_reg
    if abun > 20
        max(abun + growth(abun, self_reg, K), 0)
    else
        abun + 3
    end
end

# Define the rule for updating population abundance
abundance_rule = Neighbors(Moore(1)) do data, neighborhood, cell, I
        
    # Get the abundance, K, and self-regulation values for the current cell
    abundance = data[I...]
    K = npp_array[I...]  
    self_regulation = self_regulation_array[I...]
    
    # Update the abundance using the growth function
    # Apply a binomial event with p=0.5 to potentially reduce the abundance
    binomial_event = rand(Binomial(1, 0.5))
    new_abundance = abundance + growth(abundance, self_regulation, K) - (0.1 * abundance * binomial_event)
    
    # Ensure abundance stays non-negative
    new_abundance = max(new_abundance, 0.0)
    
    return new_abundance
end

# ArrayOutput
output_arr = ArrayOutput(init_abundances_array; tspan=1:100)
# REPLOutput
output_repl = REPLOutput(init_abundances_array; tspan=1:100)
# GifOutput
output_gif = GifOutput(init_abundances_array;
    filename="population_growth.gif",  # Renamed file to reflect content
    tspan=1:200,  # Adjusted to create a range of timesteps
    fps=0.1,
    scheme=ColorSchemes.algae,
    zerocolor=RGB24(0.0),
    minval=0.0,  # Ensures that the minimum value is zero
    maxval=maximum(npp_utm_df.npp),
    store=true  # Ensures that frames are stored
)
sim!(output_gif, abundance_rule)

######################################################################################
basedir = realpath(joinpath(dirname(Base.pathof(Rasters)), "../docs"))
humanpop_url = "https://github.com/cesaraustralia/DSuzukiiInvasionPaper/blob/master/data/population_density.tif?raw=true"
humanpop_filename = "population_density.tif"
humanpop_filepath = joinpath(basedir, humanpop_filename)
isfile(humanpop_filepath) || download(humanpop_url, humanpop_filepath);

australia_population = Rasters.Raster(humanpop_filepath)
plot(australia_population)

timestep = Month(1);
tspan = DateTime(2022, 1):timestep:DateTime(2030, 12)
aust = X(10500000 .. 15400000), Y(-5300000 .. -1000000)

australia_population_cropped = australia_population[X(10500000 .. 15400000), Y(-5300000 .. -1000000)]
spanish_population_cropped = australia_population[X(-950000 .. 370000), Y(4200000 .. 5100000)]
Rasters.dims(australia_population)
plot(spanish_population_cropped, title = "")
savefig("C:\\Users\\nicol\\.julia\\packages\\Rasters\\CrXfm\\docs\\build\\assets\\humanpop.png");


