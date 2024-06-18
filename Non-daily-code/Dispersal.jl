using Pkg
# Desktop PC
# Pkg.activate("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
# cd("C:\\Users\\MM-1\\OneDrive\\PhD\\JuliaSimulation\\simBio")
# Laptop
Pkg.activate("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio") 
cd("C:\\Users\\nicol\\OneDrive\\PhD\\JuliaSimulation\\simBio")

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
using Unitful: °C, K, cal, mol, mm
const DG = DynamicGrids
#################################################################################################
# Define parameters
rate = 0.2
carrycap = 100
nsteps = 10

ts = 1.0
tspan = 1.0:ts:100.0
# ts = Month(1);
# tspan = Date(Year(2022), Month(1)):ts:Date(Year(2030), Month(12))

# Create a LogisticGrowth rule with basic constructor
logistic_growth_basic = LogisticGrowth(rate=rate, carrycap=carrycap, timestep=timestep, nsteps_type=Int)
# Create a LogisticGrowth rule with parametric constructor
logistic_growth_custom = LogisticGrowth{Float64, Int}(rate=rate, carrycap=carrycap, timestep=timestep, nsteps_type=Int)

carrycap = 1e9
growth = LogisticGrowth(; rate=0.1, timestep=ts, carrycap=carrycap);
ruleset = Ruleset(growth; timestep=ts)

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

output = GifOutput(init_abundances_array;
    # Core keywords
    tspan=tspan,
    # Visualisation keywords
    scheme=ColorSchemes.tokyo, fps=10,
    minval=0.0, maxval=carrycap, 
    filename="sim.gif",
    print_ndims=true
);

sim!(output, ruleset)

output = REPLOutput(init_abundances_array; 
    tspan=tspan, 
    style=Block(), fps=5, color=:white, store=false
);
sim!(output, ruleset);
