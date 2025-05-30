@time begin
    PC = "nicol"
    const EXTINCTION_THRESHOLD = 1e-2
    using Pkg
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\ThePaper"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
    meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
    # using ArchGDAL #, Shapefile, NCDatasets
    using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs, Serialization
    import Base.Threads: @threads
    # include("Ladder/Scenario2/Ladder4.1.jl")
    # using CSV, DataFrames
    # using NamedArrays, StaticArrays, OrderedCollections
    # # using Rasters, RasterDataSources #, DimensionalData
    # # using DynamicGrids, Dispersal
    # using Dates, Distributions, Serialization, StatsBase, Random, GLM #JLD2
    # using ColorSchemes, Colors #Crayons, 
    using CairoMakie # ImageMagick
    # using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Logging, ForwardDiff
    # using Graphs, CategoricalArrays
    # # if end_true
    # #     using DifferentialEquations # EcologicalNetworksDynamics, 
    # # end
    # # using Plots
    # const MK, DF = Makie, DataFrames
    # const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
end

# Pkg.add("DifferentialEquations")
# Pkg.add("DiffEqCallbacks")
# Pkg.add("LinearAlgebra")
# Pkg.add("Logging")
# Pkg.add("ForwardDiff")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("NamedArrays")
# Pkg.add("StaticArrays")
# Pkg.add("OrderedCollections")
# Pkg.add("Dates")
# Pkg.add("Statistics")
# Pkg.add("Distributions")
# Pkg.add("Serialization")
# Pkg.add("StatsBase")
# Pkg.add("Random")
# Pkg.add("GLM")