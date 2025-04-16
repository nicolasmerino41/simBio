begin
    PC = "nicol"
    const EXTINCTION_THRESHOLD = 1e-6
    using Pkg
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
    meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
    using ArchGDAL #, Shapefile, NCDatasets
    using CSV, DataFrames, NaNMath
    using NamedArrays, StaticArrays, OrderedCollections
    using Rasters, RasterDataSources #, DimensionalData
    using DynamicGrids, Dispersal
    using Dates, Distributions, Serialization, StatsBase, JLD2, Random
    using ColorSchemes, Colors #Crayons, 
    using Makie, CairoMakie # ImageMagick
    using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Logging, ForwardDiff
    # if end_true
    #     using DifferentialEquations # EcologicalNetworksDynamics, 
    # end
    # using Plots
    const DG, MK, AG, RS, Disp, DF = DynamicGrids, Makie, ArchGDAL, Rasters, Dispersal, DataFrames
    const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
end
