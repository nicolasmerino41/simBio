@time begin
    PC = "nicol"
    const EXTINCTION_THRESHOLD = 1e-6
    using Pkg
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\ThePaper"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
    meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
    # using ArchGDAL #, Shapefile, NCDatasets
    using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
    using Distributions, Serialization
    import Base.Threads: @threads
    
    using CairoMakie
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