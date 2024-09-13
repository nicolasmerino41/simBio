using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
# Pkg.add(path = "C:\\Users\\MM-1\\OneDrive\\PhD\\GitHub\\Dispersal.jl")
# Packages
using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames, Serialization
using StaticArrays, Distributions, Random
using Rasters #, DimensionalData
using DynamicGrids, Dispersal
using ColorSchemes, Colors, Crayons
using Makie, WGLMakie # ImageMagick, 
#using Unitful: °C, K, cal, mol, mm
# const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile

include("HerpsVsBirmmals.jl")