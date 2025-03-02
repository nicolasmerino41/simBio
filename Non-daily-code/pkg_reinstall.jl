packages = [
    "ArchGDAL", "BenchmarkTools", "CSV", "CairoMakie", "ColorSchemes", 
    "Colors", "Crayons", "DataFrames", "Dates", "DimensionalData",
    "Distributions", "ImageMagick", "JLD2", 
    "Makie", "NCDatasets", "NamedArrays", "OrderedCollections", "PDFmerger", 
    "Plots", "ProfileView", "RasterDataSources", "Rasters", "Serialization", 
    "Shapefile", "StaticArrays", "StatsBase", "StatsFuns", "Stencils", 
    "Unitful", "WGLMakie"
]
##############################################################
# Retrieve a list of installed packages and their versions
pkg_list = Pkg.installed()
# Remove "DG" and "Dispersal" from the list if they exist
packages = filter(pkg -> pkg ∉ ["DynamicGrids", "Dispersal", "OrdinaryDiffEq"], installed_pkg_names)
# Extract only the package names
installed_pkg_names = [pkg for pkg in pkg_list]

# Removing all the packages in the installed_pkg_names list
for pkg_name in installed_pkg_names
    Pkg.add(pkg_name)
end
new_vector = append!(installed_pkg_names[1:2], installed_pkg_names[4:5])
new_vector = append!(new_vector, installed_pkg_names[7:end])
# Re-install all the packages from the list
for pkg_name in packages
    Pkg.add(pkg_name)
end
Pkg.add(PackageSpec(url="https://github.com/cesaraustralia/DynamicGrids.jl", rev="dev"))
Pkg.add(PackageSpec(url="https://github.com/cesaraustralia/Dispersal.jl", rev="dev"))

Pkg.add("C:\\Users\\nicol\\OneDrive\\PhD\\GitHub\\Stencils.jl")
Pkg.add("C:\\Users\\nicol\\OneDrive\\PhD\\GitHub\\DynamicGrids.jl")
rafs = ["Rasters", "Unitful", "ColorSchemes", "Blink", "ImageMagick", "CSV", "DataFrames", "Crayons", "Colors", "Plots"]
for pkg_name in packages
    Pkg.add(pkg_name)
end

for pkg_name in installed_pkg_names
    Pkg.free()
end


using Pkg

# Get a list of all installed packages
installed_packages = keys(Pkg.project().dependencies)
installed_packages_vector = collect(installed_packages)
installed_packages_vector = push!(installed_packages_vector, "GLM")
println(installed_packages_vector)

# List of packages to exclude
packages_to_remove = ["DynamicGrids", "Dispersal", "CUDA", "EcologicalNetworkDynamics"]

# Remove specified packages from the list
updated_packages = filter(x -> !(x in packages_to_remove), installed_packages_vector)
println(updated_packages)

for pkg_name in installed_packages_vector
    Pkg.rm(pkg_name)
end

for pkg_name in updated_packages
    Pkg.add(pkg_name)
end
