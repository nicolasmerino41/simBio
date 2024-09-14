packages = ["NCDatasets", "Shapefile", "ArchGDAL",
            "CSV", "DataFrames",
            "NamedArrays", "StaticArrays", "OrderedCollections",
            "Rasters", "RasterDataSources", "DimensionalData",
            "Dates", "Distributions", "Serialization",
            "Plots", "Colors", "Crayons", "ColorSchemes",
            "ImageMagick", "Makie", "WGLMakie",
            "Unitful", "StatsFuns"
]
##############################################################
# Retrieve a list of installed packages and their versions
pkg_list = Pkg.installed()
# Remove "DG" and "Dispersal" from the list if they exist
packages = filter(pkg -> pkg ∉ ["DynamicGrids", "Dispersal"], pkg_list)
# Extract only the package names
installed_pkg_names = [pkg for (pkg, ver) in pkg_list]

# Removing all the packages in the installed_pkg_names list
for pkg_name in installed_pkg_names
    Pkg.rm(pkg_name)
end
new_vector = append!(installed_pkg_names[1:3], installed_pkg_names[5:7])
new_vector = append!(new_vector, installed_pkg_names[9:end])
# Re-install all the packages from the list
for pkg_name in new_vector
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

