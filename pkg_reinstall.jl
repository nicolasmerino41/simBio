packages = ["Rasters", "NCDatasets", "Plots",
            "Shapefile", "ArchGDAL", "CSV", 
            "DataFrames", "Dates", "DimensionalData",
            "RasterDataSources", "IntervalSets", "GrowthMaps",
            "Dispersal", "DynamicGrids", "Unitful",
            "Plots", "ColorSchemes", "Blink",
            "Crayons"]
##############################################################
Pkg.rm("GrowthMaps")
# Retrieve a list of installed packages and their versions
pkg_list = Pkg.installed()

# Extract only the package names
installed_pkg_names = [pkg for (pkg, ver) in pkg_list]

# Removing all the packages in the installed_pkg_names list
for pkg_name in installed_pkg_names
    Pkg.rm(pkg_name)
end
# Re-install all the packages from the list
for pkg_name in installed_pkg_names
    Pkg.add(pkg_name)
end

for pkg_name in packages
    Pkg.add(pkg_name)
end
