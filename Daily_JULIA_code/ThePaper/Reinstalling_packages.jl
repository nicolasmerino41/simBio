using Pkg

# STEP 1: Get all current dependencies
project = Pkg.project()
all_packages = keys(project.dependencies)

# STEP 2: Filter out packages you want to skip
exclude_packages = Set(["GraphCommunity", "DynamicGrids", "Dispersal"])
keep_packages = [pkg for pkg in all_packages if !(pkg in exclude_packages)]

# STEP 3: Remove all packages
Pkg.activate()
Pkg.resolve()
Pkg.instantiate()
Pkg.gc()  # force garbage collection of unused deps

# STEP 4: Remove all existing packages from the environment
for pkg in all_packages
    try
        Pkg.rm(pkg)
    catch err
        @warn "Couldn't remove $pkg" err
    end
end

# STEP 5: Re-add previously kept packages
for pkg in keep_packages
    try
        Pkg.add(pkg)
    catch err
        @warn "Couldn't re-add $pkg" err
    end
end

# STEP 6: Add dev versions from GitHub
Pkg.develop(PackageSpec(url="https://github.com/cesaraustralia/DynamicGrids.jl.git"))
Pkg.develop(PackageSpec(url="https://github.com/cesaraustralia/Dispersal.git"))

# Convert all_packages to a DataFrame
df_packages = DataFrame(Package = collect(all_packages))

# Save the DataFrame to a CSV file
CSV.write("all_packages.csv", df_packages)
# Read the csv back into a DataFrame
df_packages_readback = CSV.read("all_packages.csv", DataFrame)
