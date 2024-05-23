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
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, GLMakie, WGLMakie
# using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
utmraster = Raster("C:/Users/MM-1/OneDrive/PhD/JuliaSimulation/simBio/updated_utmraster.tif") 
pr = parent(utmraster)
MK.plot(utmraster, colormap = :inferno);
RS.values(utmraster)

species_df = CSV.File("Species_spain_df.csv") |> DataFrame

variables = species_df[!, 2:5]
rename!(variables, [:ID, :Value, :sum, :UTMCODE])

species = species_df[!, unique_species_in_web]
# reorder iberian_interact_df to match the column order of species_df
# iberian_interact_df = iberian_interact_df[!, sortperm(names(iberian_interact_df), by = x -> findfirst(isequal(x), names(species_df)))]
# serialize("iberian_interact_df.jls", iberian_interact_df)

species_df = hcat(variables, species, makeunique=true)
species_df_matrix = Matrix(species_df)

# Is this necessary?
for i in axes(species_df_matrix, 1), j in axes(species_df_matrix, 2)
    if species_df_matrix[i, j] == 0
        species_df_matrix[i, j] = 0.0
    elseif species_df_matrix[i, j] == 1
        species_df_matrix[i, j] = 1.0
    end
end

utmraster_DA = DimArray(utmraster)
utmraster_da = map(x -> isnothing(x) || isnan(x) ? false : true, utmraster_DA)
MK.plot(utmraster_DA);

DA = DimArray(reshape([MyStructs256(SVector{256, Float64}(fill(0.0, 256))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            DA[j] = MyStructs256(SVector{256, Float64}(species_df_matrix[i, 5:260]))
        end
    end
end

# This is for visualising richness in the raster or for creating a boolmask
DA_sum = DimArray(reshape(fill(0.0, 125*76), 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_sum = falses(dims(DA_with_abundances))
for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            DA_sum[j] =  true #species_df_matrix[i, 3]
        # else
        #     DA_sum[j] = false
        end
    end
end
MK.plot(reverse(DA_sum_r, dims=1), colormap = :redsblues)
DA_sum_r = reverse(DA_sum, dims=1)
DA_sum_p = permutedims(DA_sum, (2, 1))

###### Let's do a small example of the simBio with the actual ATLAS data ######
DA_with_abundances = deepcopy(DA)
# for row in axes(DA, 1), col in axes(DA, 2)
#     if DA[row, col] != MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
#         new_a = SVector{256, Float64}([DA[row, col].a[i] != 0.0 ? 100 : DA[row, col].a[i] for i in 1:256])
#         DA_with_abundances[row, col] = MyStructs256(new_a)
#     end
# end
# serialize("DA_with_abundances.jls", DA_with_abundances)
# serialize("DA_with_abundances_all100.jls", DA_with_abundances)
# Load the serialized object from the file and cast to the specific type
DA_with_abundances = deserialize("DA_with_abundances_all100.jls")::DimArray{MyStructs256{Float64},2}

DA_with_abundances_r = reverse(DA_with_abundances, dims=1)
DA_with_abundances_p = permutedims(DA_with_abundances, (2, 1))
DA_with_abundances_p_masked = deepcopy(DA_with_abundances_p)
for row in axes(DA_with_abundances_p, 1), col in axes(DA_with_abundances_p, 2)
    if DA_with_abundances_p[row, col] == 0.0
        DA_with_abundances_p_masked[row, col] = MyStructs256(SVector{256, Float64}(fill(-100.0, 256)))
    end
end

######################### NPP ####################################
npp_absolute = CSV.File(joinpath(meta_path, "npp_absolute_df.csv")) |> DataFrame
rename!(npp_absolute, [:ID, :UTMCODE, :npp]) 
npp_absolute_in_kg = deepcopy(npp_absolute)
npp_absolute_in_kg.npp = npp_absolute.npp .* 1000
npp_absolute_in_kg = npp_absolute_in_kg[:, [2, 3]]

species_df = leftjoin(species_df, npp_absolute_in_kg, on = :UTMCODE)
species_df_matrix = Matrix(species_df)

npp_DA = DimArray(zeros(125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            npp_DA[j] = species_df_matrix[i, 261]
        end
    end
end
MK.plot(npp_DA_p, colormap = :darkrainbow)
npp_DA_r = reverse(npp_DA, dims=1)
npp_DA_p = permutedims(npp_DA, (2, 1))
########################## RULES #################################
################################################################## 
# Merge rule with non-negative constraint for MyStructs256
merge_rule = Cell{}() do data, state, I
    
    merged_state = state + MyStructs256(SVector{256}(
            intrinsic_growth_256(state, 0.01, npp_DA_p[I...]).a .+ # You can switch the 100-0 by
            # transposed_npp[I...] if needed 
            trophic_optimized(state, full_IM).a)
        )
    
    # if any(merged_state.a .< 0.0) || any(isnan.(merged_state.a))
    #     println("merged_state is NA", state, merged_state)
    # end
    return MyStructs256(max.(merged_state.a, 0.0))
    # if state.b > 26000.0
    #     println(state.b)
    # end
end

indisp = InwardsDispersal{}(;
    formulation=pepe_kernel(λ = 0.1),
    distancemethod=AreaToArea(30)
)

ruleset = Ruleset(merge_rule, indisp)
full_IM = results[1][0.001]
# full_IM = full_IM_list[2][10]
pepe = (grid1 = DA_with_abundances_p, grid2 = npp_DA_p)
gif_output = GifOutput(
    DA_with_abundances_p; tspan = 1:100,
    filename = "mystruct256_Atlas.gif",
    fps = 20, 
    scheme=ColorSchemes.darkrainbow,
    zerocolor=RGB24(0.0),
    minval = 0.0, maxval = 10000,
    mask = DA_sum_p
)
@time saved = sim!(gif_output, ruleset)
 
array_output = ArrayOutput(
    DA_with_abundances_p; tspan = 1:10,
    mask = DA_sum_p
)
r = sim!(array_output, ruleset)
uno = r[1]
uno_dos = r[2]
uno_diez = r[10]
uno_tres = r[3]
max_values = MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
for i in 1:length(r)
    maxu = maximum(r[i])
    if maxu > max_values #&& maxu < 40000.0
        max_values = maxu
        println(max_values.b, " in time ", i)
    end
end

Base.@kwdef struct pepe_kernel{P} <: KernelFormulation
    λ::P = Param(1.0, bounds=(0.0, 2.0))
end

(f::pepe_kernel)(d) = f.λ
(f::pepe_kernel)(d) = exp(-d / (2 * f.λ^2))
function (kernel::pepe_kernel)(distance)
    value = exp(-distance / (2 * kernelλ^2))

    return value < 1.0 ? 0.0 : value
end
