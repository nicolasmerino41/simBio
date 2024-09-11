using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")

# Packages
using NCDatasets, Shapefile, ArchGDAL
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase
using Plots
using Colors, Crayons, ColorSchemes
using ImageMagick, Makie, WGLMakie
using Unitful: °C, K, cal, mol, mm
const DG, MK, PL, AG, RS, Disp, DF, NCD, SH = DynamicGrids, Makie, Plots, ArchGDAL, Rasters, Dispersal, DataFrames, NCDatasets, Shapefile
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
#################################################################################################
utmraster = Raster(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio\\updated_utmraster.tif")) 
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

###### Let's do a small example of the simBio with the actual ATLAS data ######
DA_with_abundances = deepcopy(DA)
for row in axes(DA, 1), col in axes(DA, 2)
    if DA[row, col] != mMyStructs256(Vector{Float64}(fill(0.0, 256)))
        new_a = Vector{Float64}([DA[row, col].a[i] != 0.0 ? 100 : DA[row, col].a[i] for i in 1:256])
        DA_with_abundances[row, col] = mMyStructs256(new_a)
    end
end
# serialize("DA_with_abundances.jls", DA_with_abundances)
# serialize("DA_with_abundances_all100.jls", DA_with_abundances)
# serialize("DA_with_abundances_all100_mMyStructs256.jls", DA_with_abundances)
# Load the serialized object from the file and cast to the specific type
DA_with_abundances = deserialize("DA_with_abundances_all100.jls")::DimArray{MyStructs256{Float64},2}

# This is for visualising richness in the raster or for creating a boolmask
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
DA_richness = DimArray(reshape(fill(0.0, 125*76), 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_richness_herps = DimArray(reshape(fill(0.0, 125*76), 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
DA_richness_birmmals = DimArray(reshape(fill(0.0, 125*76), 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))

species_df_herps = species_df[:, 5:53]
species_df_herps.sum = [sum(row) for row in eachrow(species_df_herps)]

species_df_birmmals = species_df[:, 54:260]
species_df_birmmals.sum = [sum(row) for row in eachrow(species_df_birmmals)]

for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            DA_richness[j] =  species_df_matrix[i, 3]
            DA_richness_herps[j] =  species_df_herps[i, 50]
            DA_richness_birmmals[j] = species_df_birmmals[i, 208]
        # else
        #     DA_sum[j] = false
        end
    end
end
# serialize("DA_richness.jls", DA_richness)
# serialize("Objects/DA_richness_herps.jls", DA_richness_herps)
# serialize("Objects/DA_richness_birmmals.jls", DA_richness_birmmals)
DA_richness = deserialize("DA_richness.jls")::DimArray{Float64,2}
DA_richness_birmmals = deserialize("Objects/DA_richness_birmmals.jls")::DimArray{Float64,2}
DA_richness_herps = deserialize("Objects/DA_richness_herps.jls")::DimArray{Float64,2}

DA_sum_r = reverse(DA_sum, dims=1)
DA_sum_p = permutedims(DA_sum, (2, 1))
MK.plot(reverse(DA_sum_r, dims=1), colormap = :redsblues);

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

npp_DA = DimArray(zeros((125, 76)), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in 1:size(species_df, 1)
    # println(perro_cropped.Value[i])
    for j in 1:125*76
        # println(utmraster_da[j])
        if Float32(species_df.Value[i]) == Float32(utmraster_DA[j])
            npp_DA[j] = species_df_matrix[i, 261]
        end
    end
end
npp_DA_r = reverse(npp_DA, dims=1)
npp_DA_p = permutedims(npp_DA, (2, 1))
MK.plot(npp_DA_p, colormap = :darkrainbow);
########################## RULES #################################
################################################################## 
# Merge rule with non-negative constraint for MyStructs256
merge_rule = Cell{}() do data, state, I
    merged_state = state + 
        MyStructs256(SVector{256}(
            intrinsic_growth_256(state, 0.01, 200.0).a .+
            trophic_optimized(state, full_IM).a .+
            competition(state, competition_matrix).a
        ))
    return MyStructs256(max.(0.0, merged_state.a))
end
merge_rule_m = Cell{Tuple{:state, :npp}, :state}() do data, (state, npp), I
    
    merged_state = state + 
        MyStructs256(SVector{256}(
            intrinsic_growth_256(state, 0.01, 100.0).a .+
            trophic_optimized(state, full_IM).a
        ))
    return MyStructs256(max.(0.0, merged_state.a))
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

indisp = InwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.125),
    distancemethod=AreaToArea(30)
);

ruleset = Ruleset(indisp, boundary = Reflect())
full_IM = results[1][0.001]
# full_IM = full_IM_list[2][10]
# pepe = (grid1 = DA_with_abundances_p, grid2 = npp_DA_p)

gif_output = GifOutput(
    Matrix(DA_with_abundances_p); tspan = 1:100,
    filename = "custkernel1_sigma001.gif",
    fps = 20, 
    scheme=ColorSchemes.darkrainbow,
    zerocolor=RGB24(0.0),
    minval = 0.0, maxval = 20000.0,
    mask = Matrix(DA_sum_p)
)
@time saved = sim!(gif_output, ruleset)
sum(saved[1]).b
sum(saved[100]).b
#= If you get this error: ERROR: The constructor for MyStructs256{Float64}(::Float64, ::Float64)
is missing! is likely due to kernelproduct and Window{9/25}=#
DA_shity = deepcopy(DA_with_abundances_p)
for col in axes(DA_shity, 2), row in axes(DA_shity, 1) 
    DA_shity[row, col] = MyStructs256(SVector{256, Float64}(fill(100.0, 256)))
end

array_output = ArrayOutput(
    Matrix(DA_with_abundances_p); tspan = 1:1000,
    mask = Matrix(DA_sum_p)
)
array_output = ArrayOutput(
    pepe; tspan = 1:100,
    mask = Matrix(DA_sum)
)
@time r = sim!(array_output, Ruleset(indisp; boundary = Reflect()))
a = sum(r[1].state).b
b = sum(r[100].state).b
(a-b)/a
pepe = ( 
    state = Matrix(DA_with_abundances),
    npp = Matrix(npp_DA)
)

# MakieOutput
makie_output = MakieOutput(
    pepe, tspan = 1:100; 
    fps = 10, ruleset = Ruleset(merge_rule_m; boundary = Reflect()),
    mask = Matrix(DA_sum) 
) do (; layout, frame)
    ax1 = Axis(layout[1, 1])
    ax2 = Axis(layout[1, 2])
    Makie.heatmap!(ax1, frame.state, colormap = :inferno)
    Makie.heatmap!(ax2, frame.npp, colormap = :inferno) 
    ax1.yreversed[] = true
    ax2.yreversed[] = true
end

for i in 50:60
    println(r[2][i, 50])
end

max_values = MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
for i in 1:length(r)
    maxu = maximum(r[i])
    if maxu > max_values #&& maxu < 40000.0
        max_values = maxu
        println(max_values.b, " in time ", i)
    end
end

ruleset = Ruleset(Life())
# And the time-span for it to run
tspan = 1:100
# Create our own plots with Makie.jl
output = MakieOutput(rand(Bool, 200, 300); tspan, ruleset) do (; layout, frame)
    image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
end

init = fill(0.0, (100, 100))
init[2,2] = 1000000.0   
disp = InwardsDispersal(
    formulation=CustomKernel(0.125),
    distancemethod=AreaToArea(30)
)
rule = Ruleset(disp; boundary=Wrap())
out = ArrayOutput(init; tspan=1:100)
r = sim!(out, rule)
sum(r[1])
sum(r[100])

out = MakieOutput(init; tspan=1:100, ruleset) do (; layout, frame)
    image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
end    

init = fill(MyStructs256(SVector{256, Float64}(rand(256))), (100, 100))
init = fill(MyStructs256(SVector{256, Float64}(rand(256)*100)), (100,100))
out = ArrayOutput(init; tspan=1:100)
r = sim!(out, rule)
sum(r[1]).b
sum(r[100]).b

DA_trial = DimArray(reshape(rand(10.0:100.0, 100*100), 100, 100), (Dim{:a}(1:100), Dim{:b}(1:100)))
DA_mask = trues(dims(DA_trial))
for col in 1:20, row in 1:20
    DA_mask[row, col] = false
end
MK.plot(DA_mask, colormap = :redsblues);

indisp = InwardsDispersal(;
    formulation=ExponentialKernel(λ = 0.125)
)

out = ArrayOutput(DA_trial; tspan=1:100, mask = DA_mask) 
r = sim!(out, Ruleset(indisp; boundary=Wrap()))
sum(r[1])
sum(r[100])
### New NPP raster
bio_npp = Raster(joinpath(meta_path, "Rasters/CHELSA_npp_1981-2010_V.2.1.tif"))
bio_npp_croppped = crop(bio_npp, to = prec_raster)
# MK.plot(bio_npp_croppped, title = "")
reproject(bio_npp_croppped, prec_raster)
