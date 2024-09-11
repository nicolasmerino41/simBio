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
temp_raster = RS.Raster("utm_raster_temp.tif")
MK.plot(temp_raster);
prec_raster = RS.Raster("utm_raster_prec.tif")
MK.plot(prec_raster);

temp_DA = DimArray(temp_raster, (Dim{:a}(1:125), Dim{:b}(1:76)))
temp_DA = parent(temp_raster)
prec_DA = DimArray(prec_raster, (Dim{:a}(1:125), Dim{:b}(1:76)))
prec_DA = parent(prec_raster)
######## TEMPERATURE ##############
####### Species_temp
species_temp = CSV.File("DFs\\species_temp.csv") |> DataFrame
species_temp = species_temp[!, 2:4]
# Find the sorting indices based on matching names_order to species column
order_indices = indexin(names(iberian_interact_df), species_temp[:, :species])
# Reorder rows based on these indices
species_temp = species_temp[order_indices, :]
####### Species_temp_range
species_temp_range = CSV.File("DFs\\species_temp_range.csv") |> DataFrame
species_temp_range = species_temp_range[!, 2:4]
# Find the sorting indices based on matching names_order to species column
order_indices = indexin(names(iberian_interact_df), species_temp_range[:, :species])
# Reorder rows based on these indices
species_temp_range = species_temp_range[order_indices, :]
######## PRECIPITATION ##############
####### Species_prec
species_prec = CSV.File("DFs\\species_prec.csv") |> DataFrame
species_prec = species_prec[!, 2:4]
# Reorder rows based on these indices
species_prec = species_prec[order_indices, :]
####### Species_prec_range
species_prec_range = CSV.File("DFs\\species_prec_range.csv") |> DataFrame
species_prec_range = species_prec_range[!, 2:4]
# Reorder rows based on these indices
species_prec_range = species_prec_range[order_indices, :]

function int_Gr(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, temp::AbstractFloat, prec::AbstractFloat)
    return MyStructs256(self_regulation * (npp+0.1) .*
    SVector{256}((1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec))) .*
    SVector{256}((1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp))) .*
    state.a .* (1.0 - (state.b / ((npp+0.1)))))
end
function int_Gr_with_range(state::MyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, temp::AbstractFloat, prec::AbstractFloat)
    return MyStructs256(self_regulation * (npp+0.1) .*
    SVector{256}((1 ./ (1 .+ abs.(prec .- species_prec_range.mean_Prec) ./ species_prec_range.sd_Prec))) .*
    SVector{256}((1 ./ (1 .+ abs.(temp .- species_temp_range.mean_Temp) ./ species_temp_range.sd_Temp))) .*
    state.a .* (1.0 .- (state.a ./ ((npp+0.1)))))
end
dimensions = (125, 76)
idx_tupled = [(i, j) for i in 1:dimensions[1], j in 1:dimensions[2]]
function random_dimarray(dimensions::Tuple{Int64, Int64}; prevalence = 0.5)
    init_array = DimArray(zeros(dimensions, MyStructs256{Float64}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
    for i in idx
        init_array[i] = MyStructs256(100.0 .* SVector{256}(sample([1,0], Weights([prevalence, 1-prevalence]), 256)))
    end
    return init_array
end
DA_random_with_abundances = random_dimarray(dimensions; prevalence = 0.1)
# function int_Gr(state::mMyStructs256, self_regulation::AbstractFloat, npp::Union{AbstractFloat, AbstractArray}, temp::AbstractFloat, prec::AbstractFloat)
#     return mMyStructs256(self_regulation * npp .*
#     (1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
#     (1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
#     state.a .* (1.0 .- (state.a ./ (npp))))
# end
climatic_niche_rule = Cell{Tuple{:state, :npp, :temp, :prec}, :state}() do data, (state, temp, prec, npp), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @warn "state has NA values"
    # end
    # prec_factor = (1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec))
    # temp_factor = (1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp))
    # println(MyStructs256(SVector{256}((self_regulation * 1000.0) .* state.a .* (1.0 .- (state.a ./ 1000.0)))))  
    return state + int_Gr(state, self_regulation, npp, temp, prec) 
    # + MyStructs256(self_regulation * 100.0 .*
    #  SVector{256}(1 ./ (1 .+ abs.(prec .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
    #  SVector{256}(1 ./ (1 .+ abs.(temp .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
    #  state.a .* (1.0 .- (state.a ./ (100.0))))
end
function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end
# DA_with_abundances[18, 1] + MyStructs256(self_regulation * 100.0 .*
# SVector{256}(1 ./ (1 .+ abs.(30.0 .- species_prec.mean_Prec) ./ species_prec.sd_Prec)) .*
# SVector{256}(1 ./ (1 .+ abs.(20.0 .- species_temp.mean_Temp) ./ species_temp.sd_Temp)) .*
# DA_with_abundances[18, 1].a .* (1.0 .- (DA_with_abundances[18, 1].a ./ (100.0))))
indisp = InwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.1),
    distancemethod=AreaToArea(30)
);
outdisp = OutwardsDispersal{:state, :state}(;
    formulation=CustomKernel(0.1),
    distancemethod=AreaToArea(30),
    mask_flag = Dispersal.Mask()
);

pepe = ( 
    state = Matrix(DA_random_with_abundances),
    npp = Matrix(npp_DA),
    temp = temp_DA,
    prec = prec_DA
)
typeof(pepe[1]) == Matrix{MyStructs256{Float64}}
makie_output = MakieOutput(pepe, tspan = 1:1000; 
    fps = 50, ruleset = Ruleset(climatic_niche_rule, outdisp; boundary = Reflect()),
    mask = Matrix(DA_sum)) do (; layout, frame)
    
    # Setup the keys and titles for each plot
    plot_keys = [:state, :npp, :temp, :prec]  
    titles = ["Abundance", "Carrying capacity", "Temperature", "Precipitation"]
    
    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]
    
    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        Makie.heatmap!(ax, frame[key]; interpolate=false, colormap=:inferno)
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.yreversed[] = true
    end
end

array_output = ArrayOutput(
    pepe; tspan = 1:200,
    mask = Matrix(DA_sum)
)
@time r = sim!(array_output, Ruleset(climatic_niche_rule, outdisp; boundary = Reflect()))
richness_evaluation(r[200].state, DA_with_presences)
for i in 125 
    for j in 1:76
        if isone(Matrix(DA_sum)[i, j])
            println(i, " ", j)
        end
    end 
end
