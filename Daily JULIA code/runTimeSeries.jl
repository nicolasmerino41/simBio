PC = "nicol"
num_species = 256 
include("HerpsVsBirmmals.jl")
include("kernels.jl")
include("One-click code.jl")
include("human_footprint.jl")
include("Implicit competition for herbivores.jl")
tspan_intervals = Year(1)
static = false
include("2010-2100 RasterSeries.jl")

##### Parameters #####
self_regulation = 1.0
sigma = 0.5
sigma_comp = 0.1
epsilon = 1.0
beta = 3.0  
assymetry = 0.0
remove_variable(:alfa)
alfa = 0.1
##### Matrices #####
# Trophic
caca = deepcopy(iberian_interact_NA)
full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon, self_regulation, beta))
# Competition
competition_NA = deepcopy(iberian_interact_NA)
competition_NA .= 0.0
for i in names(competition_NA, 1), j in names(competition_NA, 2)
    if i in herbivore_names && j in herbivore_names
        competition_NA[i, j] = 1.0
    end
end
full_comp = turn_comp_into_inter(competition_NA, sigma_comp, assymetry)
##########################################################################
# exp(-(1^2) / (2*(alfa^2)))
m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

######### PREPARING THE MODEL ############
aux = (
    bio5_raster_series=bio5_raster_series,
    bio6_raster_series=bio6_raster_series,
    bio12_raster_series=bio12_raster_series
)

pepe_rasterTS = (
    state = raster_with_abundances,
    k_raster = k_raster_hf_multiplicative,
    npp_raster = npp_raster,
    state_richness = raster_richness
)

k_rule = let bio5_aux=Aux{:bio5_raster_series}(), bio6_aux=Aux{:bio6_raster_series}(), bio12_aux=Aux{:bio12_raster_series}()
    Cell{:k_raster, :k_raster}() do data, k_raster, I
        # bio5, bio6, bio12 = get(data, bio5_aux, I), get(data, bio6_aux, I), get(data, bio12_aux, I) # Important bit here!!
        bio5, bio6, bio12 = get(data, bio5_aux, I)[I[1], I[2]], get(data, bio6_aux, I)[I[1], I[2]], get(data, bio12_aux, I)[I[1], I[2]] 
        S_bio5 = 1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        S_bio6 = 1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)
        S_bio12 = 1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)
        return MyStructs256(
            SVector{256, Float64}(
                # S_bio5 .* S_bio6 .* S_bio12 .* herb_carv_svector
                S_bio5 .* herb_carv_svector
            )
        )   
    end
end
function GLV_raster(state::MyStructs256, k_DA::MyStructs256)
    return MyStructs256(
        SVector{256, Float64}(
            state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a) + ((full_comp * state.a) .* state.a)) 
        )
    )
end
biotic_GLV_raster = Cell{Tuple{:state, :k_raster}, :state}() do data, (state, k_raster), I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return MyStructs256(SVector{256, Float64}(max.(0.0, GLV_raster(state, k_raster).a)))
end
outdisp = OutwardsDispersal{:state, :state}(
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);
######## RUNNNING THE MODEL ############
array_output = ResultOutput(
    pepe_rasterTS; tspan = tspan,
    aux = aux,
    mask = raster_sum
)
@time s = sim!(array_output, Ruleset(k_rule, biotic_GLV_raster, outdisp); boundary = Reflect())

########### MAKIE ############
makie_output = MakieOutput(pepe_rasterTS, tspan = tspan;
    fps = 10, ruleset = Ruleset(k_rule, biotic_GLV_raster, outdisp),
    aux = aux, mask = raster_sum) do (; layout, frame)

    plot_keys = [:biomass, :simulated_richness, :npp, :real_richness]
    titles = ["Biomass", "Simulated Richness", "NPP", "Real Richness"]

    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]
    println(typeof(frame[:k_raster]))

    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
            ax.yreversed[] = true
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state], lambda_raster.multiplicative; colormap=custom_palette, colorrange = (0, 256))
            ax.yreversed[] = true
        elseif key == :npp
            Makie.heatmap!(ax, frame[:npp_raster]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
        elseif key == :real_richness
            Makie.heatmap!(ax, frame[:state_richness]; interpolate=false, colormap=custom_palette, colorrange = (0, 256))
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title
        ax.titlegap[] = 5
        ax.titlesize[] = 12
        ax.titlecolor[] = RGBA(0, 0, 0, 1)    
    end
end