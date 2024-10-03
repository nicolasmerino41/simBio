num_species = 256
include("One-click code.jl")

include("human_footprint.jl")

include("Time Series.jl")

#TODO At prevalence 0.277 or higher you get instability

caca = deepcopy(iberian_interact_NA)
self_regulation = 1.0
sigma = 1.0
epsilon = 1.0
full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon, self_regulation))
remove_variable(:alfa)
alfa = 0.2
exp(-(1^2) / (2*(alfa^2)))
m = maximum(npp_DA[.!isnan.(npp_DA)])
n = minimum(npp_DA[.!isnan.(npp_DA)])

######### HERE STARTS THE TIME SERIES PART #############
aux = (; combined_raster=combined_raster)
tspan = Date(2023,1):Month(1):Date(2025, 6)
pepe_rasterTS = (
    state = raster_with_abundances,
    k_raster = k_raster_hf_multiplicative,
    npp_raster = npp_raster,
    state_richness = raster_richness
)

self_regulation = 1.0
function int_Gr_for_timeseries(state::MyStructs256, self_regulation::AbstractFloat, combined_raster)
    return MyStructs256(SVector{256, Float64}(state.a .* self_regulation .+ combined_raster.*0.0))
end

timeseries_rule = let combined_raster_aux=Aux{:combined_raster}()
    Cell{:state, :state}() do data, state, I
        combined_raster = get(data, combined_raster_aux, I) # Important bit here!!
        merged_state = state + int_Gr_for_timeseries(state, self_regulation, combined_raster)
        return MyStructs256(SVector{256, Float64}(merged_state.a))
    end
end
k_rule = let combined_raster_aux=Aux{:combined_raster}()
    Cell{:k_raster, :k_raster}() do data, k_raster, I
        climate = get(data, combined_raster_aux, I) # Important bit here!!
        S_bio5 = 1 ./ (1 .+ abs.(climate .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        return MyStructs256(
            SVector{256, Float64}(
                S_bio5 .* herb_carv_svector
            )
        )   
    end
end

dynamic_biotic = let combined_raster_aux=Aux{:combined_raster}()
    Cell{Tuple{:state, :k_raster}, :state}() do data, (state, k_raster), I
        if any(isinf, state.a) || any(isnan, state.a)
            @error "state has NA values"
            println(I)
        end
        return MyStructs256(
            SVector{256, Float64}(
                max.(
                    0.0,
                    (state +
                    int_Gr_for_biotic_k(state, self_regulation, k_raster)  +
                    trophic_optimized(state, full_IM)).a
                )
            )
        )
    end
end
# Example usage of OutwardsDispersalRemix in a simulation
remix_outdisp = OutwardsDispersal{:state, :state}(
    formulation=CustomKernel(alfa),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
);
array_output = ArrayOutput(
    pepe_rasterTS; tspan = tspan,
    aux = aux,
    mask = raster_sum
)
@time s = sim!(array_output, Ruleset(k_rule, dynamic_biotic, remix_outdisp); boundary = Reflect())
MK.heatmap(combined_raster[Ti(2)])
f = Figure(resolution = (600, 400))
ax1 = Axis(f[1, 1])
ax2 = Axis(f[1, 2])
MK.heatmap!(ax1, s[1].state; colormap=custom_palette)
MK.heatmap!(ax2, s[end].state; colormap=custom_palette)
f
s[1].state == s[end].state

makie_output = MakieOutput(pepe_rasterTS, tspan = tspan;
    fps = 10, ruleset = Ruleset(k_rule, dynamic_biotic, remix_outdisp),
    aux = aux, mask = raster_sum) do (; layout, frame)

    plot_keys = [:biomass, :simulated_richness, :npp, :real_richness]
    titles = ["Biomass", "Simulated Richness", "NPP", "Real Richness"]

    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]
    println(typeof(frame[:k_raster]))

    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette)
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state]; interpolate=false, colormap=custom_palette)
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
        ax.yreversed[] = true
    end
end

Makie.heatmap(raster_with_abundances)