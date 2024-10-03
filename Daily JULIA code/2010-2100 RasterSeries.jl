bounds = (X(Rasters.Between(-10, 4)), Y(Rasters.Between(35.5, 44)))

# Step 2: Load your initial and final rasters (replace with your actual data paths)
bio5_initial_raster = Raster("CHELSA/BioClim/bio5_2010.tif", lazy = true)[bounds...]
bio5_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio5_2100.tif", lazy = true)[bounds...]

bio6_initial_raster = Raster("CHELSA/BioClim/bio6_2010.tif", lazy = true)[bounds...]
bio6_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio6_2100.tif", lazy = true)[bounds...]

bio12_initial_raster = Raster("CHELSA/BioClim/bio12_2010.tif", lazy = true)[bounds...]
bio12_final_raster = Raster("CHELSA/Future/Climate/SSP370/IPSLCM6ALR/bio12_2100.tif", lazy = true)[bounds...]

# Step 3: Resample all rasters to match utmraster
bio5_initial_raster = Rasters.resample(bio5_initial_raster, to = utmraster)
bio5_final_raster = Rasters.resample(bio5_final_raster, to = utmraster)

bio6_initial_raster = Rasters.resample(bio6_initial_raster, to = utmraster)
bio6_final_raster = Rasters.resample(bio6_final_raster, to = utmraster)

bio12_initial_raster = Rasters.resample(bio12_initial_raster, to = utmraster)
bio12_final_raster = Rasters.resample(bio12_final_raster, to = utmraster)

# Step 4: Apply the scaling factors and offsets
# For bio5 and bio6: Temperature in °C = (Value * Scale) + Offset
# Scale = 0.1, Offset = -273.15

bio5_initial_raster = map(x -> (Float32(x) * 0.1) + (-273.15), bio5_initial_raster)
bio5_final_raster = map(x -> (Float32(x) * 0.1) + (-273.15), bio5_final_raster)

bio6_initial_raster = Rasters.map(x -> (Float32(x) * 0.1) + (-273.15), bio6_initial_raster)
bio6_final_raster = Rasters.map(x -> (Float32(x) * 0.1) + (-273.15), bio6_final_raster)

# For bio12: Precipitation in mm/year = (Value * Scale) + Offset
# Scale = 0.1, Offset = 0

bio12_initial_raster = Rasters.map(x -> Float32(x) * 0.1, bio12_initial_raster)
bio12_final_raster = Rasters.map(x -> Float32(x) * 0.1, bio12_final_raster)

# Step 5: Define the time span and dates
t_start = Date(2010, 1, 1)
t_end = Date(2100, 1, 1)
tspan = t_start:Month(1):t_end

# Step 6: Calculate the interpolation factors
total_days = Dates.value(t_end - t_start)
factors = [Dates.value(d - t_start) / total_days for d in tspan]

# Step 7: Generate the interpolated rasters for each variable
bio5_raster_list = Raster[]
bio6_raster_list = Raster[]
bio12_raster_list = Raster[]

for f in factors
    # Interpolate bio5 raster
    bio5_interp = (1 - f) .* bio5_initial_raster + f .* bio5_final_raster
    push!(bio5_raster_list, bio5_interp)
    
    # Interpolate bio6 raster
    bio6_interp = (1 - f) .* bio6_initial_raster + f .* bio6_final_raster
    push!(bio6_raster_list, bio6_interp)
    
    # Interpolate bio12 raster
    bio12_interp = (1 - f) .* bio12_initial_raster + f .* bio12_final_raster
    push!(bio12_raster_list, bio12_interp)
end

# Step 8: Create RasterSeries for each variable
bio5_raster_series = RasterSeries(bio5_raster_list, Ti(tspan))
bio6_raster_series = RasterSeries(bio6_raster_list, Ti(tspan))
bio12_raster_series = RasterSeries(bio12_raster_list, Ti(tspan))



# Create a figure with 2 axes
fig = Figure(resolution = (800, 400))

# Plot the first raster in the bio5 series
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
heatmap!(ax1, bio5_raster_series[Ti(1)])
heatmap!(ax2, bio5_raster_series[Ti(end)])

# Display the figure
fig

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

k_rule = let bio5_aux=Aux{:bio5_raster_series}()
    Cell{:k_raster, :k_raster}() do data, k_raster, I
        # bio5, bio6, bio12 = get(data, bio5_aux, I), get(data, bio6_aux, I), get(data, bio12_aux, I) # Important bit here!!
        bio5 = get(data, bio5_aux, I)[I[1], I[2]] 
        S_bio5 = 1 ./ (1 .+ abs.(bio5 .- lax_species_niches.mean_bio5) ./ lax_species_niches.sd_bio5)
        # S_bio6 = 1 ./ (1 .+ abs.(bio6 .- lax_species_niches.mean_bio6) ./ lax_species_niches.sd_bio6)
        # S_bio12 = 1 ./ (1 .+ abs.(bio12 .- lax_species_niches.mean_bio12) ./ lax_species_niches.sd_bio12)
        return MyStructs256(
            SVector{256, Float64}(
                # S_bio5 .* S_bio6 .* S_bio12 .* herb_carv_svector
                S_bio5 .* herb_carv_svector
            )
        )   
    end
end
function GLV_raster(state::MyStructs256, k_raster::MyStructs256)
    return MyStructs256(
        SVector{256, Float64}(
            state.a + (state.a .* (k_raster.a - state.a) + ((full_IM * state.a) .* state.a)) 
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