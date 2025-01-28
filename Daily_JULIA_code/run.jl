begin
PC = "MM-1"
end_true = false
ode_true = true
num_species = 256
@time include("HerpsVsBirmmals.jl")
@time include("kernels.jl")
@time include("New_One-click code.jl")
@time include("HOfIS_intro.jl")
# @time include("human_footprint.jl")
# @time include("Implicit competition for herbivores.jl")
# include("2010-2100 RasterSeries.jl")
end
# pepe = (
#     birmmals = Matrix(DA_birmmals_with_abundances),
#     herps = Matrix(DA_herps_with_abundances),
#     k_DA = Matrix(k_DA_hf_additive),
#     birmmals_richness = Matrix(DA_richness_birmmals),
#     herps_richness = Matrix(DA_richness_herps)
# )

# DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)
# pepe_state = (
#     state = Matrix(DA_with_abundances),
#     k_DA = Matrix(k_DA_hf_additive),
#     npp_DA = Matrix(npp_DA),
#     state_richness = Matrix(DA_richness)
# )

# ##### Parameters #####
# self_regulation = 1.0
# sigma = 0.0001
# sigma_comp = 0.0
# epsilon = 0.1
# beta = 1.0  
# assymetry = 1.0
# remove_variable(:alfa)
# alfa = 0.9
# ##### Matrices #####
# # Trophic
# caca = deepcopy(iberian_interact_NA)
# full_IM = Matrix(turn_adj_into_inter(caca, sigma, epsilon, self_regulation, beta))
# count = 0
# for i in 1:256
#     if all(full_IM[:, i] == 0.0)
#         count += 1
#     end
# end
# # Competition
# competition_NA = deepcopy(iberian_interact_NA)
# competition_NA .= 0.0
# for i in names(competition_NA, 1), j in names(competition_NA, 2)
#     if i in herbivore_names && j in herbivore_names
#         competition_NA[i, j] = 1.0
#     end
# end
# full_comp = turn_comp_into_inter(competition_NA, sigma_comp, assymetry)
# ##########################################################################
# # exp(-(1^2) / (2*(alfa^2)))
# m = maximum(npp_DA[.!isnan.(npp_DA)])
# n = minimum(npp_DA[.!isnan.(npp_DA)])

# ##### RULES #####
# function GLV(state::MyStructs256, k_DA::MyStructs256)
#     return MyStructs256(
#         SVector{256, Float64}(
#             state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a) + ((full_comp * state.a) .* state.a)) 
#         )
#     )
# end

# biotic_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
#     # if any(isinf, state.a) || any(isnan, state.a)
#     #     @error "state has NA values"
#     #     println(I)
#     # end
#     return MyStructs256(SVector{256, Float64}(max.(0.0, GLV(state, k_DA).a)))
# end

# catastrophic_rule = let prob_chaos = 0.1
#     Neighbors{:cata_layer, :cata_layer}(Moore(1)) do data, neighborhood, cata_layer, I
#         event = rand() <= prob_chaos ? true : false
#         if event & isone(cata_layer)
#             return false   
#         elseif event & !iszero(cata_layer)
#             return true
#         else
#             return cata_layer
#         end
#     end
# end

# neighbors_rule = let prob_chaos = 0.01
#     Neighbors{:state, :state}(Moore(1)) do data, neighborhood, state, I
#         # println(typeof(DG.mask(data)[I[1], I[2]]))
#         # println(I)
#         event = rand() <= prob_chaos ? true : false
#         if event
#             return MyStructs256(SVector{256, Float64}(fill(0.0, 256)))
#         else
#             return state
#         # elseif any(isone, neighborhood)
#             # return state - MyStructs256(SVector{256, Float64}(fill(0.1, 256)))
#         end
#     end
# end

# ######## DISPERSAL ########
# scaling_vector = fill(1.0, 49)
# scaling_vector = append!(scaling_vector, fill(1.0, 207))
# remix_outdisp = OutwardsDispersalRemix{:state, :state}(
#     formulation=CustomKernel(alfa),
#     distancemethod=AreaToArea(30),
#     maskbehavior = Dispersal.CheckMaskEdges()
# );

# outdisp = OutwardsDispersal{:state, :state}(;
#     formulation=CustomKernel(alfa),
#     distancemethod=AreaToArea(30),
#     maskbehavior = Dispersal.CheckMaskEdges()
# );

# indisp = InwardsDispersal{:state, :state}(;
#     formulation=ExponentialKernel(),
#     distancemethod=AreaToArea(30)
# );

# ##### MAKIE STATE #####
# array_output = ResultOutput(
#     pepe_state; tspan = 1:12,
#     mask = Matrix(DA_sum)
# )

# @time c = sim!(array_output, Ruleset(biotic_GLV, outdisp; boundary = Reflect(), proc = ThreadedCPU()))

# makie_output = MakieOutput(pepe_state, tspan = 1:200;
#     fps = 10, ruleset = Ruleset(biotic_GLV, outdisp; boundary = Reflect()),
#     mask = Matrix(DA_sum)) do (; layout, frame)

#     # Setup the keys and titles for each plot
#     plot_keys = [:biomass, :simulated_richness, :npp, :real_richness]
#     titles = ["Biomass", "Simulated Richness", "NPP", "Real Richness"]

#     # Create axes for each plot and customize them
#     axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]

#     # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
#     for (ax, key, title) in zip(axes, plot_keys, titles)
#         if key == :biomass
#             Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
#         # elseif key == :simulated_richness
#         #     Makie.image!(ax, frame[:state], lambda_DA.multiplicative; colormap=custom_palette, colorrange = (0, 256))
#         elseif key == :npp
#             Makie.heatmap!(ax, frame[:npp_DA]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
#         elseif key == :real_richness
#             Makie.heatmap!(ax, frame[:state_richness]; interpolate=false, colormap=custom_palette, colorrange = (0, 256))
#         end
#         hidexdecorations!(ax; grid=false)
#         hideydecorations!(ax; grid=false)
#         ax.title = title  # Set the title for each axis
#         ax.titlegap[] = 5  # Adjust the title gap to make it smaller
#         ax.titlesize[] = 12 # Set the title font size
#         ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
#         ax.yreversed[] = true
#     end
# end

# #########################################################################
# ####################### MAKIE BIRMMALS/HERPS ############################
# #########################################################################
# makie_output = MakieOutput(pepe, tspan = 1:1000;
#     fps = 10, ruleset = Ruleset(biotic_rule_k_herps, biotic_rule_k_birmmals, outdisp_birmmals, outdisp_herps;
#     boundary = Reflect()),
#     mask = Matrix(DA_sum)) do (; layout, frame)

#     # Setup the keys and titles for each plot
#     plot_keys = [:birmmals_biomass, :herps_biomass, :birmmals_richness, :herps_richness, :real_birmmals_richness, :real_herps_richness]
#     titles = ["Birds+Mammals Biomass", "Herps Biomass", "Simulated Birds+Mammals Richness", "Simulated Herps Richness", "Real Birmmals Richness", "Real Herps Richness"]

#     # Create axes for each plot and customize them
#     axes = [Axis(layout[i, j]; title=titles[(i-1)*3 + j]) for i in 1:2, j in 1:3]

#     # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
#     for (ax, key, title) in zip(axes, plot_keys, titles)
#         if key == :birmmals_biomass
#             Makie.heatmap!(ax, frame[:birmmals]; interpolate=false, colormap=custom_palette)
#         elseif key == :herps_biomass
#             Makie.heatmap!(ax, frame[:herps]; interpolate=false, colormap=custom_palette)
#         elseif key == :birmmals_richness
#             Makie.image!(ax, frame[:birmmals]; colormap=custom_palette)
#         elseif key == :herps_richness
#             Makie.image!(ax, frame[:herps]; colormap=custom_palette)
#         elseif key == :real_birmmals_richness
#             Makie.heatmap!(ax, frame[:birmmals_richness]; interpolate=false, colormap=custom_palette)
#         elseif key == :real_herps_richness
#             Makie.heatmap!(ax, frame[:herps_richness]; interpolate=false, colormap=custom_palette)
#         end
#         hidexdecorations!(ax; grid=false)
#         hideydecorations!(ax; grid=false)
#         ax.title = title  # Set the title for each axis
#         ax.titlegap[] = 5  # Adjust the title gap to make it smaller
#         ax.titlesize[] = 12 # Set the title font size
#         ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
#         ax.yreversed[] = true
#     end
# end
