PC = "nicol"
num_species = 256
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\HerpsVsBirmmals.jl"))
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\kernels.jl"))
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\One-click code.jl"))
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\human_footprint.jl"))
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio\\Daily JULIA code\\Implicit competition for herbivores.jl"))
# include("2010-2100 RasterSeries.jl")
# pepe = (
#     birmmals = Matrix(DA_birmmals_with_abundances),
#     herps = Matrix(DA_herps_with_abundances),
#     k_DA = Matrix(k_DA_hf_additive),
#     birmmals_richness = Matrix(DA_richness_birmmals),
#     herps_richness = Matrix(DA_richness_herps)
# )
DA_with_abundances = deepcopy(DA_birmmals_with_abundances) + deepcopy(DA_herps_with_abundances)
pepe_state = (
    state = Matrix(DA_with_abundances),
    k_DA = Matrix(k_DA_hf_additive),
    npp_DA = Matrix(npp_DA),
    state_richness = Matrix(DA_richness)
)

##### Parameters #####
self_regulation = 1.0
sigma = 0.01
sigma_comp = 0.5
epsilon = 1.0
beta = 1.0  
assymetry = 1.0
remove_variable(:alfa)
alfa = 0.5
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

carnivoress = deepcopy(herb_carv_vector)
herbivores = convert(Vector{Bool}, carnivoress .== 1.0)
carnivores = convert(Vector{Bool}, carnivoress .!= 1.0)  
##### RULES #####
# Function 1: Competition
function competition(state::MyStructs256, k_DA::MyStructs256, full_comp::AbstractMatrix, carnivores)
    competition_effect = state.a .* (k_DA.a - state.a) .+ (full_comp * state.a) .* state.a
    competition_scaled = competition_effect .* herbivores  # Zero for predators
    
    return MyStructs256(SVector{256, Float64}(state.a .+ competition_scaled))
end

# Function 2: Scaling Up
function scaling_up(state::MyStructs256, lambda_DA::AbstractFloat, carnivores)
    upscaling_vector = ifelse.(carnivores, 1.0, lambda_DA)
    scaled_a = state.a .* SVector{256, Float64}(upscaling_vector)
    return MyStructs256(scaled_a)
end

# Function 3: Predation
function predation(state::MyStructs256, full_IM::AbstractArray)
    pred_effect = (full_IM * state.a) .* state.a
    new_a = state.a .+ pred_effect
    new_a = max.(0.0, new_a)  # Ensure no negative abundances
    return MyStructs256(new_a)
end

# Function 4: Scaling Down
function scaling_down(state::MyStructs256, lambda_DA::AbstractFloat, carnivores)
    downscaling_vector = ifelse.(carnivores, 1.0, 1.0 / lambda_DA)
    scaled_a = state.a .* SVector{256, Float64}(downscaling_vector)
    return MyStructs256(scaled_a)
end

new_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
    st = competition(state, k_DA, full_comp, carnivores)
    su = scaling_up(st, lambda_DA.multiplicative[I[1], I[2]], carnivores)
    pr = predation(su, full_IM)
    sd = scaling_down(pr, lambda_DA.multiplicative[I[1], I[2]], carnivores)
    return sd
end

##### MAKIE STATE #####
array_output = ResultOutput(
    pepe_state; tspan = 1:12,
    mask = Matrix(DA_sum)
)

@time c = sim!(array_output, Ruleset(new_GLV))

function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2}, new::Bool = false)
    println("new: ", new)
    scalars = map((mystruct, lambda_value) -> MyStructs256(SVector{256, Float64}((mystruct.a .* ifelse.(carnivores, 1.0, lambda_value))))b, A, lambda_grid)
    return Makie.convert_arguments(t, scalars)
end

makie_output = MakieOutput(pepe_state, tspan = 1:200;
    fps = 10, ruleset = Ruleset(new_GLV),
    mask = Matrix(DA_sum)) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:biomass, :simulated_richness, :npp, :real_richness]
    titles = ["Biomass", "Simulated Richness", "NPP", "Real Richness"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:2, j in 1:2]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state], lambda_DA.multiplicative, true; interpolate=false, colormap=custom_palette, colorrange = (0, m))
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state], lambda_DA.multiplicative; colormap=custom_palette, colorrange = (0, 256))
        elseif key == :npp
            Makie.heatmap!(ax, frame[:npp_DA]; interpolate=false, colormap=custom_palette, colorrange = (0, m))
        elseif key == :real_richness
            Makie.heatmap!(ax, frame[:state_richness]; interpolate=false, colormap=custom_palette, colorrange = (0, 256))
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.titlegap[] = 5  # Adjust the title gap to make it smaller
        ax.titlesize[] = 12 # Set the title font size
        ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
        ax.yreversed[] = true
    end
end

