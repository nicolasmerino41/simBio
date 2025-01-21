herb_carv_vector
spain_names
herbivore_names = []
for i in axes(iberian_interact_NA, 1)
    if all(x -> x == 0, iberian_interact_NA[i, :])
        push!(herbivore_names, names(iberian_interact_NA, 1)[i])
    end
end
predator_names = setdiff(spain_names, herbivore_names)
num_herbivores = length(herbivore_names)
@time include("HerbsVsPreds.jl")

if false
# Get a map of herbivore abundance
DA_herbivores = DimArray(reshape([Herbivores(SVector{num_herbivores, Float64}(fill(0.0, num_herbivores))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in idx
    v = Vector{Float64}()
    for j in 1:256
        if herb_carv_vector[j] == 1.0 && !iszero(DA_with_abundances[i].a[j])
            push!(v, 1.0)
        elseif herb_carv_vector[j] == 1.0 && iszero(DA_with_abundances[i].a[j])
            push!(v, 0.0)
        end
    end
    new = Herbivores(SVector{num_herbivores, Float64}(v))
    DA_herbivores[i] = new
end 
# Get a map of predator abundance
DA_predators = DimArray(reshape([Predators(SVector{101, Float64}(fill(0.0, 101))) for _ in 1:125*76], 125, 76), (Dim{:a}(1:125), Dim{:b}(1:76)))
for i in idx
    v = Vector{Float64}()
    for j in 1:256
        if !isone(herb_carv_vector[j]) && !iszero(DA_with_abundances[i].a[j])
            push!(v, 1.0)
        elseif !isone(herb_carv_vector[j]) && iszero(DA_with_abundances[i].a[j])
            push!(v, 0.0)
        end
    end
    new = Predators(SVector{101, Float64}(v))
    DA_predators[i] = new
end
end

function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:Herbivores, 2})
    scalars = map(herbivores -> herbivores.b, A)
    return Makie.convert_arguments(t, scalars)
end

function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:Predators, 2})
    scalars = map(predators -> predators.b, A)
    return Makie.convert_arguments(t, scalars)
end

if false
begin
    fig = Figure(resolution = (500,300))

    ax = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    MK.heatmap!(
        ax,
        Matrix(DA_herbivores),
        interpolate=false, 
        colormap=custom_palette,
        colorrange = (0, num_herbivores)
    )

    MK.heatmap!(
        ax2,
        Matrix(DA_predators),
        interpolate=false, 
        colormap=custom_palette,
        colorrange = (0, 101)
    )

    ax.yreversed[] = true
    ax2.yreversed[] = true

    display(fig)
end
end
############# PREDATOR-TO-PREY DICTIONARY ############################
# Create the predator-to-prey dictionary
using DataStructures

# Create an ordered dictionary
predator_prey_dict = OrderedDict{String, Vector{String}}()

# Iterate over all species in `spain_names`
for i in 1:length(spain_names)
    if herb_carv_vector[i] == 1e-8
        prey_indices = findall(x -> x == 1, iberian_interact_NA[i, :])
        prey_names = spain_names[prey_indices]
        predator_prey_dict[spain_names[i]] = prey_names
    end
end

# Print the ordered dictionary
# println(predator_prey_dict)

################## PREDATOR NUMBER OF PREY ############################
# Create the predator-to-prey dictionary with counts of total prey, herbivores, and predators
predator_prey_count = OrderedDict{String, Tuple{Int, Int, Int}}()

# Iterate over all species
for i in 1:length(spain_names)
    # Check if the current species is a predator
    if herb_carv_vector[i] == 1e-8
        # Initialize counts
        total_prey = 0
        herbivore_prey = 0
        predator_prey = 0

        # Iterate over the prey of this predator (row i of iberian_interact_NA)
        for j in 1:length(spain_names)
            if iberian_interact_NA[i, j] == 1
                total_prey += 1
                if herb_carv_vector[j] == 1.0
                    herbivore_prey += 1
                elseif herb_carv_vector[j] == 1e-8
                    predator_prey += 1
                end
            end
        end

        # Add predator and its prey details to the dictionary
        predator_prey_count[spain_names[i]] = (total_prey, herbivore_prey, predator_prey)
    end
end