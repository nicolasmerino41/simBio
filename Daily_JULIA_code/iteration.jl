using Pkg

# Packages
using CSV, DataFrames
using Distributions, NamedArrays, StaticArrays
using DynamicGrids, Dispersal
using DimensionalData, Rasters, Serialization, ArchGDAL
using OrderedCollections, StatsBase

# Setup code
include("HerpsVsBirmmals.jl")
@time include("efficient_setup.jl")

epsilon = 1.0
sigma = 0.001
threshold = 0.2
alpha = 0.2
init_abund = 1.0

@time for prevalence in [0.1]

println("Break 1.jl")
self_regulation = 0.001

caca = deepcopy(iberian_interact_NA)
println("Break 2.jl")
full_IM = turn_adj_into_inter(caca, sigma, epsilon)
full_IM = Matrix(full_IM)

println("Break 3.jl") ########### RULES #############
outdisp = OutwardsDispersal{:birmmals, :birmmals}(;
    formulation=CustomKernel(alpha),
    distancemethod=AreaToArea(30),
    maskbehavior = Dispersal.CheckMaskEdges()
)
biotic_rule_k_herps = Cell{Tuple{:herps, :birmmals, :k_DA}, :herps}() do data, (herps, birmmals, k_DA), I
    if any(isinf, birmmals.a) || any(isnan, birmmals.a)
        @error "state has NA values in birmmals"
        println(I)
    end
    if any(isinf, herps.a) || any(isnan, herps.a)
        @error "state has NA values in herps"
        println(I)
    end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyHerps(max.(0.0, merged_state.a[1:49]))
end
biotic_rule_k_birmmals = Cell{Tuple{:herps, :birmmals, :k_DA}, :birmmals}() do data, (herps, birmmals, k_DA), I
    if typeof(birmmals) != MyBirmmals{Float64}
        @error "birmmals is not MyBirmmals"
    end
    if typeof(herps) != MyHerps{Float64}
        @error "herps is not MyHerps"
    end
    if any(isinf, birmmals.a) || any(isnan, birmmals.a)
        @error "state has NA values in birmmals"
        println(I)
    end
    if any(isinf, herps.a) || any(isnan, herps.a)
        @error "state has NA values in herps"
        println(I)
    end
    merged_state = deepcopy(herps) + deepcopy(birmmals) +
        int_Gr_for_biotic_k(deepcopy(herps) + deepcopy(birmmals), self_regulation, k_DA)  +
        trophic_optimized(deepcopy(herps) + deepcopy(birmmals), full_IM)
    return MyBirmmals(max.(0.0, merged_state.a[50:256]))
end

println("Break 4.jl")
function random_herp_dimarray(dimensions::Tuple{Int64, Int64}, prev)
    init_array = DimArray(zeros(dimensions, MyHerps{Float64}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
    for i in idx
        init_array[i] = MyHerps((fill(init_abund, 49)) .* sample([1,0], Weights([prev, 1-prev]), 49))
    end
    return init_array
end

println("Break 5.jl")
DA_random_herps_with_abundances = random_herp_dimarray(dimensions, prevalence)

println("Break 6.jl")
function random_birmmals_dimarray(dimensions::Tuple{Int64, Int64}, prev)
    init_array = DimArray(zeros(dimensions, MyBirmmals{Float64}), (Dim{:X}(1:dimensions[1]), Dim{:Y}(1:dimensions[2])))
    for i in idx
        init_array[i] = MyBirmmals((fill(init_abund, 207)) .* sample([1,0], Weights([prev, 1-prev]), 207))
    end
    return init_array
end

println("Break 7.jl")
DA_random_birmmals_with_abundances = random_birmmals_dimarray(dimensions, prevalence)

println("Break 8.jl")
#TODO At prevalence 0.277 or higher you get instability
pepe = (
    birmmals = Matrix(DA_random_birmmals_with_abundances),
    herps = Matrix(DA_random_herps_with_abundances),
    k_DA = Matrix(k_DA)
)

println("Break 9.jl")
##### LAX NICHE #####
array_output = ArrayOutput(
    pepe; tspan = 1:10,
    mask = Matrix(DA_sum)
)

println("Break 10.jl")
println("full_IM is a ", typeof(full_IM))
println("birmmals is a ", typeof(Matrix(DA_random_birmmals_with_abundances)))

@time p = sim!(array_output, Ruleset(biotic_rule_k_herps, biotic_rule_k_birmmals, outdisp))
println("Break 11.jl")

function richness_evaluation(array_output, DA_with_presences, thresh)
    matches = 0
    for i in idx 
        if any(isnan, array_output[i].a)
            throw(ArgumentError("NaN found in array_output"))
        end
        above_threshold = [x > thresh ? 1.0 : 0.0 for x in array_output[i].a]
        matches += sum(above_threshold .== DA_with_presences[i])
    end
    return matches/(length(idx)*256)
end
richness_eval = richness_evaluation(p[length(p)].birmmals + p[length(p)].herps, DA_with_presences, threshold)
println("Break 12.jl")

####################################################################
println("At prevalence: ", prevalence)
println("Epsilon: ", epsilon)
println("Sigma: ", sigma)
println("Alpha: ", alpha)
println("Threshold: ", threshold)
println("The richness_eval was: ", richness_eval)

# Saving the results
results = DataFrame(
    prevalence = prevalence,
    epsilon = epsilon,
    sigma = sigma,
    alpha = alpha,
    threshold = threshold,
    richness_eval = richness_eval
)

CSV.write("results/results.csv", results, append=true)

end