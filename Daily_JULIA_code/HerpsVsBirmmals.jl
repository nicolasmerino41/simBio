using Pkg
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames, NaNMath
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources #, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase, JLD2, Random
using ColorSchemes, Colors #Crayons, 
using Makie, CairoMakie # ImageMagick
if ode_true
    using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Logging, ForwardDiff
end
# if end_true
#     using DifferentialEquations # EcologicalNetworksDynamics, 
# end
# using Plots
const DG, MK, AG, RS, Disp, DF = DynamicGrids, Makie, ArchGDAL, Rasters, Dispersal, DataFrames
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs256{T}) where {T <: AbstractFloat} = MyStructs256(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructsnum_species(fill(oneunit(T), num_species), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar * num_species)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar * num_species)
Base.:*(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a .* SVector{num_species, Float64}(y), sum(x.a .* SVector{num_species, Float64}(y)))
Base.:/(x::MyStructs256, y::AbstractVector) = MyStructs256(x.a ./ SVector{num_species, Float64}(y), sum(x.a ./ SVector{num_species, Float64}(y)))
Base.:*(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) .* x.a, sum(SVector{num_species, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::MyStructs256) = MyStructs256(SVector{num_species, Float64}(y) ./ x.a, sum(SVector{num_species, Float64}(y) ./ x.a))

# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyStructs256
function Base.sum(structs::MyStructs256...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs256 instance with the summed results
    return MyStructs256(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256, b::MyStructs256)
    return MyStructs256(max.(a.a, b.a))
end

# Define maximum for MyStructs256 with a scalar
function Base.maximum(a::MyStructs256, b::AbstractFloat)
    return MyStructs256(max.(a.a, b))
end

# Define maximum for a scalar with MyStructs256
function Base.maximum(a::AbstractFloat, b::MyStructs256)
    return MyStructs256(max.(a, b.a))
end

# Define maximum for MyStructs256
function Base.maximum(a::MyStructs256)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs256
function Base.maximum(a::Matrix{MyStructs256{AbstractFloat}})
    # Extract all `b` values from each MyStructs256 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

######################## DEFINING BASIC MYSTRUCTS254 METHODS ####################################
#################################################################################################
struct MyStructs254{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{254, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs254(a::SVector{254, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs254(a::SVector{254, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs254
Base.zero(::Type{MyStructs254{T}}) where {T <: AbstractFloat} = MyStructs254(SVector{254, T}(fill(zero(T), 254)), zero(T))
Base.zero(x::MyStructs254{T}) where {T <: AbstractFloat} = MyStructs254(SVector{254, T}(fill(zero(T), 254)), zero(T))
Base.oneunit(::Type{MyStructs254{T}}) where {T <: AbstractFloat} = MyStructs254(fill(oneunit(T), 254), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs254, y::MyStructs254) = isless(x.b, y.b)
Base.isless(x::MyStructs254, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs254, y::MyStructs254) = MyStructs254(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs254, y::MyStructs254) = MyStructs254(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs254, scalar::Real) = MyStructs254(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs254, scalar::Real) = MyStructs254(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs254, scalar::Real) = MyStructs254(x.a .- scalar, x.b - scalar * 254)
Base.:+(x::MyStructs254, scalar::Real) = MyStructs254(x.a .+ scalar, x.b + scalar * 254)
Base.:*(x::MyStructs254, y::AbstractVector) = MyStructs254(x.a .* SVector{254, Float64}(y), sum(x.a .* SVector{254, Float64}(y)))
Base.:/(x::MyStructs254, y::AbstractVector) = MyStructs254(x.a ./ SVector{254, Float64}(y), sum(x.a ./ SVector{254, Float64}(y)))
Base.:*(y::AbstractVector, x::MyStructs254) = MyStructs254(SVector{254, Float64}(y) .* x.a, sum(SVector{254, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::MyStructs254) = MyStructs254(SVector{254, Float64}(y) ./ x.a, sum(SVector{254, Float64}(y) ./ x.a))

# Define what a NaN is for MyStructs254
Base.isnan(x::MyStructs254) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyStructs254
function Base.sum(structs::MyStructs254...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs254 instance with the summed results
    return MyStructs254(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyStructs254
function Base.maximum(a::MyStructs254, b::MyStructs254)
    return MyStructs254(max.(a.a, b.a))
end

# Define maximum for MyStructs254 with a scalar
function Base.maximum(a::MyStructs254, b::AbstractFloat)
    return MyStructs254(max.(a.a, b))
end

# Define maximum for a scalar with MyStructs254
function Base.maximum(a::AbstractFloat, b::MyStructs254)
    return MyStructs254(max.(a, b.a))
end

# Define maximum for MyStructs254
function Base.maximum(a::MyStructs254)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs254
function Base.maximum(a::Matrix{MyStructs254{AbstractFloat}})
    # Extract all `b` values from each MyStructs254 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

################################## MYHERPS METHODS ###########################################
#################################################################################################
#################################################################################################
#################################################################################################
struct MyHerps{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{49, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyHerps(a::SVector{49, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyHerps(a::SVector{49, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyHerps
Base.zero(::Type{MyHerps{T}}) where {T <: AbstractFloat} = MyHerps(SVector{49, T}(fill(zero(T), 49)), zero(T))
Base.zero(x::MyHerps{T}) where {T <: AbstractFloat} = MyHerps(SVector{49, T}(fill(zero(T), 49)), zero(T))
Base.oneunit(::Type{MyHerps{T}}) where {T <: AbstractFloat} = MyHerps(fill(oneunit(T), 49), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyHerps, y::MyHerps) = isless(x.b, y.b)
Base.isless(x::MyHerps, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyHerps, y::MyHerps) = MyHerps(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyHerps, y::MyHerps) = MyHerps(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyHerps, scalar::Real) = MyHerps(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyHerps, scalar::Real) = MyHerps(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyHerps, scalar::Real) = MyHerps(x.a .- scalar, x.b - scalar * 49)
Base.:+(x::MyHerps, scalar::Real) = MyHerps(x.a .+ scalar, x.b + scalar * 49)

# Define what a NaN is for MyHerps
Base.isnan(x::MyHerps) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyHerps
function Base.sum(structs::MyHerps...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyHerps instance with the summed results
    return MyHerps(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyHerps
function Base.maximum(a::MyHerps, b::MyHerps)
    return MyHerps(max.(a.a, b.a))
end

# Define maximum for MyHerps with a scalar
function Base.maximum(a::MyHerps, b::AbstractFloat)
    return MyHerps(max.(a.a, b))
end

# Define maximum for a scalar with MyHerps
function Base.maximum(a::AbstractFloat, b::MyHerps)
    return MyHerps(max.(a, b.a))
end

# Define maximum for MyHerps
function Base.maximum(a::MyHerps)
    return maximum(a.a)
end

# Define maximum for a matrix of MyHerps
function Base.maximum(a::Matrix{MyHerps{AbstractFloat}})
    # Extract all `b` values from each MyHerps element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

################################## MYBIRMMALS METHODS ###########################################
#################################################################################################
#################################################################################################
#################################################################################################
struct MyBirmmals{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{205, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyBirmmals(a::SVector{205, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyBirmmals(a::SVector{205, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyBirmmals
Base.zero(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(SVector{205, T}(fill(zero(T), 205)), zero(T))
Base.zero(x::MyBirmmals{T}) where {T <: AbstractFloat} = MyBirmmals(SVector{205, T}(fill(zero(T), 205)), zero(T))
Base.oneunit(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(fill(oneunit(T), 205), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyBirmmals, y::MyBirmmals) = isless(x.b, y.b)
Base.isless(x::MyBirmmals, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .- scalar, x.b - scalar * 205)
Base.:+(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .+ scalar, x.b + scalar * 205)

# Define what a NaN is for MyBirmmals
Base.isnan(x::MyBirmmals) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for MyBirmmals
function Base.sum(structs::MyBirmmals...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyBirmmals instance with the summed results
    return MyBirmmals(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for MyBirmmals
function Base.maximum(a::MyBirmmals, b::MyBirmmals)
    return MyBirmmals(max.(a.a, b.a))
end

# Define maximum for MyBirmmals with a scalar
function Base.maximum(a::MyBirmmals, b::AbstractFloat)
    return MyBirmmals(max.(a.a, b))
end

# Define maximum for a scalar with MyBirmmals
function Base.maximum(a::AbstractFloat, b::MyBirmmals)
    return MyBirmmals(max.(a, b.a))
end

# Define maximum for MyBirmmals
function Base.maximum(a::MyBirmmals)
    return maximum(a.a)
end

# Define maximum for a matrix of MyBirmmals
function Base.maximum(a::Matrix{MyBirmmals{AbstractFloat}})
    # Extract all `b` values from each MyBirmmals element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

# Define zeros for all three types
function Base.zeros(dims::NTuple{2, Int}, type = nothing)
    if type == MyBirmmals{AbstractFloat}
        return [MyBirmmals(fill(0.0, 205)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyStructs256{AbstractFloat}
        return [MyStructs256(fill(0.0, 256)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyStructs254{AbstractFloat}
        return [MyStructs254(fill(0.0, 254)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyHerps{AbstractFloat}
        return [MyHerps(fill(0.0, 49)) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end

################# MYBIRMMALS KERNEL METHODS ################
###############################################################
###############################################################
# Define kernel product for MyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyBirmmals{AbstractFloat}}, kernel::SVector{9, AbstractFloat})
    
    result_a = SVector{205, AbstractFloat}(fill(0.0f0, 205))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end
    return MyBirmmals(result_a)
end

function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyBirmmals{AbstractFloat}}, kernel::SVector{25, AbstractFloat})
    
    result_a = SVector{205, AbstractFloat}(fill(0.0f0, 205))
    
    for (i, k) in enumerate(kernel)
        result_a += hood[i].a .* k
    end

    return MyBirmmals(result_a)
end

import Base.Broadcast: broadcastable

broadcastable(x::MyBirmmals) = Ref(x)
broadcastable(x::MyHerps) = Ref(x)

Base.:+(x::MyHerps, y::MyBirmmals) = MyStructs256(SVector{256, typeof(x.b)}(vcat(x.a, y.a)), x.b + y.b)
Base.:+(x::MyBirmmals, y::MyHerps) = MyStructs256(SVector{256, typeof(x.b)}(vcat(y.a, x.a)), y.b + x.b)