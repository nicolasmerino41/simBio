using Pkg
cd(pwd())

dir = pwd()
# Packages

using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Dates, Distributions, Serialization
using DifferentialEquations, DiffEqCallbacks, LinearAlgebra
using DimensionalData, Rasters
using Random, Logging
using Base.Threads

######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
#################################################################################################
num_species = 254
struct MyStructs254{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs254(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs254(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs254
Base.zero(::Type{MyStructs254{T}}) where {T <: AbstractFloat} = MyStructs254(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs254{T}) where {T <: AbstractFloat} = MyStructs254(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.oneunit(::Type{MyStructs254{T}}) where {T <: AbstractFloat} = MyStructsnum_species(fill(oneunit(T), num_species), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs254, y::MyStructs254) = isless(x.b, y.b)
Base.isless(x::MyStructs254, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs254, y::MyStructs254) = MyStructs254(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs254, y::MyStructs254) = MyStructs254(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs254, scalar::Real) = MyStructs254(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs254, scalar::Real) = MyStructs254(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs254, scalar::Real) = MyStructs254(x.a .- scalar, x.b - scalar * num_species)
Base.:+(x::MyStructs254, scalar::Real) = MyStructs254(x.a .+ scalar, x.b + scalar * num_species)
Base.:*(x::MyStructs254, y::AbstractVector) = MyStructs254(x.a .* SVector{num_species, Float64}(y), sum(x.a .* SVector{num_species, Float64}(y)))
Base.:/(x::MyStructs254, y::AbstractVector) = MyStructs254(x.a ./ SVector{num_species, Float64}(y), sum(x.a ./ SVector{num_species, Float64}(y)))
Base.:*(y::AbstractVector, x::MyStructs254) = MyStructs254(SVector{num_species, Float64}(y) .* x.a, sum(SVector{num_species, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::MyStructs254) = MyStructs254(SVector{num_species, Float64}(y) ./ x.a, sum(SVector{num_species, Float64}(y) ./ x.a))

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
    elseif type == MyStructs254{AbstractFloat}
        return [MyStructs254(fill(0.0, 254)) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end
