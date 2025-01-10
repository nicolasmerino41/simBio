using Pkg
cd(pwd())

dir = pwd()
# Packages

using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Dates, Distributions, Serialization
using DifferentialEquations, DiffEqCallbacks, LinearAlgebra

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

################################## MYBIRMMALS METHODS ###########################################
#################################################################################################
#################################################################################################
#################################################################################################
struct MyBirmmals{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{207, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyBirmmals(a::SVector{207, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyBirmmals(a::SVector{207, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyBirmmals
Base.zero(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(SVector{207, T}(fill(zero(T), 207)), zero(T))
Base.zero(x::MyBirmmals{T}) where {T <: AbstractFloat} = MyBirmmals(SVector{207, T}(fill(zero(T), 207)), zero(T))
Base.oneunit(::Type{MyBirmmals{T}}) where {T <: AbstractFloat} = MyBirmmals(fill(oneunit(T), 207), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyBirmmals, y::MyBirmmals) = isless(x.b, y.b)
Base.isless(x::MyBirmmals, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyBirmmals, y::MyBirmmals) = MyBirmmals(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .- scalar, x.b - scalar * 207)
Base.:+(x::MyBirmmals, scalar::Real) = MyBirmmals(x.a .+ scalar, x.b + scalar * 207)

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
        return [MyBirmmals(fill(0.0, 207)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyStructs256{AbstractFloat}
        return [MyStructs256(fill(0.0, 256)) for _ in 1:dims[1], _ in 1:dims[2]]
    elseif type == MyHerps{AbstractFloat}
        return [MyHerps(fill(0.0, 49)) for _ in 1:dims[1], _ in 1:dims[2]]
    else
        return [0.0 for _ in 1:dims[1], _ in 1:dims[2]]
    end
end

import Base.Broadcast: broadcastable

broadcastable(x::MyBirmmals) = Ref(x)
broadcastable(x::MyHerps) = Ref(x)

Base.:+(x::MyHerps, y::MyBirmmals) = MyStructs256(SVector{256, typeof(x.b)}(vcat(x.a, y.a)), x.b + y.b)
Base.:+(x::MyBirmmals, y::MyHerps) = MyStructs256(SVector{256, typeof(x.b)}(vcat(y.a, x.a)), y.b + x.b)