######################## COMPLEX RULES #############################
####################################################################
####################################################################
####################################################################
######################## DEFINING BASIC Herbivores METHODS ####################################
#################################################################################################
struct Herbivores{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{155, T}
    b::T

    # Custom constructor for automatic sum calculation
    function Herbivores(a::SVector{155, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    Herbivores(a::SVector{155, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for Herbivores
Base.zero(::Type{Herbivores{T}}) where {T <: AbstractFloat} = Herbivores(SVector{155, T}(fill(zero(T), 155)), zero(T))
Base.zero(x::Herbivores{T}) where {T <: AbstractFloat} = Herbivores(SVector{155, T}(fill(zero(T), 155)), zero(T))
Base.oneunit(::Type{Herbivores{T}}) where {T <: AbstractFloat} = MyStructs155(fill(oneunit(T), 155), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::Herbivores, y::Herbivores) = isless(x.b, y.b)
Base.isless(x::Herbivores, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::Herbivores, y::Herbivores) = Herbivores(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::Herbivores, y::Herbivores) = Herbivores(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::Herbivores, scalar::Real) = Herbivores(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::Herbivores, scalar::Real) = Herbivores(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::Herbivores, scalar::Real) = Herbivores(x.a .- scalar, x.b - scalar * 155)
Base.:+(x::Herbivores, scalar::Real) = Herbivores(x.a .+ scalar, x.b + scalar * 155)
Base.:*(x::Herbivores, y::AbstractVector) = Herbivores(x.a .* SVector{155, Float64}(y), sum(x.a .* SVector{155, Float64}(y)))
Base.:/(x::Herbivores, y::AbstractVector) = Herbivores(x.a ./ SVector{155, Float64}(y), sum(x.a ./ SVector{155, Float64}(y)))
Base.:*(y::AbstractVector, x::Herbivores) = Herbivores(SVector{155, Float64}(y) .* x.a, sum(SVector{155, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::Herbivores) = Herbivores(SVector{155, Float64}(y) ./ x.a, sum(SVector{155, Float64}(y) ./ x.a))

# Define what a NaN is for Herbivores
Base.isnan(x::Herbivores) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for Herbivores
function Base.sum(structs::Herbivores...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new Herbivores instance with the summed results
    return Herbivores(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for Herbivores
function Base.maximum(a::Herbivores, b::Herbivores)
    return Herbivores(max.(a.a, b.a))
end

# Define maximum for Herbivores with a scalar
function Base.maximum(a::Herbivores, b::AbstractFloat)
    return Herbivores(max.(a.a, b))
end

# Define maximum for a scalar with Herbivores
function Base.maximum(a::AbstractFloat, b::Herbivores)
    return Herbivores(max.(a, b.a))
end

# Define maximum for Herbivores
function Base.maximum(a::Herbivores)
    return maximum(a.a)
end

# Define maximum for a matrix of Herbivores
function Base.maximum(a::Matrix{Herbivores{AbstractFloat}})
    # Extract all `b` values from each Herbivores element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end

######################## DEFINING BASIC Predators METHODS ####################################
#################################################################################################
struct Predators{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{101, T}
    b::T

    # Custom constructor for automatic sum calculation
    function Predators(a::SVector{101, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    Predators(a::SVector{101, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for Predators
Base.zero(::Type{Predators{T}}) where {T <: AbstractFloat} = Predators(SVector{101, T}(fill(zero(T), 101)), zero(T))
Base.zero(x::Predators{T}) where {T <: AbstractFloat} = Predators(SVector{101, T}(fill(zero(T), 101)), zero(T))
Base.oneunit(::Type{Predators{T}}) where {T <: AbstractFloat} = MyStructs101(fill(oneunit(T), 101), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::Predators, y::Predators) = isless(x.b, y.b)
Base.isless(x::Predators, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::Predators, y::Predators) = Predators(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::Predators, y::Predators) = Predators(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::Predators, scalar::Real) = Predators(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::Predators, scalar::Real) = Predators(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::Predators, scalar::Real) = Predators(x.a .- scalar, x.b - scalar * 101)
Base.:+(x::Predators, scalar::Real) = Predators(x.a .+ scalar, x.b + scalar * 101)
Base.:*(x::Predators, y::AbstractVector) = Predators(x.a .* SVector{101, Float64}(y), sum(x.a .* SVector{101, Float64}(y)))
Base.:/(x::Predators, y::AbstractVector) = Predators(x.a ./ SVector{101, Float64}(y), sum(x.a ./ SVector{101, Float64}(y)))
Base.:*(y::AbstractVector, x::Predators) = Predators(SVector{101, Float64}(y) .* x.a, sum(SVector{101, Float64}(y) .* x.a))
Base.:/(y::AbstractVector, x::Predators) = Predators(SVector{101, Float64}(y) ./ x.a, sum(SVector{101, Float64}(y) ./ x.a))

# Define what a NaN is for Predators
Base.isnan(x::Predators) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for Predators
function Base.sum(structs::Predators...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new Predators instance with the summed results
    return Predators(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for Predators
function Base.maximum(a::Predators, b::Predators)
    return Predators(max.(a.a, b.a))
end

# Define maximum for Predators with a scalar
function Base.maximum(a::Predators, b::AbstractFloat)
    return Predators(max.(a.a, b))
end

# Define maximum for a scalar with Predators
function Base.maximum(a::AbstractFloat, b::Predators)
    return Predators(max.(a, b.a))
end

# Define maximum for Predators
function Base.maximum(a::Predators)
    return maximum(a.a)
end

# Define maximum for a matrix of Predators
function Base.maximum(a::Matrix{Predators{AbstractFloat}})
    # Extract all `b` values from each Predators element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end