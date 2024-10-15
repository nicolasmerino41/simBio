num_species = 1150
######################## DEFINING BASIC MyStructs1150 METHODS ####################################
#################################################################################################

# Define the struct with 1150 elements
struct MyStructs1150{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs1150(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs1150(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero for MyStructs1150
Base.zero(::Type{MyStructs1150{T}}) where {T <: AbstractFloat} = MyStructs1150(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs1150{T}) where {T <: AbstractFloat} = MyStructs1150(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))

# Define oneunit for MyStructs1150
Base.oneunit(::Type{MyStructs1150{T}}) where {T <: AbstractFloat} = MyStructs1150(SVector{num_species, T}(fill(oneunit(T), num_species)), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs1150, y::MyStructs1150) = isless(x.b, y.b)
Base.isless(x::MyStructs1150, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly

# Addition of two MyStructs1150 instances
Base.:+(x::MyStructs1150, y::MyStructs1150) = MyStructs1150(x.a .+ y.a, x.b + y.b)

# Subtraction of two MyStructs1150 instances
Base.:-(x::MyStructs1150, y::MyStructs1150) = MyStructs1150(x.a .- y.a, x.b - y.b)

# Multiplication of MyStructs1150 by a scalar
Base.:*(x::MyStructs1150, scalar::Real) = MyStructs1150(x.a .* scalar, x.b * scalar)

# Division of MyStructs1150 by a scalar
Base.:/(x::MyStructs1150, scalar::Real) = MyStructs1150(x.a ./ scalar, x.b / scalar)

# Subtraction of a scalar from MyStructs1150
Base.:-(x::MyStructs1150, scalar::Real) = MyStructs1150(x.a .- scalar, x.b - scalar * num_species)

# Addition of a scalar to MyStructs1150
Base.:+(x::MyStructs1150, scalar::Real) = MyStructs1150(x.a .+ scalar, x.b + scalar * num_species)

# Element-wise multiplication with an AbstractVector
Base.:*(x::MyStructs1150, y::AbstractVector{T}) where {T <: Real} = MyStructs1150(x.a .* SVector{num_species, T}(y), sum(x.a .* SVector{num_species, T}(y)))
Base.:*(y::AbstractVector{T}, x::MyStructs1150) where {T <: Real} = MyStructs1150(SVector{num_species, T}(y) .* x.a, sum(SVector{num_species, T}(y) .* x.a))

# Element-wise division with an AbstractVector
Base.:/(x::MyStructs1150, y::AbstractVector{T}) where {T <: Real} = MyStructs1150(x.a ./ SVector{num_species, T}(y), sum(x.a ./ SVector{num_species, T}(y)))
Base.:/(y::AbstractVector{T}, x::MyStructs1150) where {T <: Real} = MyStructs1150(SVector{num_species, T}(y) ./ x.a, sum(SVector{num_species, T}(y) ./ x.a))

# Define isnan for MyStructs1150
Base.isnan(x::MyStructs1150) = isnan(x.b) || any(isnan, x.a)

# Define sum for MyStructs1150
function Base.sum(structs::MyStructs1150...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs1150 instance with the summed results
    return MyStructs1150(summed_a, summed_b)
end

# Define maximum for MyStructs1150 with another MyStructs1150
function Base.maximum(a::MyStructs1150, b::MyStructs1150)
    return MyStructs1150(max.(a.a, b.a), maximum(a.b, b.b))
end

# Define maximum for MyStructs1150 with a scalar
function Base.maximum(a::MyStructs1150, b::AbstractFloat)
    return MyStructs1150(max.(a.a, b), maximum(a.b, b))
end

# Define maximum for a scalar with MyStructs1150
function Base.maximum(a::AbstractFloat, b::MyStructs1150)
    return MyStructs1150(max.(a, b.a), maximum(a, b.b))
end

# Define maximum for MyStructs1150
function Base.maximum(a::MyStructs1150)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs1150
function Base.maximum(a::AbstractMatrix{MyStructs1150{T}}) where {T <: AbstractFloat}
    # Extract all `b` values from each MyStructs1150 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end
