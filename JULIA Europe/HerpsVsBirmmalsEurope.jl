num_species = 1149
######################## DEFINING BASIC MyStructs1149 METHODS ####################################
#################################################################################################

# Define the struct with 1149 elements
struct MyStructs1149{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{num_species, T}
    b::T

    # Custom constructor for automatic sum calculation
    function MyStructs1149(a::SVector{num_species, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs1149(a::SVector{num_species, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero for MyStructs1149
Base.zero(::Type{MyStructs1149{T}}) where {T <: AbstractFloat} = MyStructs1149(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))
Base.zero(x::MyStructs1149{T}) where {T <: AbstractFloat} = MyStructs1149(SVector{num_species, T}(fill(zero(T), num_species)), zero(T))

# Define oneunit for MyStructs1149
Base.oneunit(::Type{MyStructs1149{T}}) where {T <: AbstractFloat} = MyStructs1149(SVector{num_species, T}(fill(oneunit(T), num_species)), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::MyStructs1149, y::MyStructs1149) = isless(x.b, y.b)
Base.isless(x::MyStructs1149, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly

# Addition of two MyStructs1149 instances
Base.:+(x::MyStructs1149, y::MyStructs1149) = MyStructs1149(x.a .+ y.a, x.b + y.b)

# Subtraction of two MyStructs1149 instances
Base.:-(x::MyStructs1149, y::MyStructs1149) = MyStructs1149(x.a .- y.a, x.b - y.b)

# Multiplication of MyStructs1149 by a scalar
Base.:*(x::MyStructs1149, scalar::Real) = MyStructs1149(x.a .* scalar, x.b * scalar)

# Division of MyStructs1149 by a scalar
Base.:/(x::MyStructs1149, scalar::Real) = MyStructs1149(x.a ./ scalar, x.b / scalar)

# Subtraction of a scalar from MyStructs1149
Base.:-(x::MyStructs1149, scalar::Real) = MyStructs1149(x.a .- scalar, x.b - scalar * num_species)

# Addition of a scalar to MyStructs1149
Base.:+(x::MyStructs1149, scalar::Real) = MyStructs1149(x.a .+ scalar, x.b + scalar * num_species)

# Element-wise multiplication with an AbstractVector
Base.:*(x::MyStructs1149, y::AbstractVector{T}) where {T <: Real} = MyStructs1149(x.a .* SVector{num_species, T}(y), sum(x.a .* SVector{num_species, T}(y)))
Base.:*(y::AbstractVector{T}, x::MyStructs1149) where {T <: Real} = MyStructs1149(SVector{num_species, T}(y) .* x.a, sum(SVector{num_species, T}(y) .* x.a))

# Element-wise division with an AbstractVector
Base.:/(x::MyStructs1149, y::AbstractVector{T}) where {T <: Real} = MyStructs1149(x.a ./ SVector{num_species, T}(y), sum(x.a ./ SVector{num_species, T}(y)))
Base.:/(y::AbstractVector{T}, x::MyStructs1149) where {T <: Real} = MyStructs1149(SVector{num_species, T}(y) ./ x.a, sum(SVector{num_species, T}(y) ./ x.a))

# Define isnan for MyStructs1149
Base.isnan(x::MyStructs1149) = isnan(x.b) || any(isnan, x.a)

# Define sum for MyStructs1149
function Base.sum(structs::MyStructs1149...)
    # Sum the 'a' vectors
    summed_a = sum([s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new MyStructs1149 instance with the summed results
    return MyStructs1149(summed_a, summed_b)
end

# Define maximum for MyStructs1149 with another MyStructs1149
function Base.maximum(a::MyStructs1149, b::MyStructs1149)
    return MyStructs1149(max.(a.a, b.a), maximum(a.b, b.b))
end

# Define maximum for MyStructs1149 with a scalar
function Base.maximum(a::MyStructs1149, b::AbstractFloat)
    return MyStructs1149(max.(a.a, b), maximum(a.b, b))
end

# Define maximum for a scalar with MyStructs1149
function Base.maximum(a::AbstractFloat, b::MyStructs1149)
    return MyStructs1149(max.(a, b.a), maximum(a, b.b))
end

# Define maximum for MyStructs1149
function Base.maximum(a::MyStructs1149)
    return maximum(a.a)
end

# Define maximum for a matrix of MyStructs1149
function Base.maximum(a::AbstractMatrix{MyStructs1149{T}}) where {T <: AbstractFloat}
    # Extract all `b` values from each MyStructs1149 element in the matrix and find the maximum
    return maximum(map(x -> x.b, a))
end
