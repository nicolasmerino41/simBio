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


# For a heatmap we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs1149, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyBirmmals, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyHerps, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_DA.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs1149, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > 1.0, 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyBirmmals, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_birds[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyHerps, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_herps[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
# WITH LAMBDA
# For MyStructs
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs1149, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# PLOT
# For MyStructs
# function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyStructs1149, 2})
#     richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_DA.multiplicative)
#     return Makie.convert_arguments(t, richness)
# end
# MK.image(Matrix(des_file_to_try); colomap = custom_palette)
# For MyBirmmals
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyBirmmals, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector_birds[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# For MyHerps
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractArray{<:MyHerps, 2}, lambda_grid::AbstractArray{<:AbstractFloat, 2})
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector_herps[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end
# FOR RASTER
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyStructs1149, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyBirmmals, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyHerps, 2})
    scalars = map(mystruct -> mystruct.b, A).*lambda_raster.multiplicative
    return Makie.convert_arguments(t, scalars)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyStructs1149, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map((mystruct, lambda_value) -> count(i -> (mystruct.a[i] * lambda_value) > body_mass_vector[i], 1:length(mystruct.a)), A, lambda_grid)
    return Makie.convert_arguments(t, richness)
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyBirmmals, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_birds[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyHerps, 2}, lambda_grid::AbstractRaster{<:AbstractFloat, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > body_mass_vector_herps[i], 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end