#################################################################################################
######################## DEFINING BASIC mMYSTRUCTS256 METHODS ####################################
#################################################################################################
struct mMyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::Vector{T}
    b::T

    # Custom constructor for automatic sum calculation
    function mMyStructs256(a::Vector{T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end

    # Explicit constructor allowing manual setting of both `a` and `b`
    mMyStructs256(a::Vector{T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for mMyStructs256
Base.zero(::Type{mMyStructs256{T}}) where {T <: AbstractFloat} = mMyStructs256(fill(zero(T), 256), zero(T))
Base.oneunit(::Type{mMyStructs256{T}}) where {T <: AbstractFloat} = mMyStructs256(fill(oneunit(T), 256), oneunit(T))

# Comparison based on 'b' field
Base.isless(x::mMyStructs256, y::mMyStructs256) = isless(x.b, y.b)
Base.isless(x::mMyStructs256, y::AbstractFloat) = isless(x.b, y)

# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::mMyStructs256, y::mMyStructs256) = mMyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::mMyStructs256, y::mMyStructs256) = mMyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::mMyStructs256, scalar::Real) = mMyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::mMyStructs256, scalar::Real) = mMyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::mMyStructs256, scalar::Real) = mMyStructs256(x.a .- scalar, x.b - scalar * 256)
Base.:+(x::mMyStructs256, scalar::Real) = mMyStructs256(x.a .+ scalar, x.b + scalar * 256)

# Define what a NaN is for mMyStructs256
Base.isnan(x::mMyStructs256) = isnan(x.b) || any(isnan, x.a)

# Adding a method in the sum function for mMyStructs256
function Base.sum(structs::mMyStructs256...)
    # Sum the 'a' vectors
    summed_a = reduce(+, [s.a for s in structs])

    # Sum the 'b' values
    summed_b = sum([s.b for s in structs])

    # Create a new mMyStructs256 instance with the summed results
    return mMyStructs256(summed_a, summed_b)
end

# Adding a method to maximum
# Define maximum for mMyStructs256
function Base.maximum(a::mMyStructs256, b::mMyStructs256)
    return mMyStructs256(max.(a.a, b.a))
end

# Define maximum for mMyStructs256 with a scalar
function Base.maximum(a::mMyStructs256, b::AbstractFloat)
    return mMyStructs256(max.(a.a, b))
end

# Define maximum for a scalar with mMyStructs256
function Base.maximum(a::AbstractFloat, b::mMyStructs256)
    return mMyStructs256(max.(a, b.a))
end
################## mMYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    α::Float64
end

abstract type AbstractKernelNeighborhood end

struct CustomDispersalKernel{N<:DG.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DG.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

# Define neighbors for custom kernel
function DynamicGrids.neighbors(kernel::CustomDispersalKernel, hood, center::mMyStructs256, I)
    result_a = zero(center.a)
    for i in 1:length(center.a)
        for (j, neighbor) in enumerate(hood)
            if center.a[i] > 0.0
                dist = distance(I, hood.coords[j])
                result_a[i] += kernel.formulation(dist) * neighbor.a[i]
            end
        end
    end
    return mMyStructs256(result_a)
end

# Define kernel product for mMyStructs256
function Dispersal.kernelproduct(hood::Window{1, 2, 9, mMyStructs256{Float64}}, kernel::SVector{9, Float64})
    result_a = fill(0.0, 256)
    
    for (i, k) in enumerate(kernel)
        result_a .= round.(hood[i].a .* k, sigdigits = 2) .+ result_a
    end
    
    return mMyStructs256(result_a)
end
