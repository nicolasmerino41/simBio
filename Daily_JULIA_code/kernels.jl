################# MYSTRUCTS256 KERNEL METHODS ################
###############################################################
###############################################################
struct CustomKernel <: KernelFormulation
    α::AbstractFloat
end

abstract type AbstractKernelNeighborhood end

struct CustomDispersalKernel{N<:DynamicGrids.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
    neighborhood::N
    formulation::F
end

function CustomDispersalKernel(; 
    neighborhood::DynamicGrids.Neighborhood=Moore(1), 
    formulation::KernelFormulation=CustomKernel(1.0)
)
    CustomDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
end

function (kernel::CustomKernel)(distance)
    return exp(-(distance^2) / (2*(kernel.α^2)))
end

# # Define neighbors for custom kernel
# function DynamicGrids.neighbors(kernel::CustomDispersalKernel, hood, center::MyStructs256, I)
#     result_a = zero(center.a)
#     for i in 1:256
#         for (j, neighbor) in enumerate(hood)
#             if center.a[i] > 0.0
#                 dist = distance(I, hood.coords[j])
#                 result_a += kernel.formulation(dist) * neighbor.a[i]
#             end
#         end
#     end
#     return MyStructs256(result_a)
# end

# # Define kernel product for MyStructs256
# function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{AbstractFloat}}, kernel::SVector{9, AbstractFloat})
    
#     result_a = SVector{256, AbstractFloat}(fill(0.0f0, 256))
    
#     for (i, k) in enumerate(kernel)
#         result_a += hood[i].a .* k
#     end
#     return MyStructs256(result_a)
# end

# function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyStructs256{AbstractFloat}}, kernel::SVector{25, AbstractFloat})
    
#     result_a = SVector{256, AbstractFloat}(fill(0.0f0, 256))
    
#     for (i, k) in enumerate(kernel)
#         result_a += hood[i].a .* k
#     end
#     return MyStructs256(result_a)
# end

# struct BodyMassKernel <: KernelFormulation
#     α::AbstractFloat
# end

# struct BodyMassDispersalKernel{N<:DynamicGrids.Neighborhood, F<:KernelFormulation} <: AbstractKernelNeighborhood
#     neighborhood::N
#     formulation::F
# end

# function BodyMassDispersalKernel(; 
#     neighborhood::DynamicGrids.Neighborhood=Moore(1), 
#     formulation::KernelFormulation=BodyMassKernel(1.0)
# )
#     BodyMassDispersalKernel{typeof(neighborhood), typeof(formulation)}(neighborhood, formulation)
# end

# # Define kernel product for MyStructs256
# function Dispersal.kernelproduct(hood::Window{1, 2, 9, MyStructs256{AbstractFloat}}, kernel::SVector{9, AbstractFloat})
    
#     result_a = SVector{256, AbstractFloat}(fill(0.0f0, 256))
    
#     for (i, k) in enumerate(kernel)
#         result_a += hood[i].a .* body_mass_vector .* k 
#     end
#     return MyStructs256(result_a)
# end

# function Dispersal.kernelproduct(hood::Window{2, 2, 25, MyStructs256{AbstractFloat}}, kernel::SVector{25, AbstractFloat})
    
#     result_a = SVector{256, AbstractFloat}(fill(0.0f0, 256))
    
#     for (i, k) in enumerate(kernel)
#         result_a += hood[i].a .* body_mass_vector .* k 
#     end
#     return MyStructs256(result_a)
# end