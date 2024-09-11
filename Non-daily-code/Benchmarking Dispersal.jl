using Pkg
PC = "MM-1"
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\JuliaSimulation\\simBio"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
using Profile
using ProfileView
using BenchmarkTools
using Dispersal
using DynamicGrids

# Define a mask
mask_data = [i == 1 || i == 100 || j == 1 || j == 100 ? false : true for i in 1:100, j in 1:100]

# Create a grid with empty borders matching the mask
init = map(x -> x ? rand(0.0:100.0) : 0.0, mask_data)

# Create OutwardsDispersal without a mask, NoMask is default
outdisp_without_mask = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)

# Run the simulation with a mask
output_with_mask = ArrayOutput(init; tspan=1:1000, mask=mask_data)
@btime sim!(output_with_mask, Ruleset(outdisp_without_mask; boundary=Reflect()))
r = sim!(output_with_mask, Ruleset(outdisp_without_mask; boundary=Reflect()))
sum(r[1]) ≈ sum(r[1000])

# Run the simulation without a mask
output_without_mask = ArrayOutput(init; tspan=1:1000)
@btime sim!(output_without_mask, Ruleset(outdisp_without_mask; boundary=Reflect()))
# @profile sim!(output_without_mask, Ruleset(outdisp_without_mask; boundary=Reflect()))
# ProfileView.view()
r = sim!(output_without_mask, Ruleset(outdisp_without_mask; boundary=Reflect()))
sum(r[1000]) ≈ sum(r[1])

## TESTING NEW IMPLEMENTATION
# Create OutwardsDispersal with a mask
outdisp_with_mask = OutwardsDispersal(
    formulation=Dispersal.ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30),
    mask_flag=Dispersal.Mask()
)

# Run the simulation with a mask
output_with_mask_and_masked_rule = ArrayOutput(init; tspan=1:1000, mask=mask_data)
# @btime sim!(output_with_mask_and_masked_rule, Ruleset(outdisp_with_mask; boundary=Reflect()))
# @profile sim!(output_with_mask_and_masked_rule, Ruleset(outdisp_with_mask; boundary=Reflect()))
# ProfileView.view()
r = sim!(output_with_mask_and_masked_rule, Ruleset(outdisp_with_mask; boundary=Reflect()))
sum(r[1]) ≈ sum(r[1000])

indisp = InwardsDispersal{}(;
    formulation=Dispersal.ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)

output_with_indisp = ArrayOutput(init; tspan=1:1000)

r = sim!(output_with_indisp, Ruleset(indisp; boundary=Wrap()))
sum(r[1]) ≈ sum(r[1000])
###################### TESTING FROM THE PACKAGE ####################
# Define a mask
mask_data = [i == 1 || i == 10 || j == 1 || j == 10 ? false : true for i in 1:10, j in 1:10]

# Create OutwardsDispersal with a mask
outdisp_with_mask = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30),
    mask_flag=Dispersal.Mask()
)

# Create OutwardsDispersal without a mask, NoMask is default
outdisp_without_mask = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)

# Create a grid with empty borders matching the mask
init = map(x -> x ? 100.0 : 0.0, mask_data)

# Create ruleset and outputs
rule_with_mask = Ruleset(outdisp_with_mask; boundary=Reflect())
rule_without_mask = Ruleset(outdisp_without_mask; boundary=Reflect())

# Run the simulation with a mask
output_with_mask = ArrayOutput(init; tspan=1:1000, mask=mask_data)
a = sim!(output_with_mask, rule_with_mask)

sum(a[1]) ≈ sum(a[1000]) # Floating error should be smaller than 1.0

# Run the simulation without a mask to check default works fine
output_without_mask = ArrayOutput(init; tspan=1:1000)
b = sim!(output_without_mask, rule_without_mask)

sum(b[1]) ≈ sum(b[1000]) # Floating error should be smaller than 1.0

# Run the simulation with a mask but outdisp_without_mask
output_without_mask = ArrayOutput(init; tspan=1:1000, mask=mask_data)
b = sim!(output_without_mask, rule_without_mask)

sum(b[1]) > sum(b[1000]) 
#= Floating error should be larger than 1.0
because this does not identify the mask properly =#

# NEW CODE

struct Mask end
struct NoMask end

struct OutwardsDispersal{R,W,S<:Stencils.AbstractKernelStencil, M} <: SetNeighborhoodRule{R,W}
    stencil::S
    mask_flag::M
end

# Constructors for OutwardsDispersal
function OutwardsDispersal{R,W}(stencil::S; mask_flag::Union{Mask, NoMask}=NoMask()) where {R,W,S<:Stencils.AbstractKernelStencil}
    OutwardsDispersal{R,W,S,typeof(mask_flag)}(stencil, mask_flag)
end

function OutwardsDispersal{R,W}(; mask_flag::Union{Mask, NoMask}=NoMask(), kw...) where {R,W}
    stencil = DispersalKernel(; kw...)
    OutwardsDispersal{R,W,typeof(stencil),typeof(mask_flag)}(stencil, mask_flag)
end

# @inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
#     N == zero(N) && return nothing
#     mask_data = rule.mask_flag === NoMask() ? nothing : DynamicGrids.mask(data)
#     sum = zero(N)

#     if isnothing(mask_data)
#         # If there is no mask
#         for (offset, k) in zip(offsets(rule), kernel(rule))
#             @inbounds propagules = N * k
#             @inbounds add!(data[W], propagules, I .+ offset...)
#             sum += propagules
#         end
#     elseif !mask_data[I...]
#         # If there is a mask and the source cell is masked
#         return nothing
#     else
#         for (offset, k) in zip(offsets(rule), kernel(rule))
#             (target_mod, inbounds) = inbounds(data, I .+ offset)
#             if inbounds && mask_data[target_mod...]
#                 @inbounds propagules = N * k  
#                 @inbounds add!(data[W], propagules, target_mod...)  
#                 sum += propagules
#             end
#         end
#     end

#     @inbounds sub!(data[W], sum, I...)
#     return nothing
# end

@inline function applyrule!(data, rule::OutwardsDispersal{R,W}, N, I) where {R,W}
    N == zero(N) && return nothing

    # Check if the current cell is masked, skip if it is
    mask_data = if rule.mask_flag === NoMask() nothing else DynamicGrids.mask(data) end
    if !isnothing(mask_data) && !mask_data[I...]
        return nothing
    end

    sum = zero(N)
    for (offset, k) in zip(offsets(rule), kernel(rule))
        target = I .+ offset
        (target_mod, inbounds) = DynamicGrids.inbounds(data, target)
        if inbounds && (isnothing(mask_data) || mask_data[target_mod...])
            @inbounds propagules = N * k  
            @inbounds add!(data[W], propagules, target_mod...)  
            sum += propagules
        end
    end
    @inbounds sub!(data[W], sum, I...)
    return nothing
end