# Make sure you have the DynamicGrids#dev version
using DynamicGrids, Makie, StaticArrays, CairoMakie, GLMakie, WGLMakie
######################## DEFINING BASIC MYSTRUCTS256 METHODS ####################################
    # This is my custom struct definition, with 'a' and 'b' fields and a custom sum function
#################################################################################################
struct MyStructs256{T <: AbstractFloat} <: FieldVector{2, T}
    a::SVector{256, T}
    b::T
    
    # Custom constructor for automatic sum calculation
    function MyStructs256(a::SVector{256, T}) where T <: AbstractFloat
        new{T}(a, sum(a))
    end
    
    # Explicit constructor allowing manual setting of both `a` and `b`
    MyStructs256(a::SVector{256, T}, b::T) where T <: AbstractFloat = new{T}(a, b)
end

# Define zero and oneunit for MyStructs256
Base.zero(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(zero(T), 256)), zero(T))
Base.oneunit(::Type{MyStructs256{T}}) where {T <: AbstractFloat} = MyStructs256(SVector{256, T}(fill(oneunit(T), 256)), oneunit(T))
# Comparison based on 'b' field
Base.isless(x::MyStructs256, y::MyStructs256) = isless(x.b, y.b)
Base.isless(x::MyStructs256, y::AbstractFloat) = isless(x.b, y)
# Element-wise arithmetic operations ensuring 'b' is recalculated correctly
Base.:+(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .+ y.a, sum(x.a .+ y.a))
Base.:-(x::MyStructs256, y::MyStructs256) = MyStructs256(x.a .- y.a, sum(x.a .- y.a))
Base.:*(x::MyStructs256, scalar::Real) = MyStructs256(x.a .* scalar, sum(x.a .* scalar))
Base.:/(x::MyStructs256, scalar::Real) = MyStructs256(x.a ./ scalar, sum(x.a ./ scalar))
Base.:-(x::MyStructs256, scalar::Real) = MyStructs256(x.a .- scalar, x.b - scalar*256)
Base.:+(x::MyStructs256, scalar::Real) = MyStructs256(x.a .+ scalar, x.b + scalar*256)
# Define what a NaN is for MyStructs256
Base.isnan(x::MyStructs256) = isnan(x.b) || any(isnan, x.a)

##############################################################################
########## MyPlot is my custom plot recipe for MyStructs256 ##################
@Makie.recipe(MyPlot, abundance) do scene
    Attributes(
        abundance_color = :inferno
    )
end

function Makie.plot!(plot::MyPlot)
    # Extract the abundance matrix from the plot object
    abundance = plot[1][]  # Use `[]` to get the value from the Observable
    
    # println("Type of abundance: ", typeof(abundance))
    # Extract the `b` elements from the DimArray
    b_matrix = [mystruct.b for mystruct in abundance]
        
    # Reshape the `b_matrix` to match the original dimensions
    reshaped_b_matrix = reshape(b_matrix, size(abundance))

    # Create the heatmap plot
    Makie.heatmap!(plot, reshaped_b_matrix, colormap = plot[:abundance_color])
    
    # println("Type of plot: ", typeof(Observable(plot)))
    # Return the updated plot
    return Observable(plot)
    
end

# This works fine, so for a single visualization the recipe works fine
example_matrix = [MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0)) for _ in 1:3, _ in 1:3]
fig, ax, pl = myplot(example_matrix, axis = (aspect = DataAspect(),))

##################################MakieOutput###################################
# This is just a random rule for the example
rule = Cell{:a, :a}() do data, state, I
    println(state.b)
    return MyStructs256(state.a .+ (state.a .* rand()))
end
# The problem rises when we try to use it in a MakieOutput. It does run the simulation cause you can see it's printing
# updated values but the plot does not update
init = (a = example_matrix,
        b = example_matrix,
        c = example_matrix,
        d = example_matrix
)


makie_output = MakieOutput(init, tspan = 1:100; fps = 50, ruleset = Ruleset(rule)) do (; layout, frame)
    ax1 = Axis(layout[1, 1])
    ax2 = Axis(layout[1, 2])
    ax3 = Axis(layout[2, 1])
    ax4 = Axis(layout[2, 2])
    
    # println(typeof(frame))
    
    myplot!(ax1, frame.a; interpolate=false, colormap=:inferno)
    myplot!(ax2, frame.b; interpolate=false)
    myplot!(ax3, frame.c; interpolate=false)
    myplot!(ax4, frame.d; interpolate=false)
end