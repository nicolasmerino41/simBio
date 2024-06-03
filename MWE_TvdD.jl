# Make sure you have the DynamicGrids#dev version
using DynamicGrids, Makie, StaticArrays
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
# Do you need this, though?? Can't we use one of the existing plot types
@Makie.recipe(MyPlot, abundance) do scene
    Attributes(
        abundance_color = :inferno
    )
end
##### We need to tell Makie how to handle arrays of mystruct for different plot types
# For a heat map we just plot the scalars
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end

# For a violin we plot each as a separate category
function Makie.convert_arguments(t::Type{<:Makie.Violin}, A::AbstractArray{<:MyStructs256, 2})
    points = mapreduce(mystruct -> Vector(mystruct.a), vcat, vec(A))
    return Makie.convert_arguments(t, repeat(1:prod(size(A)); inner = 256), points)
end

# for a scatter we just plot the first and ignore the others
function Makie.convert_arguments(t::Type{<:Makie.Scatter}, A::AbstractArray{<:MyStructs256, 2})
    points = Vector(A[1].a)
    return Makie.convert_arguments(t,points)
end

# This works fine, so for a single visualization the recipe works fine
example_matrix = [MyStructs256(SVector{256}(rand(Float32, 256) .* 90.0 .+ 10.0)) for _ in 1:3, _ in 1:3]
fig, ax, pl = myplot(example_matrix, axis = (aspect = DataAspect(),)) # for me this errors!
Makie.violin(example_matrix)
Makie.heatmap(example_matrix)
Makie.scatter(example_matrix)

##################################MakieOutput###################################
# This is just a random rule for the example
rule = Cell{:a, :a}() do data, state, I
    println(state.b)
    return MyStructs256(state.a .+ randn()) # just changed this a bit so state.a doesn't increase exponentionally
end

# The problem rises when we try to use it in a MakieOutput. It does run the simulation cause you can see it's printing
# updated values but the plot does not update
init = (a = example_matrix,
        b = example_matrix,
        c = example_matrix,
        d = example_matrix
)

# have some patience when you run this! For me it freezes for a few seconds and then catches up!
makie_output = MakieOutput(init, tspan = 1:100; fps = 3, ruleset = Ruleset(rule)) do (; layout, frame)
    ax1 = Axis(layout[1, 1])
    ax2 = Axis(layout[1, 2])
    ax3 = Axis(layout[2, 1])
    ax4 = Axis(layout[2, 2])
    Makie.heatmap!(ax1, frame.a)
    Makie.scatter!(ax2, frame.a)
   # violin!(ax3, frame.a)
   Makie.heatmap!(ax4, frame.b) # this one does not move, as we expect :)
end
