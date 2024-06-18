# Questions
## 1. What are the substitutues for these lines:
humanpop = GDALarray(humanpop_filepath; mappedcrs=EPSG(4326))[Band(1), aust...] |>
    A -> replace_missing(A, missing) |>
    A -> permutedims(A, (Lat, Lon)) |> 
    A -> reorder(A, Lat=ReverseArray, Lon=ForwardArray) 
end
rast[X(2), Y(3)]
## 2. What does this notation mean:
Dispersal.LogisticGrowth — Type
LogisticGrowth <: GrowthRule

LogisticGrowth(; rate, carrycap, timestep, [nsteps_type])
LogisticGrowth{R}(; rate, carrycap, timestep, [nsteps_type])
LogisticGrowth{R,W}(; rate, carrycap, timestep, [nsteps_type])

:hi
## 3. What does this error mean:
ERROR: MethodError: no method matching _asiterable(::Symbol)
The applicable method may be too new: running in world age 36144, while current world is 37702.

Closest candidates are:
  _asiterable(::Symbol) (method too new to be called from this world context.)
   @ DynamicGrids C:\Users\nicol\.julia\packages\DynamicGrids\FbmEN\src\utils.jl:82

## 4. What is going on with Dispersal vs DynamicGrids incompatibilities?
(simBio) pkg> add DynamicGrids#dev
    Updating git-repo `https://github.com/cesaraustralia/DynamicGrids.jl.git`
   Resolving package versions...
ERROR: Unsatisfiable requirements detected for package Dispersal [8797f018]:
 Dispersal [8797f018] log:
 ├─possible versions are: 0.1.0-0.6.0 or uninstalled
 ├─restricted to versions * by an explicit requirement, leaving only versions: 0.1.0-0.6.0
 └─restricted by compatibility requirements with DynamicGrids [a5dba43e] to versions: uninstalled — no versions left
   └─DynamicGrids [a5dba43e] log:
     ├─possible versions are: 0.21.4 or uninstalled
     └─DynamicGrids [a5dba43e] is fixed to version 0.21.4
end

## 5. MakieOutput
# Create our own plots with Makie.jl
output = MakieOutput(rand(Bool, 200, 300); tspan, ruleset) do layout, frame, t
  image!(Axis(layout[1, 1]), frame; interpolate=false, colormap=:inferno)
end
## 6. GridLayout
help?> GridLayout
search: GridLayout

  GridLayout(; kwargs...)

  Create a GridLayout without parent and with size [1, 1].

  ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────  

  GridLayout(g::Union{GridPosition, GridSubposition}, args...; kwargs...)

  Create a GridLayout at position g in the parent GridLayout of g if it is a GridPosition and in a nested child GridLayout if it is a GridSubposition. The args and   
  kwargs are passed on to the normal GridLayout constructor.
end

