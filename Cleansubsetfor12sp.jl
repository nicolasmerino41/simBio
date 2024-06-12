# Define the species symbols for convenience
species_symbols = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k]

# Create initial states for each species
pepe_for_gif = NamedTuple{species_symbols}(fill(transposed_init, 11))

# Create reversed inits for Makie visualisation
reversed_inits = [Observable(transposed_init) for _ in 1:10]
reversed_inits = push!(reversed_inits, Observable(reversed_npp))

pepe_for_makie = NamedTuple{species_symbols}(reversed_inits)

# Define cell rules with the feeding mechanism
cell_rules = []
for i in 1:11
    prey = species_symbols[i]
    predator = i < 11 ? species_symbols[i + 1] : species_symbols[1]
    resources = species_symbols[end]  # Assuming 'k' (the last one) is the resources
    cell_rule = Cell{Tuple{prey, predator, resources}, prey}() do data, (prey_val, predator_val, resources_val), I
        return max(0.0, prey_val + growth(prey_val, self_regulation, resources_val) - 0.005 * prey_val * predator_val)
    end
    push!(cell_rules, cell_rule)
end

# Define inwards dispersal rules
indisp_rules = [InwardsDispersal{sym, sym}(formulation=ExponentialKernel(λ=0.0125), distancemethod=AreaToArea(30)) for sym in species_symbols]

# Combine all rules into a ruleset
ruleset_m = Ruleset(indisp_rules..., proc=ThreadedCPU())
pepe_for_makie
# Setup for multi-grid simulation using NamedTuple of Observables
output = MakieOutput(pepe_for_makie; tspan=1:1000, ruleset=ruleset_m,
        mask=masklayer_for_makie, fps=30) do (; layout, frame)

    # Define the global color limits
    color_limits = (10.0, 2000.0)

    # Setup the keys and titles for each plot
    plot_keys = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k, :l]  # Assuming you have these keys in your NamedTuple `pepe_for_makie`
    titles = ["Prey A", "Prey B", "Prey C", "Prey D", "Prey E", "Prey F", "Prey G", "Prey H", "Prey I", "Prey J", "Prey K", "K"]  # Custom titles for each plot

# Create a 3x4 grid layout for the 12 plots and customize each axis
    axes = [Axis(layout[i, j]; title=titles[(i-1)*4 + j]) for i in 1:3, j in 1:4]
    
    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        Makie.image!(ax, frame[key]; interpolate=false, colormap=:inferno, clims=color_limits)
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
    end
    
    # Optional: Add a colorbar to one of the axes to show the color scale
    # colorbar!(layout[1:2, end], axes[1], vertical=false)
end
