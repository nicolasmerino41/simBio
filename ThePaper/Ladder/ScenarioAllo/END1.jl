using EcologicalNetworksDynamics, DifferentialEquations, ForwardDiff

function plot_simulation(sol)
    fig = Figure(; size = (600, 600))
    ax = Axis(fig[1, 1])

    for i in 1:size(sol, 1)
        lines!(ax, sol.t, sol[i, :])
    end

    display(fig)
end

# 1) Build the network & model
fw = Foodweb(:niche; S=50, C=0.1)
m = default_model(fw)
B0 = equilibrium_guess(m)   # however you compute your fixed point

# 2) Solve to get B* (or use your calibration code)
sol = simulate(m, rand(50), 10000)
Bstar = sol[:, end]

Plots.plot(Bstar; legend = false)

plot_simulation(sol)

# 3) Extract Jacobian via autodiff
using ForwardDiff

function compute_jacobian(m::Model, Bstar::AbstractVector)
    # Wrap the in-place derivative function into a pure one
    g(u) = begin
        du = similar(u)
        m.dynamics!(du, u, m.params, 0.0)   # note: dynamics! not f!
        return du
    end

    # Pre-allocate the Jacobian matrix
    n = length(Bstar)
    J = Matrix{Float64}(undef, n, n)

    # Compute in-place via ForwardDiff
    ForwardDiff.jacobian!(J, g, Bstar)

    return J
end

Jfull = compute_jacobian(m, Bstar)

# 4) run your short_transform_for_ladder_step  -> (A_s,ε_s) -> build a new `Model` and repeat...
