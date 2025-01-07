using LinearAlgebra
using DifferentialEquations
using CairoMakie  # For plotting
using Distributions

begin
    
    # Parameters to be decided
    legend = false
    S = 9  # Number of species
    M = 0.1  # Common mortality rate (M)
    H = 100.0  # Common characteristic density (H)
    NPP = 1000.0  # Net Primary Production
    mu = 0.0  # Average interaction strength
    ext_thr = 1.0
    final_time = 1000

    # Asymmetry parameters (between 0 and 1)
    asymmetry_qi_hi = 0.0  # Asymmetry for q_i and h_i (not used since m_i and H_i^0 are constants)
    asymmetry_competition = 1.0  # Asymmetry for competition coefficients
    exponent_abundance = 0.0  # (0 = even, 1 = strong power law) Exponent for species abundance distribution (power law)

    # Since m_i and H_i^0 are constants, we set:
    m_i = fill(M, S)
    H_i0 = fill(H, S)

    # Step 1: Generate p_i using species abundance distribution hat_H_i
    function generate_hat_Hi(S, exponent)
        ranks = collect(1:S)
        abundances = ranks .^ (-exponent)
        abundances /= sum(abundances)
        return abundances
    end

    hat_H = generate_hat_Hi(S, exponent_abundance)

    # Compute hat_p_i
    function compute_hat_p_i(hat_H, mu_matrix)
        hat_p = hat_H + mu_matrix * hat_H
        return hat_p
    end

    # Step 2: Generate competition matrix V with asymmetry
    function generate_competition_matrix(S, mu, asymmetry)
        # Symmetric competition matrix
        V_symmetric = ones(S, S)
        for i in 1:S
            for j in 1:S
                if i == j
                    V_symmetric[i, j] = 1.0
                else
                    V_symmetric[i, j] = (1 + mu)^(-1)
                end
            end
        end
        # Random competition matrix
        V_random = ones(S, S)
        for i in 1:S
            for j in 1:S
                if i != j
                    V_random[i, j] = rand()
                else
                    V_random[i, j] = 1.0
                end
            end
        end
        # Normalize off-diagonal elements to have desired mean
        off_diag_values = [ V_random[i, j] for i in 1:S, j in 1:S if i != j ]
        current_mean = mean(off_diag_values)
        desired_mean = (1 + mu)^(-1)
        scaling_factor = desired_mean / current_mean
        # Scale off-diagonal elements
        for i in 1:S
            for j in 1:S
                if i != j
                    V_random[i, j] *= scaling_factor
                end
            end
        end
        # Interpolate between symmetric and random matrices
        V = asymmetry * V_symmetric + (1.0 - asymmetry) * V_random
        return V
    end

    V = generate_competition_matrix(S, mu, asymmetry_competition)

    # Since V = (I + μ)^{-1}, we can compute μ_matrix accordingly
    mu_matrix = zeros(S, S)
    for i in 1:S
        for j in 1:S
            if i != j
                mu_matrix[i, j] = mu
            else
                mu_matrix[i, j] = 0.0
            end
        end
    end

    # Compute hat_p_i
    hat_p = compute_hat_p_i(hat_H, mu_matrix)

    # Normalize p_i
    p_i = hat_p / sum(hat_p)

    # Compute ⟨p | V p⟩
    Vp = V * p_i
    inner_product = dot(p_i, Vp)

    # Compute x using the analytical expression
    x = sqrt(NPP / (M * H * inner_product))

    println("Computed x: ", x)

    # Compute r_i and g_i
    r_i = x * p_i
    g_i = x * M .* p_i

    # Compute H_i
    H_i = x * H * (V * p_i)

    # Set up Species struct
    mutable struct Species
        id::Int
        name::String
        m::Float64        # Mortality rate (m_i)
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (g_i)
        p::Float64        # Proportion of resources allocated to species i (p_i)
    end

    # Constructor for Species
    Species(; id, name, m, H0, H_init, g, p) = Species(id, name, m, H0, H_init, g, p)

    # Create species list
    species_list = Species[]
    for i in 1:S
        sp = Species(
            id = i,
            name = "Species $i",
            m = M,
            H0 = H,
            H_init = 100.0, # H_i[i],
            g = g_i[i],
            p = p_i[i]
        )
        push!(species_list, sp)
    end

    # Step 3: Define the ecosystem dynamics function
    function ecosystem_dynamics!(du, u, p, t)
        species_list, V = p
        S = length(species_list)
        du .= 0.0  # Initialize derivatives

        for i in 1:S
            sp = species_list[i]
            H_i = u[i]
            m_i = sp.m
            g_i = sp.g
            H_i0 = sp.H0

            # Compute interaction term
            interaction = 0.0
            for j in 1:S
                H_j = u[j]
                V_ij = V[i, j]
                interaction += (1 - V_ij) * H_j  # Interaction effect
            end

            # Compute the derivative
            du[i] = H_i * m_i * ( (g_i / m_i - 1) - (H_i + interaction) / H_i0 )
        end
    end

    # Define extinction callbacks and positive domain
    callbacks = []
    push!(callbacks, PositiveDomain())
    # For herbivores
    for i in 1:length(species_list)
        condition(u, t, integrator) = u[i] - ext_thr
        affect!(integrator) = integrator.u[i] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    # Combine all callbacks
    cb = CallbackSet(callbacks...)

    # Initial abundances
    u0 = [sp.H_init for sp in species_list]

    # Time span for the simulation
    tspan = (0.0, final_time)

    # Parameters for the ODEProblem
    p = (species_list, V)

    # Define the ODE problem
    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, p)

    # Solve the ODE with callbacks
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    # Extract times and solutions
    times = sol.t
    u_array = hcat(sol.u...)  # Each column is the state at a time point
    surviving_species = length(findall(u_array[:, end] .> 0.0))
    println(surviving_species, " / ", S, " species survived.")
    
    # Plot the abundances over time
    fig = Figure(; size = (500, 400))
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Abundance", title = "Species Abundances")

    # Plot each species
    for i in 1:S
        lines!(ax, times, u_array[i, :], label = species_list[i].name)
    end
    lines!(ax, times, fill(0.0, length(times)), color = :red, linewidth = 0.5)
    if legend
        axislegend(ax)
    end

    display(fig)
end

