using LinearAlgebra
using DifferentialEquations
using CairoMakie 
using Roots
using Distributions

begin
    
    # Parameters to be decided
    S = 5  # Number of species
    bar_M = 0.1  # Mean mortality rate
    bar_H = 100.0  # Mean characteristic density
    NPP = 1000.0  # Net Primary Production
    mu = 0.5  # Average interaction strength

    # Asymmetry parameters (between 0 and 1)
    asymmetry_qi_hi = 1.0  # Asymmetry for q_i and h_i (0: random, 1: uniform)
    asymmetry_competition = 0.5  # Asymmetry for competition coefficients
    exponent_abundance = 1.0  # Exponent for species abundance distribution (power law)

    # Step 1: Generate q_i and h_i with asymmetry
    function generate_qi_hi(S, asymmetry)
        if asymmetry == 1.0
            qi = fill(1.0 / S, S)
            hi = fill(1.0 / S, S)
        else
            alpha = asymmetry * 10.0  # Scale alpha to control concentration
            if asymmetry == 0.0
                alpha = 0.1  # Small value to allow variability
            end
            dirichlet_dist = Dirichlet(alpha * ones(S))
            qi = rand(dirichlet_dist)
            hi = rand(dirichlet_dist)
        end
        return qi, hi
    end

    q_i, h_i = generate_qi_hi(S, asymmetry_qi_hi)
    
    # Step 2: Generate species abundance distribution hat_H_i
    function generate_hat_Hi(S, exponent)
        ranks = collect(1:S)
        abundances = ranks .^ (-exponent)
        abundances /= sum(abundances)
        return abundances
    end

    hat_H = generate_hat_Hi(S, exponent_abundance)

    # Step 3: Compute m_i and H_i0
    m_i = S * bar_M * q_i
    H_i0 = S * bar_H * h_i

    # Step 4: Generate competition matrix V with asymmetry
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

    # Step 5: Compute hat_p_i
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

    hat_p = hat_H + mu_matrix * hat_H

    # Normalize p_i
    p_i = hat_p / sum(hat_p)

    # # Step 6: Define a function to compute NPP given x
    # function compute_NPP(x)
    #     # Compute r_i and g_i
    #     r_i = x * p_i
    #     g_i = x * m_i .* p_i

    #     # Compute H_i
    #     H_i = x * sum(H_i0) * (V * (h_i .* p_i))

    #     # Compute NPP
    #     NPP_computed = sum(g_i .* H_i)

    #     return NPP_computed
    # end

    # # Function to compute residual between computed NPP and given NPP
    # function NPP_residual(x)
    #     return compute_NPP(x) - NPP
    # end

    # NEW
    x = r_i = x * p_i
    g_i = x * m_i .* p_i
    f = 
    x = sqrt((NPP/(bar_M*bar_H))* 1/f)
    
    
    # Find x that makes NPP_residual(x) = 0
    x_initial = 1.0
    x_solution = find_zero(NPP_residual, x_initial, method=Roots.Brent())

    println("Found x: ", x_solution)

    # Step 7: Compute final parameters
    x = x_solution
    r_i = x * p_i
    g_i = x * m_i .* p_i

    # Compute H_i
    H_i = x * sum(H_i0) * (V * (h_i .* p_i))

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
    # Outer constructor to accept keyword arguments
    Species(; id::Int, name::String, m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0, p::Float64=0.0) = 
        Species(id, name, m, H0, H_init, g, p)
    # Create species list
    species_list = Species[]
    for i in 1:S
        sp = Species(
            id = i,
            name = "Species $i",
            m = m_i[i],
            H0 = H_i0[i],
            H_init = 100.0,
            g = g_i[i],
            p = p_i[i]
        )
        push!(species_list, sp)
    end

    # Step 8: Define the ecosystem dynamics function
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

    # Initial abundances
    u0 = [sp.H_init for sp in species_list]

    # Time span for the simulation
    tspan = (0.0, 100.0)

    # Parameters for the ODEProblem
    p = (species_list, V)

    # Define the ODE problem
    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, p)

    # Solve the ODE
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract times and solutions
    times = sol.t
    u_array = hcat(sol.u...)  # Each column is the state at a time point

    # Plot the abundances over time
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Abundance", title = "Species Abundances")

    # Plot each species
    for i in 1:S
        lines!(ax, times, u_array[i, :], label = species_list[i].name)
    end

    axislegend(ax)

    display(fig)
end
