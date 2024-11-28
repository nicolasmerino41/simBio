begin
    using DifferentialEquations
    using Distributions
    using CairoMakie

    approximate = false
    plot = true
    # Set parameters
    legend = false
    num_herbivores = 5
    NPP = 100.0
    μ = 0.1  # Competition coefficient (mu)
    H0_mean_aprox = NPP / num_herbivores + 10.0  # Average characteristic density
    H0_sd = 0.0000001  # Standard deviation of characteristic density
    last_year = 1000.0
    # Herbivores:
    m_mean_h = 0.1  # Mean mortality rate
    m_sd_h = 0.00000000000002  # Standard deviation of mortality rate

    const EXTINCTION_THRESHOLD = 1e-6

    # Define the Herbivores struct
    mutable struct Herbivore
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
    end

    # Outer constructor to accept keyword arguments
    Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(m, H0, H_init, g)

    # Function to create herbivores_list
    function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=m_sd_h,
                                    H0_mean::Float64=H0_mean_aprox, H0_sd::Float64=H0_sd)
        herbivores_list = Herbivore[]
        for i in 1:num_herbivores
            m = rand(Normal(m_mean, m_sd))          # Mortality rate
            H0 = abs(rand(Normal(H0_mean, H0_sd)))  # Characteristic density
            H_init = H0  # Initial abundance set to characteristic density
            push!(herbivores_list, Herbivore(m=m, H0=H0, H_init=H_init))
        end
        return herbivores_list
    end

    # Function to calculate growth rates based on NPP
    function calculate_growth_rates(herbivores_list, NPP, μ)
        S_star = length(herbivores_list)
        # Calculate Fi for each herbivore
        F_list = [sp.H0 * sp.m for sp in herbivores_list]
        # Calculate the numerator of the competition term
        competition_numerator = 1 + μ * (S_star - 1)
        # Calculate gi for each herbivore
        for (i, sp) in enumerate(herbivores_list)
            Fi = F_list[i]
            if approximate
                sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))  # The approximation
            else
                sp.g = sp.m * ((1 + sqrt(1 + ((4 * competition_numerator * NPP) / (S_star * Fi)))) / 2)  # The exact
            end
        end
    end

    # Function to generate the competition matrix β
    function create_beta_matrix(S_star::Int, μ::Float64)
        β = fill(μ, S_star, S_star)
        for i in 1:S_star
            β[i, i] = 1.0  # Self-competition
        end
        return β
    end

    # Ecosystem dynamics function without predators
    function ecosystem_dynamics!(du, u, p, t)
        herbivores_list, β = p
        S_star = length(herbivores_list)
        H = u[1:S_star]  # Herbivore densities
        du_H = zeros(S_star)

        # Herbivore dynamics
        for i in 1:S_star
            if H[i] > EXTINCTION_THRESHOLD  # Only update if species is alive
                sp = herbivores_list[i]
                # Compute competition term
                competition = sum(β[i, :] .* H) / sp.H0
                # Compute derivative for herbivores using Eq.1
                du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition)
            else
                du_H[i] = 0.0  # Keep derivative at zero if extinct
            end
        end

        # Assign derivatives
        du[1:S_star] = du_H
    end

    # Create herbivores_list and calculate growth rates
    herbivores_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean_aprox)
    calculate_growth_rates(herbivores_list, NPP, μ)
    growth_rates = [sp.g for sp in herbivores_list]

    # Create beta_matrix
    S_star = length(herbivores_list)
    β = create_beta_matrix(S_star, μ)

    # Initial conditions for herbivores
    H_init_values = Float64[sp.H_init for sp in herbivores_list]

    # Combined initial conditions for the system
    u_init = H_init_values

    # Define the time span for simulation
    tspan = (0.0, last_year)

    # Define the ODE problem
    p = (herbivores_list, β)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Define extinction callbacks and positive domain
    callbacks = []
    push!(callbacks, PositiveDomain())
    # For herbivores
    for i in 1:length(herbivores_list)
        condition(u, t, integrator) = u[i] - EXTINCTION_THRESHOLD
        affect!(integrator) = integrator.u[i] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    # Combine all callbacks
    cb = CallbackSet(callbacks...)

    # Solve the ODE with callbacks
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    # Extract time series data
    times = sol.t
    herbivore_data = sol[1:length(herbivores_list), :]  # Herbivore dynamics

    # Calculate biomasses, excluding extinct species
    final_H = herbivore_data[:, end]
    herbivore_biomass = sum(final_H[final_H .> EXTINCTION_THRESHOLD])

    # Calculate NPP and sum of g_i * H_i
    NPP_calculated = sum([sp.g * H_i for (sp, H_i) in zip(herbivores_list, final_H)])
    holding = isapprox(NPP, NPP_calculated; atol=1e-6)

    # Print results
    num_survived_herbivores = count(final_H .> EXTINCTION_THRESHOLD)
    println("$num_survived_herbivores/$num_herbivores herbivore(s) survived.")
    println("Herbivore biomass: ", round(herbivore_biomass, digits=2))
    println("Is NPP = ∑g_i H_i? ", holding)
    println("NPP = $NPP & ∑g_i H_i = ", NPP_calculated)
    println("Difference = ", NPP - NPP_calculated)

    if plot
        # Create a figure
        fig = Figure(; size = (600, 500))
        ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Density",
                  title = "Herbivore Dynamics Over Time")

        # Plot herbivore dynamics
        for i in 1:length(herbivores_list)
            lines!(ax, times, herbivore_data[i, :], label = "Herbivore $(i)")
        end

        # Add a legend if desired
        if legend
            axislegend(ax; position = :rt)
        end

        Makie.ylims!(ax, 0, 1.1 * maximum(herbivore_data))
        # Display the figure
        display(fig)
    end
end
