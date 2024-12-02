using GLM
using CairoMakie
########## NEW TRY WITH IM ############################
begin
    approximate = nothing
    plot = true
    # Set parameters
    legend = false
    num_herbivores = 1
    num_predators = 0
    NPP = 1000.0
    mu = 0.8
    H0_mean_aprox = NPP/num_herbivores  # Average characteristic density
    H0_sd = 0.0000001  # Standard deviation of characteristic density
    connectivity = 0.75  # Connectivity for interaction matrix IM
    last_year = 1000
    # Herbivores:
    m_mean_h = 0.1  # Mean mortality rate
    m_sd_h = 0.01  # Standard deviation of mortality rate
    H_init = NPP/num_herbivores
    p_i_mean = 1/num_herbivores
    p_i_sd = 0.0 # Standard deviation of p_i
    # Predator:
    m_mean_p = 0.3
    a_mean_p = 0.01
    h_mean_p = 0.1
    e_mean_p = 0.1
    c_mean_p = 1.0  # Self-regulation coefficient mean
    P_init_mean = H0_mean_aprox

    is_pi_NPP_larger_than_fi_over_4 = p_i_mean*NPP > m_mean_h*H0_mean_aprox/4
    println("is_pi_NPP_larger_than_fi_over_4: ", is_pi_NPP_larger_than_fi_over_4)
    println(p_i_mean*NPP, " Vs ", m_mean_h*H0_mean_aprox/4)
    const EXTINCTION_THRESHOLD = 1e-6
    ext_thr = 10.0

    # Define the Herbivores struct
    mutable struct Herbivores
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
        p_i::Float64      # Portion
    end

    # Outer constructor to accept keyword arguments
    Herbivores(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0, p_i::Float64) = Herbivores(m, H0, H_init, g, p_i)

    # Function to create herbivores_list
    function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=m_sd_h,
                                    H0_mean::Float64=H0_mean_aprox, H0_sd::Float64=H0_sd,
                                    H_init::Float64=H0_mean_aprox,
                                    p_i_mean::Float64=1/num_herbivores, p_i_sd::Float64=p_i_sd)
        herbivores_list = Herbivores[]
        for i in 1:num_herbivores
            m = abs(rand(Normal(m_mean, m_sd)))          # Mortality rate
            H0 = abs(rand(Normal(H0_mean, H0_sd)))        # Characteristic density
            p_i = abs(rand(Normal(p_i_mean, p_i_sd)))
            push!(herbivores_list, Herbivores(m=m, H0=H0, H_init=H_init, p_i=p_i))
        end
        return herbivores_list
    end

    # Function to calculate growth rates based on NPP
    function calculate_growth_rates(herbivores_list, NPP, mu)
        S_star = length(herbivores_list)
        # Calculate Fi for each herbivore
        F_list = [sp.H0 * sp.m for sp in herbivores_list]
        # Calculate the numerator of the competition term
        competition_numerator = 1 + mu * (S_star - 1)
        # Calculate gi for each herbivore
        for (i, sp) in enumerate(herbivores_list)
            Fi = F_list[i]
            # if approximate == true
            #     sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi)) # The approximation
            # elseif approximate == false  
            #     sp.g = sp.m * ((1+sqrt(1+((4*competition_numerator*NPP)/(S_star*Fi))))/2) # The exact
            # elseif isnothing(approximate)
                sp.g = sp.m * ((1+sqrt(1+((4*competition_numerator*sp.p_i*NPP)/(Fi))))/2) # The new
            # end
        end
    end

    # Predator struct definition
    mutable struct Predator
        m::Float64        # Mortality rate
        a::Float64        # Attack rate
        h::Float64        # Handling time
        e::Float64        # Conversion efficiency
        P_init::Float64   # Initial abundance
        c::Float64        # Self-regulation coefficient
    end

    # Predator constructor
    Predator(; m::Float64, a::Float64, h::Float64, e::Float64, P_init::Float64, c::Float64) = Predator(m, a, h, e, P_init, c)

    # Function to create a list of predators with adjusted parameters
    function create_predator_list(num_predators::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                  a_mean::Float64=0.01, a_sd::Float64=0.0001,
                                  h_mean::Float64=0.1, h_sd::Float64=0.01,
                                  e_mean::Float64=0.1, e_sd::Float64=0.01,
                                  P_init_mean::Float64=P_init_mean, P_init_sd::Float64=P_init_mean/10,
                                  c_mean::Float64=0.1, c_sd::Float64=0.01)
        predator_list = Predator[]
        for _ in 1:num_predators
            m = rand(Normal(m_mean, m_sd))              # Mortality rate
            a = rand(Normal(a_mean, a_sd))            # Attack rate
            h = rand(Normal(h_mean, h_sd))              # Handling time
            e = rand(Normal(e_mean, e_sd))              # Conversion efficiency
            P_init = rand(Normal(P_init_mean, P_init_sd))  # Initial abundance
            c = rand(Normal(c_mean, c_sd))              # Self-regulation coefficient
            push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init, c=c))
        end
        return predator_list
    end

    # Function to generate the interaction matrix
    function generate_interaction_matrix(num_predators::Int, num_prey::Int, connectivity::Float64)
        IM = zeros(Bool, num_predators, num_prey)
        for k in 1:num_predators
            for i in 1:num_prey
                IM[k, i] = rand() < connectivity  # Assign 1 with probability equal to connectivity
            end
        end
        return IM
    end

    # Ecosystem dynamics function with IM and self-regulation
    function ecosystem_dynamics!(du, u, p, t)
        herbivores_list, beta_matrix, predator_list, IM = p
        S_star = length(herbivores_list)
        num_predators = length(predator_list)
        H = u[1:S_star]  # Herbivore densities
        P = u[S_star+1:end]  # Predator densities
        du_H = zeros(S_star)
        du_P = zeros(num_predators)

        # Herbivore dynamics
        for i in 1:S_star
            if H[i] > 0  # Only update if species is alive
                sp = herbivores_list[i]
                # Compute competition term
                competition = sum(beta_matrix[i, :] .* H) / sp.H0
                # Compute predation term
                predation = 0.0
                for k in 1:num_predators
                    if IM[k, i] && P[k] > 0
                        pred = predator_list[k]
                        f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                        predation += P[k] * f_ki
                    end
                end
                # Compute derivative for herbivores
                du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
            else
                du_H[i] = 0.0  # Keep derivative at zero if extinct
            end
        end

        # Predator dynamics with self-regulation
        for k in 1:num_predators
            if P[k] > 0  # Only update if species is alive
                pred = predator_list[k]
                ingestion = 0.0
                for i in 1:S_star
                    if IM[k, i] && H[i] > 0
                        f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                        ingestion += f_ki
                    end
                end
                # Compute derivative for predators with self-regulation
                du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
            else
                du_P[k] = 0.0  # Keep derivative at zero if extinct
            end
        end

        # Assign derivatives
        du[1:S_star] = du_H
        du[S_star+1:end] = du_P
    end

    # Create herbivores_list and calculate growth rates
    herbivores_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean_aprox, H_init = H_init, p_i_mean=p_i_mean)
    calculate_growth_rates(herbivores_list, NPP, mu)  # Use positional arguments
    growth_rates = [sp.g for sp in herbivores_list]
    # Create beta_matrix
    S_star = length(herbivores_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0  # Self-competition
    end

    # Create predator list
    predator_list = create_predator_list(
        num_predators;
        m_mean=m_mean_p, a_mean=a_mean_p,
        h_mean=h_mean_p, e_mean=e_mean_p,
        c_mean=c_mean_p
    )

    # Generate the interaction matrix
    IM = generate_interaction_matrix(num_predators, S_star, connectivity)

    # Initial conditions for herbivores and predators
    H_init_values = Float64[sp.H_init for sp in herbivores_list]
    P_init_values = Float64[pred.P_init for pred in predator_list]

    # Combined initial conditions for the system
    u_init = vcat(H_init_values, P_init_values)

    # Define the time span for simulation
    tspan = (0.0, last_year)

    # Define the ODE problem
    p = (herbivores_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Define extinction callbacks and positive domain
    callbacks = []
    push!(callbacks, PositiveDomain())
    # For herbivores
    for i in 1:length(herbivores_list)
        condition(u, t, integrator) = u[i] - ext_thr
        affect!(integrator) = integrator.u[i] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    # For predators
    offset = length(herbivores_list)
    for i in 1:length(predator_list)
        idx = offset + i  # Adjust index for predators
        condition(u, t, integrator) = u[idx] - ext_thr
        affect!(integrator) = integrator.u[idx] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    # Combine all callbacks
    cb = CallbackSet(callbacks...)

    # Solve the ODE with callbacks
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    # Extract time series data
    times = sol.t
    herbivore_data = sol[1:length(herbivores_list), :]  # Herbivore dynamics
    predator_data = sol[length(herbivores_list)+1:end, :]  # Predator dynamics

    # Calculate biomasses, excluding extinct species
    herbivore_biomass = sum(herbivore_data[:, end][herbivore_data[:, end] .> ext_thr])
    predator_biomass = sum(predator_data[:, end][predator_data[:, end] .> ext_thr])

    # Equation holding true?
    holding = (NPP/sum(growth_rates .* herbivore_data[:, end]) > 0.99) && (NPP/sum(growth_rates .* herbivore_data[:, end]) < 1.001)

    pred_herb_ratio = predator_biomass / herbivore_biomass
    total_biomass = herbivore_biomass + predator_biomass

    # Count surviving species
    num_survived_herbivores = count(herbivore_data[:, end] .> ext_thr)
    println("$num_survived_herbivores/$num_herbivores herbivore(s) survived.")

    num_survived_predators = count(predator_data[:, end] .> ext_thr)
    println("$num_survived_predators/$num_predators predator(s) survived.")

    println("Herbivore biomass ", round(herbivore_biomass, digits=2))
    println("Predator biomass ", round(predator_biomass, digits=2))
    println("Predator/herbivore ratio ", pred_herb_ratio)
    println("Total_biomass/NPP = ", round(total_biomass/NPP, digits=2))
    println("Is NPP ≈ ∑g_iH_i? ", holding)
    println("NPP = $NPP & ∑g_iH_i = ", round(sum(growth_rates .* herbivore_data[:, end]), digits=2))
    println("∑g_iH_i/NPP = ", round(sum(growth_rates .* herbivore_data[:, end])/NPP, digits=2))
    
    if plot
    # Create a single figure
    fig = Figure(; size = (500, 400))

    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Density",
                 title = "Herbivore and Predator Dynamics Over Time")

    # Plot herbivore dynamics with solid lines
    for i in 1:length(herbivores_list)
        lines!(ax, times, herbivore_data[i, :], label = "Herbivore $(i)")
    end

    # Plot predator dynamics with dashed lines
    for k in 1:length(predator_list)
        lines!(ax, times, predator_data[k, :], label = "Predator $(k)", linestyle = :dash)
    end

    # Add a legend
    if legend
        axislegend(ax; position = :rt)
    end
    Makie.ylims!(ax, 0, 1.1 * maximum(vcat(herbivore_data, predator_data)))
    # Display the figure
    display(fig)
    end 
end

