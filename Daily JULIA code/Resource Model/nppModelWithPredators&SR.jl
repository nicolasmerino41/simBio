using GLM
using CairoMakie
########## NEW TRY WITH IM ############################
begin
    
    # Set parameters
    legend = false
    num_herbivores = 80
    num_predators = 40
    NPP = 1000.0
    mu = 0.5
    H0_mean_aprox = NPP / num_herbivores
    connectivity = 0.4  # Connectivity for interaction matrix IM
    last_year = 2000
    # Herbivores:
    m_mean_h = 0.1
    # Predator:
    m_mean_p = 0.1
    a_mean_p = 0.1
    h_mean_p = 0.1
    e_mean_p = 0.1
    c_mean_p = 1.0

    # Define the Herbivores struct
    mutable struct Herbivores
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
    end

    # Outer constructor to accept keyword arguments
    Herbivores(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivores(m, H0, H_init, g)

    # Function to create herbivores_list
    function create_herbivores_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                    H0_mean::Float64=H0_mean_aprox, H0_sd::Float64=H0_mean_aprox/10,
                                    H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
        herbivores_list = []
        for i in 1:num_herbivores
            m = rand(Normal(m_mean, m_sd))          # Mortality rate
            H0 = rand(Normal(H0_mean, H0_sd))        # Characteristic density
            H_init = H0  # Initial abundance set to characteristic density
            push!(herbivores_list, Herbivores(m=m, H0=H0, H_init=H_init))
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
            sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
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
                                  P_init_mean::Float64=5.0, P_init_sd::Float64=1.0,
                                  c_mean::Float64=0.1, c_sd::Float64=0.01)
        predator_list = []
        for _ in 1:num_predators
            m = rand(Normal(m_mean, m_sd))              # Mortality rate
            a = rand(Normal(a_mean, a_sd))            # Attack rate
            h = rand(Normal(h_mean, h_sd))              # Handling time
            e = rand(Normal(e_mean, e_sd))              # Conversion efficiency
            P_init = rand(Normal(P_init_mean, P_init_mean/10))  # Initial abundance
            c = rand(Normal(c_mean, c_sd))  # Define mean and sd for c
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

    # Ecosystem dynamics function with IM
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
            sp = herbivores_list[i]
            # Compute competition term
            competition = sum(beta_matrix[i, :] .* H) / sp.H0
            # Compute predation term
            predation = 0.0
            for k in 1:num_predators
                if IM[k, i]
                    pred = predator_list[k]
                    f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                    predation += P[k] * f_ki
                end
            end
            # Compute derivative for herbivores
            du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
        end

        # Predator dynamics with self-regulation
        for k in 1:num_predators
            pred = predator_list[k]
            ingestion = 0.0
            for i in 1:S_star
               if IM[k, i]
                  f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                  ingestion += f_ki
              end
            end
            # Self-regulation term
            du_P[k] = P[k] * (pred.e * ingestion - pred.m - pred.c * P[k])
        end
        
        # Assign derivatives
        du[1:S_star] = du_H
        du[S_star+1:end] = du_P
    end

    # Create herbivores_list and calculate growth rates
    herbivores_list = create_herbivores_list(num_herbivores; m_mean=m_mean_h, H0_mean=H0_mean_aprox)
    calculate_growth_rates(herbivores_list, NPP, mu)  # Use positional arguments

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
    H_init_values = [sp.H_init for sp in herbivores_list]
    P_init_values = [pred.P_init for pred in predator_list]

    # Combined initial conditions for the system
    u_init = vcat(H_init_values, P_init_values)

    # Define the time span for simulation
    tspan = (0.0, last_year)

    # Define the ODE problem
    p = (herbivores_list, beta_matrix, predator_list, IM)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Solve the ODE
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract time series data
    times = sol.t
    herbivore_data = sol[1:length(herbivores_list), :]  # Herbivore dynamics
    predator_data = sol[length(herbivores_list)+1:end, :]  # Predator dynamics
    herbivore_biomass = sum(sol[1:length(herbivores_list), end][sol[1:length(herbivores_list), end] .> 1.0])
    predator_biomass = sum(sol[length(herbivores_list)+1:end, end][sol[length(herbivores_list)+1:end, end] .> 1.0])
    pred_herb_ratio = predator_biomass / herbivore_biomass
    total_biomass = herbivore_biomass + predator_biomass

    if true #any(herbivore_data[:, end] .<= 1.0) 
        # num_extinct_herbivores = count(herbivore_data[:, end] .<= 1.0)
        num_survived_herbivores = count(herbivore_data[:, end] .> 1.0)
        println("$num_survived_herbivores/$num_herbivores herbivore(s) survived.")
    end
    if true #any(predator_data[:, end] .<= 1.0)
        # num_extinct_predators = count(predator_data[:, end] .<= 1.0)
        num_survived_predators = count(predator_data[:, end] .> 1.0)
        println("$num_survived_predators/$num_predators predator(s) survived.")
    end
    println("Herbivore biomass ", round(herbivore_biomass, digits=2))
    println("Predator biomass ", round(predator_biomass, digits=2))
    println("Predator/herbivore ratio ", pred_herb_ratio)
    println("Total_biomass/NPP = ", round(total_biomass/NPP, digits=2))

    # Create a single figure
    fig = MK.Figure(; size = (600, 500))

    ax = MK.Axis(fig[1, 1], xlabel = "Time", ylabel = "Density",
                 title = "Herbivore and Predator Dynamics Over Time")

    # Plot herbivore dynamics with solid lines
    for i in 1:length(herbivores_list)
        MK.lines!(ax, times, herbivore_data[i, :], label = "Herbivore $(i)")
    end

    # Plot predator dynamics with dashed lines
    for k in 1:length(predator_list)
        MK.lines!(ax, times, predator_data[k, :], label = "Predator $(k)", linestyle = :dash)
    end

    # Add a legend
    if legend
        MK.axislegend(ax; position = :rt)
    end

    # Display the figure
    display(fig)
end

################### LOOPING OVER AVERAGE INTERACTION STRENGTH AND NUMBER OF SPECIES ###################
#######################################################################################################
#######################################################################################################
#######################################################################################################
begin
    ########### Make the loop ############
    # Define the Herbivore struct
    mutable struct Herbivore
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
    end

    # Constructor for Herbivore
    Herbivore(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Herbivore(m, H0, H_init, g)

    # Function to create a list of herbivores
    function create_herbivore_list(num_herbivores::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                   H0_mean::Float64=10.0, H0_sd::Float64=2.0)
        herbivore_list = []
        for i in 1:num_herbivores
            m = max(0.01, rand(Normal(m_mean, m_sd)))          # Mortality rate
            H0 = max(1.0, rand(Normal(H0_mean, H0_sd)))        # Characteristic density
            H_init = H0                                        # Initial abundance set to H0
            push!(herbivore_list, Herbivore(m=m, H0=H0, H_init=H_init))
        end
        return herbivore_list
    end

    # Function to calculate growth rates based on NPP
    function calculate_growth_rates(herbivore_list, NPP, mu)
        S_star = length(herbivore_list)
        # Calculate Fi for each herbivore
        F_list = [sp.H0 * sp.m for sp in herbivore_list]
        # Calculate the numerator of the competition term
        competition_numerator = 1 + mu * (S_star - 1)
        # Calculate gi for each herbivore
        for (i, sp) in enumerate(herbivore_list)
            Fi = F_list[i]
            sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
        end
    end

    # Predator struct definition
    mutable struct Predator
        m::Float64        # Mortality rate
        a::Float64        # Attack rate
        h::Float64        # Handling time
        e::Float64        # Conversion efficiency
        P_init::Float64   # Initial abundance
    end

    # Predator constructor
    Predator(; m::Float64, a::Float64, h::Float64, e::Float64, P_init::Float64) = Predator(m, a, h, e, P_init)

    # Function to create a list of predators
    function create_predator_list(num_predators::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                  a_mean::Float64=0.001, a_sd::Float64=0.0001,
                                  h_mean::Float64=0.1, h_sd::Float64=0.01,
                                  e_mean::Float64=0.1, e_sd::Float64=0.01)
        predator_list = []
        for _ in 1:num_predators
            m = max(0.01, rand(Normal(m_mean, m_sd)))              # Mortality rate
            a = max(0.0001, rand(Normal(a_mean, a_sd)))            # Attack rate
            h = max(0.01, rand(Normal(h_mean, h_sd)))              # Handling time
            e = max(0.01, rand(Normal(e_mean, e_sd)))              # Conversion efficiency
            P_init = max(0.1, rand(Normal(5.0, 1.0)))              # Initial abundance
            push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init))
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

    # Ecosystem dynamics function with IM
    function ecosystem_dynamics!(du, u, p, t)
        herbivore_list, beta_matrix, predator_list, IM = p
        S_star = length(herbivore_list)
        num_predators = length(predator_list)
        H = u[1:S_star]  # Herbivore densities
        P = u[S_star+1:end]  # Predator densities
        du_H = zeros(S_star)
        du_P = zeros(num_predators)

        # Herbivore dynamics
        for i in 1:S_star
            sp = herbivore_list[i]
            # Compute competition term
            competition = sum(beta_matrix[i, :] .* H) / sp.H0
            # Compute predation term
            predation = 0.0
            for k in 1:num_predators
                if IM[k, i]
                    pred = predator_list[k]
                    f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                    predation += P[k] * f_ki
                end
            end
            # Compute derivative for herbivores
            du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
        end

        # Predator dynamics
        for k in 1:num_predators
            pred = predator_list[k]
            ingestion = 0.0
            for i in 1:S_star
                if IM[k, i]
                    f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                    ingestion += f_ki
                end
            end
            # Compute derivative for predators
            du_P[k] = P[k] * (pred.e * ingestion - pred.m)
        end

        # Assign derivatives
        du[1:S_star] = du_H
        du[S_star+1:end] = du_P
    end

    # Set up parameters for simulations
    NPP_values = collect(range(1000.0, stop=10000.0, length=20))  # NPP from 1000 to 10000
    mu_values = [0.2, 0.4, 0.6, 0.8, 1.0]  # Different mu values
    species_numbers = [3, 5, 7, 10, 20, 100]    # Different numbers of herbivore species

    # Fixed parameters
    num_predators_fixed = 5  # Number of predators fixed
    connectivity_fixed = 0.5  # Connectivity fixed

    # Lists to store results
    results_mu = []       # For varying mu
    results_species = []  # For varying number of species

    # Simulation loops (first for mu values, then for species numbers)
    num_herbivores_fixed = 10  # Fix number of herbivores for this part
    for mu in mu_values
        for NPP in NPP_values
            herbivore_list = create_herbivore_list(num_herbivores_fixed)
            calculate_growth_rates(herbivore_list, NPP, mu)
            S_star = length(herbivore_list)
            beta_matrix = fill(mu, S_star, S_star)
            for i in 1:S_star
                beta_matrix[i, i] = 1.0
            end
            predator_list = create_predator_list(num_predators_fixed)
            IM = generate_interaction_matrix(num_predators_fixed, S_star, connectivity_fixed)
            H_init_values = [sp.H_init for sp in herbivore_list]
            P_init_values = [pred.P_init for pred in predator_list]
            u_init = vcat(H_init_values, P_init_values)
            tspan = (0.0, 200.0)
            p = (herbivore_list, beta_matrix, predator_list, IM)
            prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)
            sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
            H_array = sol[1:S_star, :]
            total_biomass_end = sum(H_array[:, end])
            push!(results_mu, (NPP=NPP, total_biomass_end=total_biomass_end, mu=mu))
        end
    end

    mu_fixed = 0.5
    for num_herbivores in species_numbers
        for NPP in NPP_values
            herbivore_list = create_herbivore_list(num_herbivores)
            calculate_growth_rates(herbivore_list, NPP, mu_fixed)
            S_star = length(herbivore_list)
            beta_matrix = fill(mu_fixed, S_star, S_star)
            for i in 1:S_star
                beta_matrix[i, i] = 1.0
            end
            predator_list = create_predator_list(num_predators_fixed)
            IM = generate_interaction_matrix(num_predators_fixed, S_star, connectivity_fixed)
            H_init_values = [sp.H_init for sp in herbivore_list]
            P_init_values = [pred.P_init for pred in predator_list]
            u_init = vcat(H_init_values, P_init_values)
            tspan = (0.0, 200.0)
            p = (herbivore_list, beta_matrix, predator_list, IM)
            prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)
            sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
            H_array = sol[1:S_star, :]
            total_biomass_end = sum(H_array[:, end])
            push!(results_species, (NPP=NPP, total_biomass_end=total_biomass_end, num_species=num_herbivores))
        end
    end

    # Prepare data
    df_mu = DataFrame(results_mu)
    df_species = DataFrame(results_species)
    
    ########### Plot total biomass vs NPP for both variables ###########

    # Create a combined figure
    fig = MK.Figure(resolution = (1200, 500))  # Adjust resolution for side-by-side plots

    # Axis for the first plot
    ax1 = MK.Axis(fig[1, 1], xlabel = "NPP", ylabel = "Total Herbivore Biomass",
                  title = "Total Herbivore Biomass vs NPP (by μ)")

    # Axis for the second plot
    ax2 = MK.Axis(fig[1, 2], xlabel = "NPP", ylabel = "Total Herbivore Biomass",
                  title = "Total Herbivore Biomass vs NPP (by Number of Species)")

    # Plot data for μ (first subplot)
    unique_mu = sort(unique(df_mu.mu))
    num_colors_mu = length(unique_mu)
        
    for (i, mu_value) in enumerate(unique_mu)
        mask = df_mu.mu .== mu_value
        npp_data = df_mu.NPP[mask]
        biomass_data = df_mu.total_biomass_end[mask]

        # Scatter plot
        MK.scatter!(
            ax1, npp_data, biomass_data,
            label = "μ = $(mu_value)", markersize = 10
        )

        # Regression line
        regression_data = DataFrame(NPP = npp_data, Biomass = biomass_data)
        model = lm(@formula(Biomass ~ NPP), regression_data)
        line_x = range(minimum(npp_data), maximum(npp_data), length=100)
        line_y = coef(model)[1] .+ coef(model)[2] .* line_x
        MK.lines!(ax1, line_x, line_y, linewidth = 2)
    end

    MK.axislegend(ax1; position = :rb)

    # Plot data for number of species (second subplot)
    unique_species = sort(unique(df_species.num_species))
    num_colors_species = length(unique_species)
        
    for (i, species_num) in enumerate(unique_species)
        mask = df_species.num_species .== species_num
        npp_data = df_species.NPP[mask]
        biomass_data = df_species.total_biomass_end[mask]

        # Scatter plot
        MK.scatter!(
            ax2, npp_data, biomass_data,
            label = "$(species_num) Species", markersize = 10
        )

        # Regression line
        regression_data = DataFrame(NPP = npp_data, Biomass = biomass_data)
        model = lm(@formula(Biomass ~ NPP), regression_data)
        line_x = range(minimum(npp_data), maximum(npp_data), length=100)
        line_y = coef(model)[1] .+ coef(model)[2] .* line_x
        MK.lines!(ax2, line_x, line_y, linewidth = 2)
    end

    MK.axislegend(ax2; position = :rb)

    # Adjust layout
    fig.layout[1, 1] = ax1
    fig.layout[1, 2] = ax2

    # Display the combined figure
    display(fig)
end
