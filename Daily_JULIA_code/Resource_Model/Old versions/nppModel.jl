##### IMPORTANT! if sum_over_S(H_i^0*mi) = NPP, then total biomass can be above NPP #####
begin
    # Define the Species struct with an outer constructor
    mutable struct Species
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
    end

    # Outer constructor to accept keyword arguments
    Species(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Species(m, H0, H_init, g)

    # Set NPP and competition coefficient
    NPP = 1000.0
    mu = 0.0
    
    # Function to create species_list
    function create_species_list(num_species::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                 H0_mean::Float64=10.0, H0_sd::Float64=2.0,
                                 H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
        species_list = []
        for i in 1:num_species
            m = max(0.01, rand(Normal(m_mean, m_sd)))          # Mortality rate
            H0 = 1/num_species * NPP / m
            # H0 = max(1.0, rand(Normal(H0_mean, H0_sd)))        # Characteristic density
            H_init = H0
            # H_init = max(0.1, rand(Normal(H_init_mean, H_init_sd)))  # Initial abundance
            push!(species_list, Species(m=m, H0=H0, H_init=H_init))
        end
        return species_list
    end

    # Specify the number of species
    num_species = 5

    # Create the species list
    species_list = create_species_list(num_species)

    # Function to calculate growth rates based on NPP
    function calculate_growth_rates(species_list, NPP, mu)
        S_star = length(species_list)
        # Calculate Fi for each species
        F_list = [sp.H0 * sp.m for sp in species_list]
        # Calculate the numerator of the competition term
        competition_numerator = 1 + mu * (S_star - 1)
        # Calculate gi for each species
        for (i, sp) in enumerate(species_list)
            Fi = F_list[i]
            sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
        end
    end

    # Calculate growth rates based on NPP
    calculate_growth_rates(species_list, NPP, mu)

    # Create beta matrix (competition coefficients)
    S_star = length(species_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0  # Self-competition
    end
    # Initial conditions
    H_init_values = [sp.H_init for sp in species_list]

    # Define the ODE function
    function herbivore_dynamics!(du, u, p, t)
        species_list, beta_matrix = p
        S_star = length(species_list)
        H = u
        for i in 1:S_star
            sp = species_list[i]
            # Compute competition term
            competition = 0.0
            for j in 1:S_star
                competition += beta_matrix[i, j] * H[j]
            end
            competition /= sp.H0
            # Compute the derivative
            du[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition)
        end
    end

    # Define the problem
    tspan = (0.0, 100.0)
    prob = ODEProblem(herbivore_dynamics!, H_init_values, tspan, (species_list, beta_matrix))

    # Solve the ODEs
    sol = solve(prob, Tsit5())

    # Prepare data for plotting
    t = sol.t
    H_array = Array(sol)  # Each row corresponds to a species

    # Create a plot object
    plt = PL.plot(title="Herbivore Biomasses Over Time", xlabel="Time", ylabel="Biomass")

    # Plot herbivore biomasses
    for i in 1:S_star
        PL.plot!(plt, t, H_array[i, :], label="Species $i", linewidth=2)
    end

    # Customize plot attributes
    PL.xlims!(plt, 0, tspan[2])
    PL.ylims!(plt, 0, maximum(H_array) * 1.1)

    # Display the plot
    display(plt)
end

total_biomass = sum(H_array[:, end])
NPP > total_biomass

################### LOOPING OVER AVERAGE INTERACTION STRENGTH AND NUMBER OF SPECIES ###################
#######################################################################################################
#######################################################################################################
#######################################################################################################
########### Make the loop ############
begin
    # Define the Species struct with an outer constructor
    mutable struct Species
        m::Float64        # Mortality rate
        H0::Float64       # Characteristic density (H_i^0)
        H_init::Float64   # Initial abundance (H_i(0))
        g::Float64        # Growth rate (to be calculated)
    end

    # Outer constructor to accept keyword arguments
    Species(; m::Float64, H0::Float64, H_init::Float64, g::Float64=0.0) = Species(m, H0, H_init, g)

    # Function to create species_list
    function create_species_list(num_species::Int; m_mean::Float64=0.1, m_sd::Float64=0.02,
                                 H0_mean::Float64=10.0, H0_sd::Float64=2.0,
                                 H_init_mean::Float64=5.0, H_init_sd::Float64=1.0)
        species_list = []
        for i in 1:num_species
            m = max(0.01, rand(Normal(m_mean, m_sd)))          # Mortality rate
            H0 = max(1.0, rand(Normal(H0_mean, H0_sd)))        # Characteristic density
            H_init = max(0.1, rand(Normal(H_init_mean, H_init_sd)))  # Initial abundance
            push!(species_list, Species(m=m, H0=H0, H_init=H_init))
        end
        return species_list
    end

    # Function to calculate growth rates based on NPP
    function calculate_growth_rates(species_list, NPP, mu)
        S_star = length(species_list)
        # Calculate Fi for each species
        F_list = [sp.H0 * sp.m for sp in species_list]
        # Calculate the numerator of the competition term
        competition_numerator = 1 + mu * (S_star - 1)
        # Calculate gi for each species
        for (i, sp) in enumerate(species_list)
            Fi = F_list[i]
            sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
        end
    end

    # Define the ODE function
    function herbivore_dynamics!(du, u, p, t)
        species_list, beta_matrix = p
        S_star = length(species_list)
        H = u
        for i in 1:S_star
            sp = species_list[i]
            # Compute competition term
            competition = 0.0
            for j in 1:S_star
                competition += beta_matrix[i, j] * H[j]
            end
            competition /= sp.H0
            # Compute the derivative
            du[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition)
        end
    end

    # Set up parameters for simulations
    NPP_values = collect(range(1000.0, stop=10000.0, length=20))  # NPP from 1000 to 10000
    mu_values = [0.2, 0.4, 0.6, 0.8, 1.0]  # Different mu values
    species_numbers = [3, 5, 7, 10, 20]        # Different numbers of species

    # Lists to store results
    results_mu = []       # For varying mu
    results_species = []  # For varying number of species

    # Loop over mu values and NPP (keeping number of species constant)
    num_species_fixed = 5  # Fix number of species for this part

    for mu in mu_values
        for NPP in NPP_values
            # Create species list
            species_list = create_species_list(num_species_fixed)
            # Calculate growth rates
            calculate_growth_rates(species_list, NPP, mu)
            # Create beta matrix
            S_star = length(species_list)
            beta_matrix = fill(mu, S_star, S_star)
            for i in 1:S_star
                beta_matrix[i, i] = 1.0  # Self-competition
            end
            # Initial conditions
            H_init_values = [sp.H_init for sp in species_list]
            # Define the problem
            tspan = (0.0, 100.0)
            prob = ODEProblem(herbivore_dynamics!, H_init_values, tspan, (species_list, beta_matrix))
            # Solve the ODEs
            sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
            # Get total biomass at the end
            H_array = Array(sol)
            total_biomass = sum(H_array[:, end])
            # Store the results
            push!(results_mu, (NPP=NPP, total_biomass=total_biomass, mu=mu))
        end
    end

    # Loop over number of species and NPP (keeping mu constant)
    mu_fixed = 0.5  # Fix mu for this part

    for num_species in species_numbers
        for NPP in NPP_values
            # Create species list
            species_list = create_species_list(num_species)
            # Calculate growth rates
            calculate_growth_rates(species_list, NPP, mu_fixed)
            # Create beta matrix
            S_star = length(species_list)
            beta_matrix = fill(mu_fixed, S_star, S_star)
            for i in 1:S_star
                beta_matrix[i, i] = 1.0  # Self-competition
            end
            # Initial conditions
            H_init_values = [sp.H_init for sp in species_list]
            # Define the problem
            tspan = (0.0, 100.0)
            prob = ODEProblem(herbivore_dynamics!, H_init_values, tspan, (species_list, beta_matrix))
            # Solve the ODEs
            sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
            # Get total biomass at the end
            H_array = Array(sol)
            total_biomass = sum(H_array[:, end])
            # Store the results
            push!(results_species, (NPP=NPP, total_biomass=total_biomass, num_species=num_species))
        end
    end

    # Prepare data
    df_mu = DataFrame(results_mu)
    df_species = DataFrame(results_species)
end

########### Plot total biomass vs NPP, colored by mu ###########
begin
    fig1 = Figure(resolution = (600, 500))

    ax1 = Axis(fig1[1, 1], xlabel = "NPP", ylabel = "Total Biomass",
               title = "Total Biomass vs NPP (by μ)")

    # Get unique mu values and assign colors
    unique_mu = sort(unique(df_mu.mu))
    num_colors_mu = length(unique_mu)

    # Plot data and regression line for each mu value
    for (i, mu_value) in enumerate(unique_mu)
        mask = df_mu.mu .== mu_value
        npp_data = df_mu.NPP[mask]
        biomass_data = df_mu.total_biomass[mask]

        # Scatter plot
        MK.scatter!(
            ax1, npp_data, biomass_data,
            label = "μ = $(mu_value)", markersize = 10
        )

        # Regression line using GLM
        regression_data = DataFrame(NPP = npp_data, Biomass = biomass_data)
        model = lm(@formula(Biomass ~ NPP), regression_data)
        line_x = range(minimum(npp_data), maximum(npp_data), length=100)
        line_y = coef(model)[1] .+ coef(model)[2] .* line_x
        lines!(ax1, line_x, line_y, linewidth = 2)
    end

    axislegend(ax1; position = :rb)
    fig1
end

########### Plot total biomass vs NPP, colored by number of species ###########
begin
    fig2 = Figure(resolution = (600, 500))

    ax2 = Axis(fig2[1, 1], xlabel = "NPP", ylabel = "Total Biomass",
               title = "Total Biomass vs NPP (by Species)")

    unique_species = sort(unique(df_species.num_species))
    num_colors_species = length(unique_species)

    # Plot data and regression line for each number of species
    for (i, species_num) in enumerate(unique_species)
        mask = df_species.num_species .== species_num
        npp_data = df_species.NPP[mask]
        biomass_data = df_species.total_biomass[mask]

        # Scatter plot
        MK.scatter!(
            ax2, npp_data, biomass_data,
            label = "$(species_num) Species", markersize = 10
        )

        # Regression line using GLM
        regression_data = DataFrame(NPP = npp_data, Biomass = biomass_data)
        model = lm(@formula(Biomass ~ NPP), regression_data)
        line_x = range(minimum(npp_data), maximum(npp_data), length=100)
        line_y = coef(model)[1] .+ coef(model)[2] .* line_x
        lines!(ax2, line_x, line_y, linewidth = 2)
    end

    axislegend(ax2; position = :rb)
    fig2
end

########### Plot total biomass vs NPP for both variables ###########
begin
    # Create a combined figure
    fig = Figure(resolution = (600, 500))  # Adjust resolution for side-by-side plots

    # Axis for the first plot
    ax1 = Axis(fig[1, 1], xlabel = "NPP", ylabel = "Total Biomass",
               title = "Total Biomass vs NPP (by μ)")

    # Axis for the second plot
    ax2 = Axis(fig[1, 2], xlabel = "NPP", ylabel = "Total Biomass",
               title = "Total Biomass vs NPP (by Species)")

    # Plot data for μ (first subplot)
    unique_mu = sort(unique(df_mu.mu))
    num_colors_mu = length(unique_mu)
    
    for (i, mu_value) in enumerate(unique_mu)
        mask = df_mu.mu .== mu_value
        npp_data = df_mu.NPP[mask]
        biomass_data = df_mu.total_biomass[mask]

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
        lines!(ax1, line_x, line_y, linewidth = 2)
    end

    axislegend(ax1; position = :rb)

    # Plot data for number of species (second subplot)
    unique_species = sort(unique(df_species.num_species))
    num_colors_species = length(unique_species)
    
    for (i, species_num) in enumerate(unique_species)
        mask = df_species.num_species .== species_num
        npp_data = df_species.NPP[mask]
        biomass_data = df_species.total_biomass[mask]

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
        lines!(ax2, line_x, line_y, linewidth = 2)
    end

    axislegend(ax2; position = :rb)

    # Display the combined figure
    fig
end

##########################################################################################################
