include("DA_birmmals_with_pi.jl")
include("Functions/generate_competition_matrix.jl")
include("Functions/species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")

##############################################################################################
# We apply the dynamics using the species, variables and parameters from a chosen cell.
##############################################################################################
begin

    # Choose a cell
    i, j = 22, 3

    # Extract the NPP from the cell
    NPP = Float64(npp_DA[i, j])
    
    S, R, species_names, H_i0, m_i, p_vec, x_final, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = setup_community_from_cell(
        i, j;
        NPP = NPP,
        M_mean = 0.1,
        mu = 0.5,
        symmetrical_competition = false,
        mean_m_alpha = 0.1,
        epsilon_val = 0.00008,
        mu_predation = 0.00001,
        iberian_interact_NA=iberian_interact_NA,
        species_dict=species_dict,
        m_standard_deviation = 0.0,
        h_standard_deviation = 0.0,
        artificial_pi = true
    )

    include("Callbacks.jl")

    params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
    H_init = H_i0
    P_init = H_i0[1:R] ./ 10 # fill(10.0, R)
    u0 = vcat(H_init, P_init)
    tspan = (0.0, 500.0)
    EXTINCTION_THRESHOLD = 1e-6

    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    times = sol.t
    H_data = sol[1:S, :]
    P_data = sol[S+1:S+R, :]

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time", ylabel="Density", title="Dynamics")
    for x in 1:S
        lines!(ax, times, H_data[x, :], label="H$x")
    end
    for α in 1:R
        lines!(ax, times, P_data[α, :], label="P$α", linestyle=:dash)
    end
    # axislegend(ax; position=:rt)
    display(fig)

    # Summary stats
    herb_survivors = count(H_data[:, end] .> EXTINCTION_THRESHOLD)
    pred_survivors = count(P_data[:, end] .> EXTINCTION_THRESHOLD)

    println("Herbivores survived: $herb_survivors/$S")
    println("Predators survived: $pred_survivors/$R")

    H_biomass = sum(H_data[:, end][H_data[:, end] .> EXTINCTION_THRESHOLD])
    P_biomass = sum(P_data[:, end][P_data[:, end] .> EXTINCTION_THRESHOLD])
    println("Herbivore biomass at end: ", H_biomass)
    println("Predator biomass at end: ", P_biomass)
    println("Total biomass: ", H_biomass + P_biomass)
    println("herb_pred_ratio: ", P_biomass / H_biomass)
    println("NPP == ∑g_iH_i? ", isapprox(NPP, sum(g_i .* H_data[:, end]), atol=50.0))
    println("∑g_iH_i = ", sum(g_i .* H_data[:, end]))
end
