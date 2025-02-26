#################### SAME THING BUT AS A FUNCTION ######################
function single_run(
    cell, mu_val, mu_pred_val, eps_val, sym_competition;
    plot=false, sp_removed_name=nothing,
    artificial_pi=false, NPP=nothing,
    herbivore_m = 0.1,
    predator_m = 0.1
)
    
    # Placeholder for results
    AAAA = DataFrame()

    @info "Processing cell $cell..."
    local_i, local_j = idx[cell][1], idx[cell][2]
    
    # Gather cell data
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey     = check_predator_has_prey(sp_nm)

    if !predator_has_prey[1]
       local_R -= predator_has_prey[2]
       filter!(name -> !(name in predator_has_prey[3]), sp_nm)
       @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    end
    
    if !isnothing(sp_removed_name) 
        filter!(name -> (name != sp_removed_name), sp_nm)
    end
    
    localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j]) #1000.0
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)
    if !isnothing(NPP)
        localNPP = NPP
    end

    # Attempt setup
    results = attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, sym_competition;
        localNPP      = localNPP,
        localH0_vector= localH0_vector,
        artificial_pi = artificial_pi,
        species_names = sp_nm,
        herbivore_m = herbivore_m,
        predator_m = predator_m
    )
    
    if results === nothing
        @error "Error: results === nothing"
        return nothing
    end
    
    # Destructure the results
    (S2, R2, H_i0, m_i, g_i, G, M_modified,
     a_matrix, A, epsilon_vector, m_alpha) = (
        results.S, results.R, results.H_i0, results.m_i,
        results.g_i, results.G, results.M_modified,
        results.a_matrix, results.A, results.epsilon_vector,
        results.m_alpha
    )

    if (S2 + R2) == 0 || R2 > length(H_i0)
        @error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
        return nothing
    end

    # Build initial conditions
    H_init = H_i0
    P_init = H_init[1:R2] ./ 10.0
    u0 = vcat(H_init, P_init)

    params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
    prob = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)

    # Solve the ODE problem with error suppression
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end

    # Skip if the integration didn't complete successfully
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        @error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
        return nothing
    end

    # --- Plotting the Dynamics using Makie ---
    fig = Figure(resolution=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass", title="Dynamics for cell $cell")
    times_combined = sol.t

    # Plot herbivore dynamics (indices 1:S2) as blue solid lines.
    for i in 1:S2
        lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
    end

    # Plot predator dynamics (indices S2+1:S2+R2) as red dashed lines.
    for α in 1:R2
        lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
    end

    if plot
        display(fig)
    end

    # Evaluate outputs
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
    H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
    P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S2 + R2
    survival_rate = total_surv / total_species
    giHi = sum(g_i .* H_end)
    herbivore_survival_rate = survived_herb / S2
    predator_survival_rate = survived_pred / R2
    H_biomass = sum(H_end)
    P_biomass = sum(P_end)
    biomass_at_the_end = H_biomass + P_biomass
    ratio = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
    real_abundances = filter(x -> !iszero(x), DA_birmmals_with_pi_corrected[local_i, local_j].a)
    similarity = cor(real_abundances, vcat(H_end, P_end))

    # Store the results in a DataFrame row
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
        similarity              = similarity, 
        mu                      = mu_val,
        mu_predation            = mu_pred_val,
        epsilon_val             = eps_val,
        symmetrical_competition = sym_competition,
        NPP                     = localNPP,
        g_iH_i                  = giHi,
        g_iH_i_over_NPP         = round(giHi / localNPP, digits=4), 
        survived_herbivores     = survived_herb,
        survived_predators      = survived_pred,
        total_survivors         = total_surv,
        total_species           = total_species,
        herbivore_survival_rate = herbivore_survival_rate,
        predator_survival_rate  = predator_survival_rate,
        H_biomass               = H_biomass,
        P_biomass               = P_biomass,
        biomass_at_the_end      = biomass_at_the_end,
        herb_pred_ratio         = ratio,
        H_vector                = [H_end],
        P_vector                = [P_end]
    )
    
    # Append the results to the DataFrame
    append!(AAAA, single_run_results)

    @info "The survival rate is $(round(survival_rate, digits=4))"
    return AAAA
end

AAAA = single_run(
    1, 0.3, 0.02, 0.85, true;
    plot=true, sp_removed_name=nothing,
    artificial_pi=true, NPP=840.0
)

local_i, local_j = idx[1][1], idx[1][2]
local_NPP = npp_DA_relative_to_1000[local_i, local_j]
# Gather cell data
real_abundances = filter(x -> !iszero(x), DA_birmmals_with_pi_corrected[local_i, local_j].a)
cor_real_pred = round(cor(real_abundances, vcat(AAAA.H_vector[1], AAAA.P_vector[1])), digits=4)

begin
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title="Real abundances", xlabel="Species", ylabel="Abundance")
    barplot!(ax, real_abundances, color = :lightblue)
    ax2 = Axis(fig[1, 2], title="Predicted abundances", xlabel="Species", ylabel="Abundance")
    barplot!(ax2, vcat(AAAA.H_vector[1], AAAA.P_vector[1]), color = :lightblue)
    display(fig)
end
