
include("prior.jl")
include("prior2.jl")
include("DA_birmmals_with_pi.jl")
include("generate_competition_matrix.jl")
include("species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
#= This script will iterate over all iberian cells and save the parameter configuration 
that maximises the number of survived species in the cell. If there is more than one configuration, 
it will save all of them =#

# Define your parameter ranges
mu_vals = range(0.1, 0.9, length=2)
mu_predation_vals = range(0.0, 0.1, length=100)
epsilon_vals = range(0.1, 1.0, length=2)
sym_competition_vals = [true, false]
EXTINCTION_THRESHOLD = 1e-6
T_ext = 250.0

# Build all parameter combinations
param_combinations = [
    (mu, mu_predation, epsilon_val, sym_comp)
    for mu in mu_vals
    for mu_predation in mu_predation_vals
    for epsilon_val in epsilon_vals
    for sym_comp in sym_competition_vals
]

# Prepare a DataFrame to store the best configurations for each cell
best_params_all_cells = DataFrame(
    cell_id                 = Int[],
    i                       = Int[],
    j                       = Int[],
    mu                      = Float64[],
    mu_predation            = Float64[],
    epsilon_val             = Float64[],
    symmetrical_competition = Bool[],
    NPP                     = Float64[],
    g_iH_i                  = Float64[],
    survived_herbivores     = Int[],
    survived_predators      = Int[],
    total_survivors         = Int[],
    total_species           = Int[],
    survival_rate           = Float64[],
    herbivore_survival_rate = Float64[],
    predator_survival_rate  = Float64[],
    H_biomass               = Float64[],
    P_biomass               = Float64[],
    biomass_at_the_end      = Float64[],
    herb_pred_ratio         = Float64[]
)

# A global lock for pushing results in a threadsafe way
global_lock = ReentrantLock()

# Loop over all Iberian cells
@time for cell in 1
    i, j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell (i=$i, j=$j)..."
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[i, j])
    global S, R = identify_n_of_herbs_and_preds(sp_nm)
    # Include expensive callbacks once per cell
    include("Callbacks.jl")
    println("length callbacks = ", length(callbacks))
    # Gather cell-specific data
    NPP = Float64(npp_DA[i, j])  # or a fixed test value
    H0_vector = Vector{Float64}(H0_DA[i, j].a)

    # We collect results for *all* combinations in this cell
    scan_results = DataFrame(
        mu                      = Float64[0],
        mu_predation            = Float64[0],
        epsilon_val             = Float64[0],
        symmetrical_competition = Bool[0],
        NPP                     = Float64[0],
        g_iH_i                  = Float64[0],
        survived_herbivores     = Int[0],
        survived_predators      = Int[0],
        total_survivors         = Int[0],
        total_species           = Int[0],
        survival_rate           = Float64[0],
        herbivore_survival_rate = Float64[0],
        predator_survival_rate  = Float64[0],
        H_biomass               = Float64[0],
        P_biomass               = Float64[0],
        biomass_at_the_end      = Float64[0],
        herb_pred_ratio         = Float64[0]
    )

    # Multi-threaded loop over parameter combinations
    for p_idx in 1:length(param_combinations)
        mu, mu_predation, epsilon_val, sym_competition = param_combinations[p_idx]

        # Wrap in try/catch so we skip if there's an error (e.g. "No real solution")
        try
            global (S, R, species_names, herbivore_list, predator_list,
             H_i0, m_i, p_vec, x_final, g_i, localHatH, G,
             M_modified, a_matrix, A, epsilon_vector, m_alpha) =
                setup_community_from_cell(
                    i, j;
                    NPP = NPP,
                    M_mean = 0.1,
                    mu = mu,
                    symmetrical_competition = sym_competition,
                    mean_m_alpha = 0.1,
                    epsilon_val = epsilon_val,
                    mu_predation = mu_predation,
                    iberian_interact_NA = iberian_interact_NA,
                    species_dict = species_dict,
                    m_standard_deviation = 0.0,
                    h_standard_deviation = 0.0,
                    artificial_pi = false,
                    real_H0 = true,
                    H0_vector = H0_vector,
                )
                # println("S = $S, R = $R, length(H0_vector) = ", length(H_i0))
                # println("u0 length = ", length(u0))
                # println("epsilon_val = ", epsilon_val)

        catch e
            @warn "Cell $cell, mu=$mu, mu_predation=$mu_predation, epsilon=$epsilon_val, sym_comp=$sym_competition failed: $e. Skipping combo."
            continue
        end

        println("The combination ",  param_combinations[p_idx] , " seems to work")
        # Solve the ODE
        H_init = H_i0
        P_init = H_i0[1:R] ./ 10.0
        u0 = vcat(H_init, P_init)

        params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
        tspan = (0.0, 500.0)

        extinction_trigger = false
        cb = extinction_trigger ? cb_trigger : cb_no_trigger

        prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
        sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

        # Validate the solution
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            @warn "Unstable solution for Cell $cell, parameters $param_combinations[p_idx]. Skipping."
            continue
        end

        # Count survivors
        EXTINCTION_THRESHOLD = 1e-6
        H_end = sol[1:S, end]
        P_end = sol[S+1:S+R, end]
        survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv = survived_herb + survived_pred
        total_species = S + R
        survival_rate = total_surv / total_species
        herbivore_survival_rate = survived_herb / S
        predator_survival_rate = survived_pred / R
        H_biomass = sum(H_end[H_end .> EXTINCTION_THRESHOLD])
        P_biomass = sum(P_end[P_end .> EXTINCTION_THRESHOLD])
        biomass_at_the_end = H_biomass + P_biomass
        herb_pred_ratio = P_biomass/H_biomass
        g_iH_i = sum(g_i .* H_end)

        # Thread-safe insertion into scan_results
        @lock global_lock begin
            push!(scan_results, (
                mu,
                mu_predation,
                mean(epsilon_vector), #TODO improve this fix
                sym_competition,
                NPP,
                g_iH_i,
                survived_herb,
                survived_pred,
                total_surv,
                total_species,
                survival_rate,
                herbivore_survival_rate,
                predator_survival_rate,
                H_biomass,
                P_biomass,
                biomass_at_the_end,
                herb_pred_ratio
            ))
        end
    end  # end threaded loop

    # Now find the best result(s) in this cell
    max_survivors = maximum(scan_results.total_survivors)
    best_in_cell = scan_results[scan_results.total_survivors .== max_survivors, :]

    # Append to the global DataFrame
    @lock global_lock begin
        for row in eachrow(best_in_cell)
            push!(best_params_all_cells, (
                cell_id                 = cell,
                i                       = i,
                j                       = j,
                mu                      = row.mu,
                mu_predation            = row.mu_predation,
                epsilon_val             = row.epsilon_val,
                symmetrical_competition = row.symmetrical_competition,
                NPP                     = row.NPP,
                g_iH_i                  = row.g_iH_i,
                survived_herbivores     = row.survived_herbivores,
                survived_predators      = row.survived_predators,
                total_survivors         = row.total_survivors,
                total_species           = row.total_species,
                survival_rate           = row.survival_rate,
                herbivore_survival_rate = row.herbivore_survival_rate,
                predator_survival_rate  = row.predator_survival_rate,
                H_biomass               = row.H_biomass,
                P_biomass               = row.P_biomass,
                biomass_at_the_end      = row.biomass_at_the_end,
                herb_pred_ratio         = row.herb_pred_ratio
            ))
        end
    end
end

# Finally, write out the best parameters for each cell
CSV.write("best_params_per_cell.csv", best_params_all_cells)
@info "Done! Best parameter configurations for each cell have been saved."