include("DA_birmmals_with_pi.jl")
include("Functions/generate_competition_matrix.jl")
include("Functions/species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
include("Functions/attempt_setup_community.jl")
include("Functions/Callbacks_function.jl")
#= This script will iterate over all iberian cells and save the parameter configuration 
that maximises the number of survived species in the cell. If there is more than one configuration, 
it will save all of them =#

# Define your parameter ranges
mu_vals = range(0.1, 0.9, length=10)
mu_predation_vals = range(0.0, 0.1, length=100)
epsilon_vals = range(0.1, 1.0, length=50)
sym_competition_vals = [true]
EXTINCTION_THRESHOLD = 1e-6
T_ext = 250.0

# Build all param combos
param_combinations = [
    (mu, mu_predation, epsilon_val, sym_comp)
    for mu in mu_vals
    for mu_predation in mu_predation_vals
    for epsilon_val in epsilon_vals
    for sym_comp in sym_competition_vals
]

# Prepare global DataFrame
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

global_lock = ReentrantLock()

global EXTINCTION_THRESHOLD = 1e-6
global T_ext = 250.0

for cell in 2
    local_i, local_j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell (i=$local_i, j=$local_j)..."

    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

    localNPP      = 1000.0 #Float64(npp_DA[local_i, local_j])
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # Prepare a DataFrame for this cell
    scan_results = DataFrame(
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

    cb_no_trigger, cb_trigger = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)

    Threads.@threads for combo in param_combinations
        # println("vaig fent")
        mu_val, mu_pred_val, eps_val, sym_competition = combo

        # Try to set up the community
        results = attempt_setup_community(
            local_i, local_j, mu_val, mu_pred_val, eps_val, sym_competition;
            localNPP = localNPP,
            localH0_vector = localH0_vector
        )
        # If nothing, skip
        if results === nothing
            continue
        end

        # Otherwise destructure the NamedTuple
        (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha) = 
            (results.S, results.R, results.H_i0, results.m_i, 
             results.g_i, results.G, results.M_modified, 
             results.a_matrix, results.A, results.epsilon_vector, 
             results.m_alpha)

        # If S=0 or R=0 skip
        if S2 + R2 == 0
            continue
        end

        if R2 > length(H_i0)
            # skip, because P_init from H_i0[1:R2] won't work
            continue
        end

        # Build initial conditions
        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0 = vcat(H_init, P_init)

        params = (
            S2, R2,
            H_i0, m_i, 
            g_i, G,
            M_modified, a_matrix, A, 
            epsilon_vector, m_alpha
        )

        prob = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
        sol  = solve(prob, Tsit5(); callback=cb_no_trigger, reltol=1e-6, abstol=1e-6)

        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            continue
        end

        # Evaluate results
        H_end = sol[1:S2, end]
        P_end = sol[S2+1:S2+R2, end]

        survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv    = survived_herb + survived_pred
        total_species = S2 + R2

        survival_rate            = round(total_surv / total_species, digits=2)
        herbivore_survival_rate  = round(survived_herb / S2, digits=2)
        predator_survival_rate   = round(survived_pred / R2, digits=2)
        H_biomass                = round(sum(H_end[H_end .> EXTINCTION_THRESHOLD]), digits=2)
        P_biomass                = round(sum(P_end[P_end .> EXTINCTION_THRESHOLD]), digits=2)
        biomass_at_the_end       = round(H_biomass + P_biomass, digits=2)
        herb_pred_ratio          = ifelse(H_biomass == 0.0, NaN, round(P_biomass / H_biomass, digits=5))
        g_iH_i                   = round(sum(g_i .* H_end), digits=0)

        # push! to local scan_results
        @lock global_lock begin
            push!(scan_results, (
                mu_val,
                mu_pred_val,
                eps_val,
                sym_competition,
                localNPP,
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
    end  # end param_combinations

    # If we got any solutions, pick best, push! to best_params_all_cells
    if !isempty(scan_results)
        max_survivors = maximum(scan_results.total_survivors)
        best_in_cell  = scan_results[scan_results.total_survivors .== max_survivors, :]

        @lock global_lock begin
            for row in eachrow(best_in_cell)
                push!(best_params_all_cells, (
                    cell_id                 = cell,
                    i                       = local_i,
                    j                       = local_j,
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
    else
        println("No results found for cell $cell.")
    end
end  # end @threads

# Save final
# CSV.write("best_params_per_cell_9to50.csv", best_params_all_cells)
# @info "Done! Best parameter configurations for each cell have been saved."

# perro = CSV.File("best_params_per_cell.csv") |> DataFrame