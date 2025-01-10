include("DA_birmmals_with_pi.jl")
include("Functions/generate_competition_matrix.jl")
include("Functions/species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
#= This script will iterate over all iberian cells and save the parameter configuration 
that maximises the number of survived species in the cell. If there is more than one configuration, 
it will save all of them =#

# Define your parameter ranges
mu_vals = range(0.1, 0.9, length=2)
mu_predation_vals = range(0.0, 0.1, length=2)
epsilon_vals = range(0.1, 1.0, length=2)
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

EXTINCTION_THRESHOLD = 1e-6
T_ext = 250.0

Threads.@threads for cell in 1:7
    local_i, local_j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell (i=$local_i, j=$local_j)..."

    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)

    local_NPP      = Float64(npp_DA[local_i, local_j])
    local_H0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

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

    for (p_idx, combo) in enumerate(param_combinations)
        mu_val, mu_pred_val, eps_val, sym_competition = combo
        local_H_i0 = Float64[]
        local_m_i  = Float64[]

        try
            local_S, local_R,
            species_names, herb_list, pred_list,
            H_i0, m_i,
            p_vec, x_final, g_i,
            localHatH, G,
            M_modified, a_matrix, A,
            epsilon_vector, m_alpha = setup_community_from_cell(
                local_i, local_j;
                NPP = local_NPP,
                M_mean = 0.1,
                mu = mu_val,
                symmetrical_competition = sym_competition,
                mean_m_alpha = 0.1,
                epsilon_val = eps_val,
                mu_predation = mu_pred_val,
                iberian_interact_NA = iberian_interact_NA,
                species_dict = species_dict,
                m_standard_deviation = 0.0,
                h_standard_deviation = 0.0,
                artificial_pi = false,
                real_H0 = true,
                H0_vector = local_H0_vector
            )

            # only if successful, we set local_H_i0, local_m_i, etc.
            local local_H_i0 = H_i0
            local local_m_i  = m_i
            local local_herb_list = herb_list
            local local_pred_list = pred_list
            local local_S, local_R = local_S, local_R
            local local_hatH = localHatH,
            local local_g_i = g_i,
            local local_M_modified = M_modified,
            local local_a_matrix = a_matrix,
            local local_A = A,
            local local_epsilon_vector = epsilon_vector,
            local local_m_alpha = m_alpha
            println("Thread $(Threads.threadid()): Configuration $p_idx for cell $cell (i=$local_i, j=$local_j) was successful.")
        catch e
            # If an error occurred inside `setup_community_from_cell`, skip
            println("Thread $(Threads.threadid()): Configuration $p_idx for cell $cell (i=$local_i, j=$local_j) failed with error: $e")
            continue
        end

        # If S=0 (no herbivores), or R=0, can skip or continue as needed
        if local_S + local_R == 0
            println("Thread $(Threads.threadid()): Configuration $p_idx for cell $cell (i=$local_i, j=$local_j) has no herbivores or predators, skipping...")
            continue
        end

        # Construct initial conditions
        H_init = local_H_i0
        println("H_init: ", H_init)
        # If local_R > length(H_i0), handle it or skip
        if local_R > length(H_init)
            println("Thread $(Threads.threadid()): Configuration $p_idx for cell $cell (i=$local_i, j=$local_j) has more predators than herbivores, skipping...")
            println("local_R: $local_R, length(H_init): $(length(H_init))")
            # means more predators than size of H_i0 => skip
            continue
        end
        P_init = H_init[1:local_R] ./ 10.0
        u0 = vcat(H_init, P_init)

        params = (
            local_S, local_R,
            local_H_i0, local_m_i,
            local_g_i, local_G,
            local_M_modified,
            local_a_matrix, local_A,
            local_epsilon_vector, local_m_alpha
        )
        tspan = (0.0, 500.0)

        # Choose callback
        # extinction_trigger = false
        # cb = extinction_trigger ? cb_trigger : cb_no_trigger

        # Solve ODE
        prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
        sol  = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6) #callback=cb, reltol=1e-6, abstol=1e-6)
        println("Break1")
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            println("Thread $(Threads.threadid()): Configuration $p_idx for cell $cell (i=$local_i, j=$local_j) has unstable solution, skipping...")
            continue
        end
        
        # Evaluate survivors
        H_end = sol[1:local_S, end]
        P_end = sol[local_S+1:local_S+local_R, end]

        survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv    = survived_herb + survived_pred
        total_species = local_S + local_R

        survival_rate            = round(total_surv / total_species, digits=2)
        herbivore_survival_rate  = round(survived_herb / local_S, digits=2)
        predator_survival_rate   = round(survived_pred / local_R, digits=2)
        H_biomass                = round(sum(H_end[H_end .> EXTINCTION_THRESHOLD]), digits=2)
        P_biomass                = round(sum(P_end[P_end .> EXTINCTION_THRESHOLD]), digits=2)
        biomass_at_the_end       = round(H_biomass + P_biomass, digits=2)
        herb_pred_ratio          = ifelse(H_biomass == 0.0, NaN, round(P_biomass / H_biomass, digits=5))
        g_iH_i                   = round(sum(g_i .* H_end), digits=0)

        # Threadsafe push! to scan_results
        @lock global_lock begin
            push!(scan_results, (
                mu_val,
                mu_pred_val,
                eps_val,
                sym_competition,
                local_NPP,
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
    end  # for p_idx in eachindex(param_combinations)

    # find best result(s) in this cell
    if !isempty(scan_results)
        max_survivors = maximum(scan_results.total_survivors)
        best_in_cell  = scan_results[scan_results.total_survivors .== max_survivors, :]

        # Append to best_params_all_cells
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
    elseif isempty(scan_results)
        println("No results found for cell $cell.")
    end
end  # end Threads.@threads for cell in 1:7

# Save final data to CSV
CSV.write("best_params_per_cell.csv", best_params_all_cells)
@info "Done! Best parameter configurations for each cell have been saved."
perro = CSV.File("best_params_per_cell.csv") |> DataFrame