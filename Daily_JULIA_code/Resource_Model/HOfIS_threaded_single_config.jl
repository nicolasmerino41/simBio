################################################################################
# Instead of storing the best configurations (max survivors), we stop as soon as
# we find a parameter configuration for which all species survive.
################################################################################
include("DA_birmmals_with_pi.jl")
include("Functions/generate_competition_matrix.jl")
include("Functions/species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
include("Functions/attempt_setup_community.jl")
include("Functions/Callbacks_function.jl")
include("npp_DA_relative_to_1000.jl")
include("Functions/Computing_metrics.jl")
include("../HerbivoresVsPredators/Exploring HerbPred metaweb composition.jl")

global EXTINCTION_THRESHOLD = 1e-6
global T_ext               = 250.0

# const iberian_interact_NA = iberian_interact_NA
# const species_dict = species_dict
# const predator_names = predator_names
# const herbivore_names = herbivore_names
# const spain_names = spain_names

# 1) Define your parameter ranges
mu_vals               = [0.5]
mu_predation_vals     = range(0.0, 0.1, length=100)
epsilon_vals          = range(0.1, 1.0, length=50)
sym_competition_vals  = [true]

# mu_vals = [0.811111111]
# mu_predation_vals = [0.012121212]
# epsilon_vals = [0.8897959]
# sym_competition_vals = [true]

# 2) Build the parameter combinations
param_combinations = [
    (mu, mu_predation, epsilon_val, sym_comp) 
    for mu in mu_vals
    for mu_predation in mu_predation_vals
    for epsilon_val in epsilon_vals
    for sym_comp in sym_competition_vals
]

# Shuffle the order of param_combinations so that we don't always test the same ones first
param_combinations = Random.shuffle!(param_combinations)
param_combinations = param_combinations[1:2000]

# 3) Prepare DataFrame for storing best results
best_params_24_cells_fixed_mu05_2000iter = DataFrame(
    cell_id                 = Int[],
    i                       = Int[],
    j                       = Int[],
    mu                      = Float64[],
    mu_predation            = Float64[],
    epsilon_val             = Float64[],
    symmetrical_competition = Bool[],
    NPP                     = Float64[],
    g_iH_i                  = Float64[],
    g_iH_i_over_NPP         = Float64[],
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

# Constants
MAX_ITERS                  = 100      # Up to 2000 combos
SURVIVAL_THRESHOLD         = 0.0      # We'll store the best if it reaches at least 0.75

@time Threads.@threads for cell in 1:24  # or whichever cells you want
    local_i, local_j = idx[cell][1], idx[cell][2]
    @info "Processing cell $cell (i=$local_i, j=$local_j)..."

    println("Thread $(Threads.threadid()) is processing cell $cell")

    # Gather cell data
    sp_nm                 = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
    local_S, local_R      = identify_n_of_herbs_and_preds(sp_nm)
    predator_has_prey     = check_predator_has_prey(sp_nm)

    if !predator_has_prey[1]
       local_R -= predator_has_prey[2]
       filter!(name -> !(name in predator_has_prey[3]), sp_nm)
       @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
       if predator_has_prey[2] > 0
        species_names = sp_nm
       else
        species_names = nothing # If we can find all predators, we pass species_names = nothing cause setup_community_from_cell will generate it
       end
    end
    
    localNPP              = Float64(npp_DA_relative_to_1000[local_i, local_j]) 
    localH0_vector        = Vector{Float64}(H0_DA[local_i, local_j].a)
    # cb_no_trigger, cb_trg = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)

    # We'll track the best result so far in local variables
    best_survival_rate  = 0.0
    best_result         = nothing
    found_full_survival = false

    # 4) Up to 2000 combos (or fewer, if param_combinations < 2000)
    for (p_idx, combo) in enumerate(param_combinations[1:min(end,MAX_ITERS)])
        mu_val, mu_pred_val, eps_val, sym_competition = combo
        
        if p_idx % 10 == 0
            println("Thread $(Threads.threadid()) is processing cell $cell, combo $p_idx")
        end  
        # Attempt setup
        results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_pred_val, eps_val, sym_competition;
            localNPP      = localNPP,
            localH0_vector= localH0_vector,
            species_names = species_names
        )
        if results === nothing
            continue
        end

        # Destructure
        (S2, R2, H_i0, m_i, g_i, G, M_modified,
         a_matrix, A, epsilon_vector, m_alpha) = (
            results.S, results.R, results.H_i0, results.m_i,
            results.g_i, results.G, results.M_modified,
            results.a_matrix, results.A, results.epsilon_vector,
            results.m_alpha
        )

        if (S2 + R2) == 0 || R2 > length(H_i0)
            continue
        end

        # Build initial conditions
        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0     = vcat(H_init, P_init)

        params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
        prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
        
        #### THIS APPROACH ELIMINATES THE WARNING ####
        # Create a logger that only shows errors (i.e. skips warnings)
        logger = SimpleLogger(stderr, Logging.Error)

        # Wrap the solver call in the logging context:
        sol = with_logger(logger) do
            solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
        ###############################################

        # If it didn't integrate to 500 or is inf/nan, skip
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            continue
        end

        # Evaluate
        H_end          = sol[1:S2, end]
        P_end          = sol[S2+1:S2+R2, end]
        survived_herb  = count(H_end .> EXTINCTION_THRESHOLD)
        survived_pred  = count(P_end .> EXTINCTION_THRESHOLD)
        total_surv     = survived_herb + survived_pred
        total_species  = S2 + R2
        survival_rate  = total_surv / total_species
        giHi           = sum(g_i .* H_end)
        ratio_ok       = (giHi / localNPP > 0.5) && (giHi / localNPP < 10.0)
        

        # If better than what's stored, update best result
        if survival_rate > best_survival_rate && ratio_ok
            # compute the rest:
            herbivore_survival_rate = survived_herb / S2
            predator_survival_rate  = survived_pred / R2
            H_biomass               = sum(H_end)
            P_biomass               = sum(P_end)
            biomass_at_the_end      = H_biomass + P_biomass
            ratio                   = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
            
            best_survival_rate = survival_rate
            best_result = (
                cell_id                 = cell,
                i                       = local_i,
                j                       = local_j,
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
                survival_rate           = survival_rate,
                herbivore_survival_rate = herbivore_survival_rate,
                predator_survival_rate  = predator_survival_rate,
                H_biomass               = H_biomass,
                P_biomass               = P_biomass,
                biomass_at_the_end      = biomass_at_the_end,
                herb_pred_ratio         = ratio
            )
        end

        # Early stop if full survival + the ratio check
        if isapprox(survival_rate, 1.0; atol=1e-10)
            ratio_ok = (giHi / localNPP > 0.5) && (giHi / localNPP < 10.0)
            if ratio_ok
                @info "Cell $cell => full survival with (mu=$mu_val, mu_pred=$mu_pred_val, eps=$eps_val). Stopping early."
                found_full_survival = true
                break
            end
        end
    end

    # after combos (or if we broke early)
    if found_full_survival
        # best_result *should* be 1.0
        @lock global_lock begin
            push!(best_params_24_cells_fixed_mu05_2000iter, best_result)
        end
    else
        # If not found full survival, check if best >= 0.75
        if best_survival_rate >= SURVIVAL_THRESHOLD  # e.g. 0.75
            @lock global_lock begin
                push!(best_params_24_cells_fixed_mu05_2000iter, best_result)
            end
        else
            @info "Cell $cell: best survival was $(round(best_survival_rate, digits=2)) < $SURVIVAL_THRESHOLD => Not storing."
        end
    end
end  # end Threads.@threads

# done
CSV.write("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/best_params_24_cells_fixed_mu05_2000iter.csv", best_params_24_cells_fixed_mu05_2000iter)
@info "Finished. We stored combos with survival >= 0.75, or full survival."
best_params_24_cells_fixed_mu05 = CSV.read("best_results_per_cell_limited_to_2000.csv", DataFrame)

# Plot total biomass vs mu
fig = Figure(resolution = (600, 500))
ax = Axis(fig[1, 1], xlabel = "mu", ylabel = "Total Biomass", title = "Total Biomass vs mu")

# Get unique mu values and assign colors
unique_mu = sort(unique(best_params_all_cells_limited_to_2000.mu))
num_colors_mu = length(unique_mu)

# Plot data and regression line for each mu value
for (i, mu_value) in enumerate(unique_mu)
    mask = best_params_all_cells_limited_to_2000.mu .== mu_value
    mu_data = best_params_all_cells_limited_to_2000.mu[mask]
    biomass_data = best_params_all_cells_limited_to_2000.biomass_at_the_end[mask]

    # Scatter plot
    MK.scatter!(
        ax, mu_data, biomass_data,
        label = "mu = $(mu_value)", markersize = 10
    )
end
display(fig)
