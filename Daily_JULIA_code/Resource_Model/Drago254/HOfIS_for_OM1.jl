# ---------------------------------------------------------------------------
# 1) Load necessary modules and source scripts
# ---------------------------------------------------------------------------
include("DFs/prior.jl")                  # If needed
include("DFs/prior2.jl")                 # If needed
include("DFs/DA_birmmals_with_pi.jl")
include("DFs/generate_competition_matrix.jl")
include("DFs/species_dict.jl")

include("DFs/ecosystem_dynamics!.jl")
include("DFs/FI_functions.jl")
include("DFs/extract_H0_DA.jl")
include("DFs/attempt_setup_community.jl")
include("DFs/npp_DA_relative_to_1000.jl")
# ---------------------------------------------------------------------------
# 2) Set parameters from environment (or hardcode for testing)
# In a production OpenMole run these would be provided externally.
# ---------------------------------------------------------------------------
# Other fixed parameters
M_mean               = 0.1
mean_m_alpha         = 0.1
m_standard_deviation = 0.0
h_standard_deviation = 0.0
artificial_pi        = false
real_H0              = true

# Define constants
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 50000      # Up to 2000 combos
const SURVIVAL_THRESHOLD  = 0.0       # For example, store best if survival rate >= threshold

# ---------------------------------------------------------------------------
# 3) Define the cell to process
# For this run, we pick one cell—OpenMole will call this script repeatedly.
# ---------------------------------------------------------------------------
cell = 2
local_i, local_j = idx[cell][1], idx[cell][2]

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

# Attempt setup
results = attempt_setup_community(
    local_i, local_j,
    mu, mu_predation, epsilon, true;
    localNPP      = localNPP,
    localH0_vector= localH0_vector,
    species_names = sp_nm
)

if results === nothing
    survival_rate = 0.0
    herb_pred_ratio = 0.0
else
    # Destructure
    (S2, R2, H_i0, m_i, g_i, G, M_modified,
    a_matrix, A, epsilon_vector, m_alpha) = (
       results.S, results.R, results.H_i0, results.m_i,
       results.g_i, results.G, results.M_modified,
       results.a_matrix, results.A, results.epsilon_vector,
       results.m_alpha
   )
   
    # Build initial conditions
    H_init = H_i0
    P_init = H_init[1:R2] ./ 10.0
    u0     = vcat(H_init, P_init)

    params = (S2, R2, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
    prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), params)
    sol  = solve(prob, Tsit5();
                #  callback = cb_no_trigger, # or cb_trigger if you want forced extinction
                reltol=1e-6, abstol=1e-6
            )

    # If it didn't integrate to 500 or is inf/nan, skip
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        survival_rate = 0.0
        herb_pred_ratio = 0.0
        return survival_rate, herb_pred_ratio
    else
    # Evaluate
    EXTINCTION_THRESHOLD = 1e-6
    H_end          = sol[1:S2, end]
    P_end          = sol[S2+1:S2+R2, end]
    survived_herb  = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred  = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv     = survived_herb + survived_pred
    total_species  = S2 + R2
    survival_rate  = total_surv / total_species
    P_biomass      = sum(P_end[P_end .> EXTINCTION_THRESHOLD])
    H_biomass      = sum(H_end[H_end .> EXTINCTION_THRESHOLD])
    herb_pred_ratio = (H_biomass == 0.0) ? NaN : (P_biomass / H_biomass)
    
    return survival_rate, herb_pred_ratio
    end

end

# ---------------------------------------------------------------------------
# 5) Final outputs: these variables will be read by OpenMole.
# ---------------------------------------------------------------------------
@info "Final outputs: survival_rate = $survival_rate, herb_pred_ratio = $herb_pred_ratio"