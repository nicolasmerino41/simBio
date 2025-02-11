# ---------------------------------------------------------------------------
# 1) Load necessary modules and source scripts
# ---------------------------------------------------------------------------
include("Scripts/prior.jl")
include("Scripts/prior2.jl")
include("Scripts/DA_birmmals_with_pi.jl")
include("Scripts/generate_competition_matrix.jl")
include("Scripts/species_dict.jl")

include("Scripts/ecosystem_dynamics!.jl")
include("Scripts/FI_functions.jl")
include("Scripts/extract_H0_DA.jl")
include("Scripts/attempt_setup_community.jl")
# include("Scripts/Callbacks_function.jl")
include("Scripts/npp_DA_relative_to_1000.jl")
include("Scripts/attempt_feasibility.jl")
# ---------------------------------------------------------------------------
# 2) Set parameters from environment (or hardcode for testing)
# In a production OpenMole run these would be provided externally.
# ---------------------------------------------------------------------------
# Define constants
const EXTINCTION_THRESHOLD = 1e-6
const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0
const art_pi              = true

# ---------------------------------------------------------------------------
# 3) Define the cell to process
# For this run, we pick one cell—OpenMole will call this script repeatedly.
# ---------------------------------------------------------------------------
cell = 2
local_i, local_j = idx[cell][1], idx[cell][2]

# Gather cell data
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
local_S, local_R      = identify_n_of_herbs_and_preds(sp_nm)
predator_has_prey     = check_predator_has_prey(sp_nm)

if !predator_has_prey[1]
    local_R -= predator_has_prey[2]
    filter!(name -> !(name in predator_has_prey[3]), sp_nm)
    @info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
    species_names = sp_nm
else
    species_names = sp_nm
end

localNPP              = Float64(npp_DA_relative_to_1000[local_i, local_j]) 
localH0_vector        = Vector{Float64}(H0_DA[local_i, local_j].a)

# Attempt setup
results = attempt_setup_community(
    local_i, local_j,
    mu, mu_predation, epsilon, true;
    localNPP      = localNPP,
    localH0_vector= localH0_vector,
    species_names = sp_nm,
    artificial_pi = art_pi
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
    
    #### THIS APPROACH ELIMINATES THE SOLVER WARNING ####
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    end

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
    
        giHi = sum(g_i .* H_end)
        ratio_ok = (giHi / localNPP > 0.5) && (giHi / localNPP < 1.5)
        if ratio_ok
            return survival_rate, herb_pred_ratio
        else 
            survival_rate = 0.0
            herb_pred_ratio = 0.0
            return survival_rate, herb_pred_ratio
        end
    end

end

# ---------------------------------------------------------------------------
# 5) Final outputs: these variables will be read by OpenMole.
# ---------------------------------------------------------------------------
@info "Final outputs: survival_rate = $survival_rate, herb_pred_ratio = $herb_pred_ratio"