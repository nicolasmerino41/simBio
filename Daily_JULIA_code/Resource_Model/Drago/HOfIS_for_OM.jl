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
include("attempt_setup_community.jl")
# ---------------------------------------------------------------------------
# 2) Set parameters from environment (or hardcode for testing)
# In a production OpenMole run these would be provided externally.
# ---------------------------------------------------------------------------
mu           = parse(Float64, get(ENV, "MU", "0.5"))
mu_predation = parse(Float64, get(ENV, "MU_PRED", "0.01"))
epsilon_val  = parse(Float64, get(ENV, "EPSILON", "0.8"))
sym_competition = (get(ENV, "SYM_COMP", "true") == "true")

# Other fixed parameters
M_mean               = 0.1
mean_m_alpha         = 0.1
m_standard_deviation = 0.0
h_standard_deviation = 0.0
artificial_pi        = false
real_H0              = true

# ---------------------------------------------------------------------------
# 3) Define the cell to process
# For this run, we pick one cell—OpenMole will call this script repeatedly.
# ---------------------------------------------------------------------------
cell = 1
local_i, local_j = idx[cell][1], idx[cell][2]
localNPP = 1000.0  # You may replace this with real cell data if available.
localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

# ---------------------------------------------------------------------------
# 4) Run the simulation in a try/catch block so that any error results in outputs 0.0
# ---------------------------------------------------------------------------
# Initialize outputs
survival_rate = 0.0
herb_pred_ratio = 0.0

try
    # 4A) Set up the community for the chosen cell:
    results = setup_community_from_cell(local_i, local_j;
                NPP              = localNPP,
                M_mean           = M_mean,
                mu               = mu,
                symmetrical_competition = sym_competition,
                mean_m_alpha     = mean_m_alpha,
                epsilon_val      = epsilon_val,
                mu_predation     = mu_predation,
                iberian_interact_NA = iberian_interact_NA,
                species_dict     = species_dict,
                m_standard_deviation = m_standard_deviation,
                h_standard_deviation = h_standard_deviation,
                artificial_pi    = artificial_pi,
                real_H0          = real_H0,
                H0_vector        = localH0_vector
            )

    # Destructure the returned NamedTuple:
    (S, R, species_names, herbivore_list, predator_list,
     H_i0, m_i, p_vec, x_final, g_i, localHatH,
     G, M_modified, a_matrix, A, epsilon_vector, m_alpha) = results

    # If no species found or there are more predators than available herbivore data, throw an error.
    if (S + R) == 0 || R > length(H_i0)
         error("Inadequate community data in this cell.")
         @info "shit happened"
    end

    # 4B) Construct initial conditions
    H_init = H_i0
    P_init = H_init[1:R] ./ 10.0
    u0 = vcat(H_init, P_init)
    params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon_vector, m_alpha)
    tspan = (0.0, 500.0)
    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # If integration did not reach t = 500 or the solution is unstable, throw an error.
    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
         error("Integration unstable or incomplete")
         @info "shit happened"
    end

    # 4C) Evaluate results: count surviving species and compute biomass ratio.
    EXTINCTION_THRESHOLD = 1e-6
    H_end = sol[1:S, end]
    P_end = sol[S+1:S+R, end]
    survived_herb = count(H_end .> EXTINCTION_THRESHOLD)
    survived_pred = count(P_end .> EXTINCTION_THRESHOLD)
    total_surv = survived_herb + survived_pred
    total_species = S + R
    global survival_rate = (total_species > 0) ? total_surv / total_species : 0.0

    H_biomass = sum(H_end[H_end .> EXTINCTION_THRESHOLD])
    P_biomass = sum(P_end[P_end .> EXTINCTION_THRESHOLD])
    herb_pred_ratio = (H_biomass > 0) ? (P_biomass / H_biomass) : 0.0
    @info "hey, the try worked"
catch e
    @error "Error in simulation: $e. Returning 0.0 outputs."
    @info "shit happened"
    global survival_rate = 0.0
    global herb_pred_ratio = 0.0
end

# ---------------------------------------------------------------------------
# 5) Final outputs: these variables will be read by OpenMole.
# ---------------------------------------------------------------------------
@info "Final outputs: survival_rate = $survival_rate, herb_pred_ratio = $herb_pred_ratio"

# # Optionally, assign them as global variables for OpenMole to retrieve:
# global survival_rate = survival_rate
# global herb_pred_ratio = herb_pred_ratio
