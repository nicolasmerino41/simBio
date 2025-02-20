# Define your parameter ranges
mu_vals               = 0.2777777778
mu_predation_vals     = 0.02713178
epsilon_vals          = 0.81052
sym_competition_vals  = [true]

cell = 2

# Build a single parameter combination
param_combinations = [
    (mu, mu_predation, epsilon_val, sym_comp) 
    for mu in mu_vals
    for mu_predation in mu_predation_vals
    for epsilon_val in epsilon_vals
    for sym_comp in sym_competition_vals
]

# Define constants
# const EXTINCTION_THRESHOLD = 1e-6
# const T_ext               = 250.0
const MAX_ITERS           = 2000
const SURVIVAL_THRESHOLD  = 0.0

# Placeholder for results
AAAA = DataFrame()

# Run the model for the single parameter configuration
for (mu_val, mu_pred_val, eps_val, sym_competition) in param_combinations
    
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
    
    localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j]) #1000.0
    localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

    # Attempt setup
    results = attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, sym_competition;
        localNPP      = localNPP,
        localH0_vector= localH0_vector,
        species_names = sp_nm
    )
    
    if results === nothing
        continue
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
        continue
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
        continue
    end

    # Evaluate outputs
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
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

    # Store the results in a DataFrame row
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
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
        herb_pred_ratio         = ratio
    )
    
    # Append the results to the DataFrame
    append!(AAAA, single_run_results)
end

# Print the results DataFrame
println(AAAA)
println(AAAA.survival_rate)

#################### SAME THING BUT AS A FUNCTION ######################
function single_run(cell, mu_val, mu_pred_val, eps_val, sym_competition; sp_removed_name=nothing, artificial_pi=false, NPP=nothing)
    
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
        species_names = sp_nm,
        artificial_pi = artificial_pi
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

    # Evaluate outputs
    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]
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

    # Store the results in a DataFrame row
    single_run_results = DataFrame(
        cell_id                 = cell,
        i                       = local_i,
        j                       = local_j,
        survival_rate           = survival_rate,
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
        herb_pred_ratio         = ratio
    )
    
    # Append the results to the DataFrame
    append!(AAAA, single_run_results)

    @info "The survival rate is $(round(survival_rate, digits=4))"
    return AAAA
end
# single_run(1, 0.9, 0.16, 1.0, true; sp_removed_name=nothing, artificial_pi=true, NPP=nothing)

function new_single_run(cell, mu_val, mu_pred_val, eps_val, sym_competition; 
    sp_removed_name=nothing, artificial_pi=false, NPP=nothing)

# Placeholder for results
AAAA = DataFrame()

@info "Processing cell $cell..."
local_i, local_j = idx[cell][1], idx[cell][2]

# Gather cell data
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
predator_has_prey = check_predator_has_prey(sp_nm)

if !predator_has_prey[1]
local_R -= predator_has_prey[2]
filter!(name -> !(name in predator_has_prey[3]), sp_nm)
@info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
end

if !isnothing(sp_removed_name) 
filter!(name -> (name != sp_removed_name), sp_nm)
end

localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)
if !isnothing(NPP)
localNPP = NPP
end

# Attempt to set up the community using the new parametrisation
results = new_attempt_setup_community(
local_i, local_j,
mu_val, mu_pred_val, eps_val, sym_competition;
localNPP      = localNPP,
localH0_vector= localH0_vector,
species_names = sp_nm,
artificial_pi = artificial_pi
)

if results === nothing
@error "Error: results === nothing"
return nothing
end

# Destructure the results.
# The new attempt_setup_community returns a NamedTuple with:
# S, R, species_names, herbivore_list, predator_list,
# H_i0, m_i, g_i, x_final, beta, G, M_modified, a_matrix, A, epsilon, m_alpha, raw_g.
S2           = results.S
R2           = results.R
H_i0         = results.H_i0
m_i          = results.m_i
g_i          = results.g_i
x_final      = results.x_final   # scaling parameter (if you need it later)
beta         = results.beta      # niche parameters
G            = results.G
M_modified   = results.M_modified
a_matrix     = results.a_matrix
A            = results.A
epsilon_vec  = results.epsilon   # renamed for clarity
m_alpha      = results.m_alpha
# raw_g is also returned, if needed: raw_g = results.raw_g

if (S2 + R2) == 0 || R2 > length(H_i0)
@error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
return nothing
end

# Build initial conditions.
# Herbivore initial conditions come from the empirical abundances (H_i0).
# For predators, we initialize them as a fraction of the herbivore abundances.
H_init = H_i0
P_init = H_init[1:R2] ./ 10.0
u0 = vcat(H_init, P_init)

# Package parameters for the new ODE function.
# Note: new_dynamics! expects parameters in the following order:
# (S, R, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon, m_alpha)
params = (S2, R2, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon_vec, m_alpha)

# Set up the ODE problem using the new dynamics function.
prob = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params)

# Solve the ODE problem with error suppression
logger = SimpleLogger(stderr, Logging.Error)
sol = with_logger(logger) do
solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
end

# Skip if the integration did not complete successfully.
if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
@error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
return nothing
end

# Evaluate outputs.
H_end = sol[1:S2, end]
P_end = sol[S2+1:S2+R2, end]
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

# Store the results in a DataFrame row.
single_run_results = DataFrame(
cell_id                 = cell,
i                       = local_i,
j                       = local_j,
survival_rate           = survival_rate,
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
herb_pred_ratio         = ratio
)

# Append the results to the DataFrame.
append!(AAAA, single_run_results)

@info "The survival rate is $(round(survival_rate, digits=4))"
return AAAA
end

# new_single_run(1, 0.5, 0.01, 0.8, true; sp_removed_name=nothing, NPP=nothing)

function new_single_run_with_plot(cell, mu_val, mu_pred_val, eps_val, sym_competition; 
    sp_removed_name=nothing, artificial_pi=false, NPP=nothing)
# Placeholder for results
AAAA = DataFrame()

@info "Processing cell $cell..."
local_i, local_j = idx[cell][1], idx[cell][2]

# Gather cell data
sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
predator_has_prey = check_predator_has_prey(sp_nm)
println("here S = $local_S, R = $local_R")
if !predator_has_prey[1]
local_R -= predator_has_prey[2]
filter!(name -> !(name in predator_has_prey[3]), sp_nm)
@info("In cell $cell, we removed $(predator_has_prey[2]) predators: $(predator_has_prey[3]).")
end

if !isnothing(sp_removed_name)
filter!(name -> (name != sp_removed_name), sp_nm)
end

localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)
if !isnothing(NPP)
localNPP = NPP
end

# Attempt setup using the updated parametrisation
results = new_attempt_setup_community(
local_i, local_j,
mu_val, mu_pred_val, eps_val, sym_competition;
localNPP       = localNPP,
localH0_vector = localH0_vector,
species_names  = sp_nm,
artificial_pi  = artificial_pi
)

if results === nothing
@error "Error: results === nothing"
return nothing
end

# Destructure the returned NamedTuple.
S2          = results.S
R2          = results.R
H_i0        = results.H_i0
m_i         = results.m_i
g_i         = results.g_i
x_final     = results.x_final   # scaling parameter (if needed)
beta        = results.beta      # niche parameters
G           = results.G
M_modified  = results.M_modified
a_matrix    = results.a_matrix
A           = results.A
epsilon_vec = results.epsilon   # predator conversion vector
m_alpha     = results.m_alpha
println("here S = $S2, R = $R2")
if (S2 + R2) == 0 || R2 > length(H_i0)
@error "Error: (S2 + R2) == 0 || R2 > length(H_i0)"
return nothing
end

# Build initial conditions.
# Herbivore initial conditions are set to the empirical abundances (H_i0).
# Predator initial conditions are initialized as a fraction of the herbivore abundances.
H_init = H_i0
P_init = H_init[1:R2] ./ 10.0
u0 = vcat(H_init, P_init)

# Package parameters for new_dynamics!
# new_dynamics! expects: (S, R, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon, m_alpha)
params = (S2, R2, H_i0, m_i, g_i, beta, G, M_modified, a_matrix, A, epsilon_vec, m_alpha)
prob = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params)

# Solve the ODE problem with error suppression
logger = SimpleLogger(stderr, Logging.Error)
sol = with_logger(logger) do
solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
end

if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
@error "Error: sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])"
return nothing
end

# --- Plotting the Dynamics using Makie ---
# Create a figure with the specified resolution.
fig = Figure(resolution=(1200, 600))
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Biomass", 
title="Dynamics for cell $cell")

# Extract the time vector from the solution.
times_combined = sol.t

# Plot herbivores (indices 1:S2) as blue solid lines.
for i in 1:S2
lines!(ax, times_combined, sol[i, :], label="H$(i)", color=:blue)
end

# Plot predators (indices S2+1:S2+R2) as red dashed lines.
for α in 1:R2
lines!(ax, times_combined, sol[S2+α, :], label="P$(α)", linestyle=:dash, color=:red)
end

# Add a legend to the plot.
# axislegend(ax, position=:rt)
display(fig)

# --- Evaluate Outputs ---
H_end = sol[1:S2, end]
P_end = sol[S2+1:S2+R2, end]
Hobs = H_i0
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

# Store the results in a DataFrame row.
single_run_results = DataFrame(
cell_id                 = cell,
i                       = local_i,
j                       = local_j,
survival_rate           = survival_rate,
mu                      = mu_val,
mu_predation            = mu_pred_val,
epsilon_val             = eps_val,
symmetrical_competition = sym_competition,
NPP                     = localNPP,
g_iH_i                  = giHi,
g_iH_i_over_NPP         = round(giHi / localNPP, digits=4),
g_iHobs                = round(sum(g_i .* Hobs), digits=4),
survived_herbivores     = survived_herb,
survived_predators      = survived_pred,
total_survivors         = total_surv,
total_species           = total_species,
herbivore_survival_rate = herbivore_survival_rate,
predator_survival_rate  = predator_survival_rate,
H_biomass               = H_biomass,
P_biomass               = P_biomass,
biomass_at_the_end      = biomass_at_the_end,
herb_pred_ratio         = ratio
)

append!(AAAA, single_run_results)
@info "The survival rate is $(round(survival_rate, digits=4))"

# Return both the DataFrame with summary results and the plot figure.
return AAAA, fig
end

# # old_DA_birmmals_with_pi_corrected = deepcopy(DA_birmmals_with_pi_corrected)
# DA_birmmals_with_pi_corrected[18,1] = MyBirmmals(SVector{205, Float64}(rand(205)))
# DA_birmmals_with_pi_corrected = deepcopy(old_DA_birmmals_with_pi_corrected)

# # Call the function with the specified arguments.
# cucu, cucut = new_single_run_with_plot(1, 0.3, 0.03, 0.1, true; sp_removed_name=nothing, NPP=nothing, artificial_pi=false)

# println(cucu)