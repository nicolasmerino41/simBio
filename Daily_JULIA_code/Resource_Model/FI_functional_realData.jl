##############################################################################################
##############################################################################################
############################ SAME CODE BUT INSIDE A FUNCTION #################################

# Assuming the following globals:
# - birmmals_biomass_fixed::DataFrame with columns:
#   species (String), mean_density (Float64), bodyMass (Float64), biomass (Float64)
# - herbivore_names::Vector{String} containing names of herbivores
# - predator_names::Vector{String} containing names of predators
# - DA_birmmals::Matrix{MyBirmmals} containing the DA cells
# - MyBirmmals is a struct with a field 'a' that indicates presence/absence (or abundance)
#   of each species (length matches rows in birmmals_biomass_fixed)

# Example:
# struct MyBirmmals
#     a::SVector{207,Float64}
# end

predator_names = setdiff(spain_names, herbivore_names)

function extract_species_names_from_a_cell(cell::MyBirmmals)
    names = birmmals_biomass_fixed[:, :species]
    species_names = String[]
    for i in 1:length(cell.a)
        if !iszero(cell.a[i])
            push!(species_names, names[i])
        end
    end
    return species_names
end

function identify_n_of_herbs_and_preds(species_names::Vector{String})
    S = 0
    R = 0
    for name in species_names
        if name in herbivore_names
            S += 1
        elseif name in predator_names
            R += 1
        end
    end
    return S, R
end

function parametrise_the_community(
    species_names::Vector{String};
    NPP::Float64 = 1000.0,
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    asymmetry_competition::Bool = false,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    iberian_interact_NA::NamedMatrix{Float64}=iberian_interact_NA,   # Accept NamedMatrix{Float64} here
    species_dict::Dict{String,Int}=species_dict,
    connectivity_hp::Float64 = 1.0,
    m_standard_deviation::Float64 = 0.1,
    h_standard_deviation::Float64 = 10.0
)

    ############################################################################
    # A) Identify which species in cell => which are herbivores vs predators
    ############################################################################
    herbivore_list = [sp for sp in species_names if sp in herbivore_names]
    predator_list  = [sp for sp in species_names if sp in predator_names]
    S = length(herbivore_list)
    R = length(predator_list)

    # We assume you want all herbivores to get H_i^0 = NPP/S, m_i = M_mean
    # and your standard approach to competition among herbivores.
    if S>0
        H0_mean = NPP / S
        # H_i0 = fill(NPP / S, S)  
        # m_i  = fill(M_mean, S)
        H_i0 = [abs(rand(Normal(H0_mean, h_standard_deviation))) for _ in 1:S]
        m_i = [abs(rand(Normal(M_mean, m_standard_deviation))) for _ in 1:S]
    else
        H_i0 = Float64[]
        m_i  = Float64[]
    end

    ############################################################################
    # B) Build the standard mu_matrix among herbivores => M_modified base
    ############################################################################
    if S>0
        # e.g. (V, mu_matrix) = generate_competition_matrix(S, mu, asymmetry_competition)
        V, mu_matrix = generate_competition_matrix(S, mu, asymmetry_competition; check_condition=true)
    else
        # No herbivores => no competition matrix
        V = zeros(0,0)
        mu_matrix = zeros(0,0)
    end

    ############################################################################
    # C) Build submatrix from iberian_interact_NA for predator-predator & 
    #    herb-pred feeding (a_matrix).
    #    1 means "row-species can feed on col-species."
    ############################################################################

    # For each predator sp in predator_list, find row index in iberian_interact_NA
    # For each herbivore in herbivore_list, find row index, etc.
    # We'll define "all_sp = herbivore_list... + predator_list..." in that order
    all_sp = vcat(herbivore_list, predator_list)
    total_in_cell = S + R

    # local_indices: row/col in iberian_interact_NA
    local_indices = [ species_dict[sp] for sp in all_sp ]

    # We can skip building a full submatrix if we just parse: 
    # - (S x R) for herb->pred = a_matrix
    # - (R x R) for pred->pred = A
    # We'll define a_matrix, A:

    a_matrix = zeros(S, R)
    A        = Matrix{Float64}(I, R, R)  # We'll turn diag negative later

    # fill a_matrix: row=herb_i, col=pred_j
    # Suppose herb_i is row i in [1..S], pred_j is row j in [S+1..S+R]
    # in the sub-block. If iberian_interact_NA says "pred_j eats herb_i" => 1
    # => we define a certain feeding strength, e.g. mu_predation or scaled by "connectivity_hp"
    for herb_i in 1:S
        # global index for that herb in iberian_interact_NA
        global_herb_idx = species_dict[herbivore_list[herb_i]]
        for pred_j in 1:R
            global_pred_idx = species_dict[predator_list[pred_j]]
            if iberian_interact_NA[global_pred_idx, global_herb_idx] == 1
                # pred_j can eat herb_i
                # scale by mu_predation & connectivity or define a random approach
                # for demonstration, we just set a_matrix[herb_i, pred_j] = mu_predation
                if rand() < connectivity_hp
                    a_matrix[herb_i, pred_j] = mu_predation
                end
            end
        end
    end

    # fill predator-predator matrix A: row=pred_j, col=pred_k
    # if iberian_interact_NA says "pred_j eats pred_k", we interpret as A[j,k] > 0 => row-pred j 
    # has feeding on col-pred k. Usually for a standard LV approach, we store negative diag for self-limit
    for row_pred in 1:R
        # negative diag
        A[row_pred, row_pred] = -1.0  
        global_pred_row = species_dict[predator_list[row_pred]]
        for col_pred in 1:R
            if row_pred != col_pred
                global_pred_col = species_dict[predator_list[col_pred]]
                if iberian_interact_NA[global_pred_row, global_pred_col] == 1
                    # row_pred eats col_pred => define some + or -? Typically you'd interpret "feeding" 
                    # in a predator-pred matrix as predation, so A[row_pred,col_pred] is beneficial to row_pred.
                    # In many LV setups, we store negative diag, plus 0 or minimal for other off-diag. 
                    # Let's define if row_pred "eats" col_pred => A[row_pred,col_pred] = some pos or just 0. 
                    # As a naive example, set 0.01:
                    A[row_pred, col_pred] = 0.01
                end
            end
        end
    end

    ############################################################################
    # D) Predator mortality, assimilation
    ############################################################################
    m_alpha = fill(mean_m_alpha, R)
    epsilon = fill(epsilon_val, R)

    ############################################################################
    # E) Compute C_{ij}, G_i, then final M_modified among herbivores
    ############################################################################
    # invert A => A_inv
    if R>0
        A_inv = inv(A)
    else
        A_inv = zeros(0,0)
    end

    SxS = S
    C = zeros(SxS,SxS)
    G = zeros(SxS)

    for r in 1:S
        for c in 1:S
            val = 0.0
            for α in 1:R
                for β in 1:R
                    val += epsilon[α] * a_matrix[r, α]*A_inv[α, β]*a_matrix[c, β]
                end
            end
            C[r,c] = val
        end
    end

    for r in 1:S
        val = 0.0
        for β in 1:R
            for α in 1:R
                val += a_matrix[r, α]*A_inv[α, β]*m_alpha[β]
            end
        end
        G[r] = val
    end

    # M_modified = mu_matrix + C_ij * H_i^0 / m_i
    M_modified = copy(mu_matrix)
    for r in 1:S, c in 1:S
        M_modified[r,c] += C[r,c] * (H_i0[r]/m_i[r])
    end

    ############################################################################
    # F) Solve for x from doc formula => final g_i
    ############################################################################
    # We'll replicate your typical quadratic approach with (I + M_modified), etc.
    # or do a simpler approach if you prefer. 
    # The doc approach:
    IplusM = I + M_modified
    IM_inv = inv(IplusM)


    barH = mean(H_i0)
    
    # Now define h_i with length S
    h_i = zeros(S)
    for i in 1:S
        # According to the doc formula: H_i^0 = S*bar(H)*h_i => h_i = H_i^0 / (S*bar(H))
        h_i[i] = H_i0[i] / (S * barH)
    end

    # define A_vec, B_vec
    A_vec = zeros(S)
    B_vec = zeros(S)
    for r in 1:S
        A_val = 0.0
        B_val = 0.0
        for c in 1:S
            A_val += IM_inv[r,c]* ( h_i[c] )  # or p[c], but you're using doc eq
            B_val += IM_inv[r,c]* ( G[c]/m_i[c] - 1.0 )
        end
        A_vec[r] = A_val
        B_vec[r] = B_val
    end

    # Then define NPP eq in terms of x with g_i = x * p_i * m_i. If your doc eq 
    # uses hat_H => p => etc., adapt accordingly. The details can mirror your "virtual code".

    # For brevity, let's do the standard A_coef, B_coef, C_coef approach:
    A_coef = 0.0
    B_coef = 0.0
    C_coef = 0.0

    # let p_i = h_i or we define p from hat_H? 
    # Suppose p_i = hat_H, you have p = hat_H. Then g_i = x p_i m_i
    p_vec = copy(h_i)  # or hat_H

    for r in 1:S
        # x^2 term
        A_coef += p_vec[r]*m_i[r]*H_i0[r]* A_vec[r]
        # x^1 term
        B_coef += p_vec[r]*m_i[r]*H_i0[r]* B_vec[r] + G[r]*H_i0[r]* A_vec[r]
        # x^0
        C_coef += G[r]*H_i0[r]* B_vec[r]
    end
    # subtract NPP from constant
    C_coef -= NPP

    disc = B_coef^2 - 4*A_coef*C_coef
    if disc < 0
        error("No real solution for x found!")
    end
    x_candidates = [(-B_coef + sqrt(disc))/(2A_coef), (-B_coef - sqrt(disc))/(2A_coef)]
    pos_x = filter(x->x>0, x_candidates)
    if isempty(pos_x)
        error("No positive x found!")
    end
    x_final = maximum(pos_x)

    # define g_i for herbivores
    g_i = [x_final * p_vec[r] * m_i[r] for r in 1:S]

    return (
        S, R, # for clarity
        H_i0, m_i,
        h_i,
        x_final, g_i,
        G, # from pred interactions
        M_modified,
        a_matrix,
        A,
        epsilon,
        m_alpha
    )
end


function setup_community_from_cell(i::Int, j::Int; 
    NPP::Float64 = 1000.0,
    M_mean::Float64 = 0.1,
    mu::Float64 = 0.5,
    asymmetry_competition::Bool = false,
    mean_m_alpha::Float64 = 0.1,
    epsilon_val::Float64 = 1.0,
    mu_predation::Float64 = 0.01,
    iberian_interact_NA::NamedMatrix{Float64}=iberian_interact_NA,   # Accept NamedMatrix{Float64} here
    species_dict::Dict{String,Int}=species_dict,
    connectivity_hp::Float64 = 1.0,
    m_standard_deviation::Float64 = 0.1,
    h_standard_deviation::Float64 = 10.0
    )

    # 1) Extract species from cell
    cell = DA_birmmals[i,j]
    species_names = extract_species_names_from_a_cell(cell)
    # 2) Identify S,R
    S, R = identify_n_of_herbs_and_preds(species_names)
    # 3) Parametrise community
    S, R, H_i0, m_i, h_i, x_final, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = parametrise_the_community(species_names; 
        NPP, M_mean, mu, asymmetry_competition, mean_m_alpha, epsilon_val, mu_predation, 
        iberian_interact_NA, species_dict, connectivity_hp, m_standard_deviation, h_standard_deviation)

    # Return all parameters, so the calling code can build the ODE problem
    return S, R, species_names, H_i0, m_i, h_i, x_final, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha
end

# With these functions, you can now pick any cell (i,j) from DA_birmmals and quickly build the model parameters:
# (S, R, names, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha) = setup_community_from_cell(20,20)
CELL = idx[20]
spu = extract_species_names_from_a_cell(DA_birmmals[CELL])
S, R = identify_n_of_herbs_and_preds(spu)
H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = parametrise_the_community()

##############################################################################################
# Now we apply the dynamics as done before, but using the parameters from a chosen cell.
##############################################################################################
begin
    
# Choose a cell, for example (20,20)
i, j = 22, 3
S, R, species_names, H_i0, m_i, h_i, x_final, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = setup_community_from_cell(i, j; epsilon_val=0.01)

# Here you can define NPP and other parameters if needed, or assume they are global or computed before
NPP = 1000.0

# ODE definition:
function ecosystem_dynamics!(du, u, p, t)
    S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = p

    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics
    for i in 1:S
        if H[i] > 0.0
            m_ii = m_i[i]
            numerator = (g_i[i] + G[i])/m_ii - 1.0
            interaction_sum = H[i]
            for x in 1:S
                interaction_sum += M_modified[i,x]*H[x]
            end
            duH[i] = H[i]*m_ii*(numerator - interaction_sum/H_i0[i])
        else
            duH[i] = 0.0
        end
    end

    # Predator dynamics
    for α in 1:R
        if P[α] > 0.0
            predation_sum = 0.0
            for j in 1:S
                predation_sum += a_matrix[j, α]*H[j]
            end
            predator_interactions = 0.0
            for β in 1:R
                predator_interactions += A[α, β]*P[β]
            end
            duP[α] = P[α]*(epsilon[α]*predation_sum - m_alpha[α] + predator_interactions)
        else
            duP[α] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end

params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
H_init = H_i0
P_init = fill(10.0, R)
u0 = vcat(H_init, P_init)
tspan = (0.0, 500.0)
EXTINCTION_THRESHOLD = 1e-6

# Callbacks
callbacks = []
push!(callbacks, PositiveDomain())

for x in 1:S
    condition(u, t, integrator) = u[x] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[x] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

offset = S
for α in 1:R
    ind = offset + α
    condition(u, t, integrator) = u[ind] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[ind] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

cb = CallbackSet(callbacks...)

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
println("herb_pred_ratio: ", P_biomass/H_biomass)
println("NPP == ∑g_iH_i? ", isapprox(NPP, sum(g_i.*H_data[:, end]), atol=50.0))
println("∑g_iH_i = ", sum(g_i.*H_data[:, end]))
end