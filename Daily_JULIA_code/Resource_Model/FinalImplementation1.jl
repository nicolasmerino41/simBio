using Statistics, Distributions, LinearAlgebra
using DifferentialEquations, DiffEqCallbacks
using CairoMakie
#############################
# Parameters
#############################
begin
legend = false
S = 6   # Number of herbivore species
R = 3  # Number of predator species

exponent_abundance = 0.0  # exponent for SAD (power law)

NPP = 1000.0
H0_mean = NPP / S
M_mean = 1.5
H0_sd = 0.0 #H0_mean/10
M_sd = 0.0 #M_mean/10

condition_limit_number = 100.0
# For dimensionless scaling:
# We will compute bar(H) and bar(M) as the averages of H_i^0 and m_i after we draw them.

# Interaction strengths
mu = 0.5                # average herbivore-herbivore interaction strength
mu_pred = 0.01          # average predator-predator interaction strength (negative for stable)
mu_predation = 0.01      # average herbivore-predator interaction strength
asymmetry_competition = false # Very low assyemtries will lead to unstable solutions
asymmetry_predators = 0.7
asymmetry_predation = 1.0
epsilon = 0.1 # assimilation efficiency of predators

connectivity_hp = 1.0  # Herbivore-Predator connectivity
connectivity_pp = 0.4     # Predator-Predator connectivity

#############################
# 1. Generate SAD (hat_H)
#############################
function generate_hat_Hi(S::Int, exponent::Float64)
    ranks = collect(1:S)
    abundances = ranks .^ (-exponent)
    abundances ./= sum(abundances)
    return abundances
end

hat_H = generate_hat_Hi(S, exponent_abundance)

#############################
# 2. Draw H_i^0 and m_i from normal dist, ensure positivity
#############################
H_i0 = [abs(rand(Normal(H0_mean, H0_sd))) for _ in 1:S]
m_i = [abs(rand(Normal(M_mean, M_sd))) for _ in 1:S]

# Compute bar(H) and bar(M)
barH = mean(H_i0)
barM = mean(m_i)

# Compute h_i and q_i dimensionless factors:
# h_i = H_i^0/(S*barH) and q_i = m_i/(S*barM) in the original doc definition, but here we can just set:
# Actually, doc sets H_i^0 = S bar(H) h_i => h_i = H_i^0/(S*bar(H))
# and m_i = S bar(M) q_i => q_i = m_i/(S*bar(M))

h_i = [H0/(S*barH) for H0 in H_i0]
q_i = [mm/(S*barM) for mm in m_i]

# Check sums of h_i and q_i:
# doc requires sum h_i = 1 and sum q_i = 1
# If they don't sum to 1, we can normalize them:
h_sum = sum(h_i)
q_sum = sum(q_i)
h_i = h_i ./ h_sum
q_i = q_i ./ q_sum

#############################
# 3. Generate herbivore-herbivore competition matrix V
#############################
function generate_competition_matrix(S::Int, mu::Float64, symmetric::Bool; check_condition=true)
    # This function returns (V, μ_matrix) given S, mu, and whether we want a symmetric scenario or not.
    #
    # If symmetric=true:
    #   (I+μ) has diagonal = 1.0, off-diagonal = μ (perfectly uniform scenario).
    #
    # If symmetric=false:
    #   (I+μ) has diagonal = 1.0, off-diagonal random in [0,1], then scaled so mean=μ.
    #   If condition number > condition_limit_number, reduce asymmetry by blending off-diagonal values towards μ.
    #
    # After constructing (I+μ), we invert it to get V, and compute μ_matrix = (I+μ)-I.
    # Check condition number if requested.

    I_plus_mu = Matrix{Float64}(undef, S, S)

    if symmetric
        # Symmetric scenario: diagonal=1.0, off-diagonal=μ
        for i in 1:S, j in 1:S
            I_plus_mu[i,j] = (i == j) ? 1.0 : mu
        end
    else
        # Asymmetric scenario:
        # diagonal=1.0
        # off-diagonal: random values scaled to have mean=μ
        for i in 1:S, j in 1:S
            if i == j
                I_plus_mu[i,j] = 1.0
            else
                I_plus_mu[i,j] = rand()
            end
        end

        # Scale off-diagonal to have mean=μ
        off_diag_vals = [I_plus_mu[i,j] for i in 1:S, j in 1:S if i!=j]
        current_mean = mean(off_diag_vals)
        if current_mean != 0.0
            scaling_factor = mu / current_mean
            for i in 1:S, j in 1:S
                if i != j
                    I_plus_mu[i,j] *= scaling_factor
                end
            end
        else
            # If current_mean=0 (extremely unlikely), set all off-diag to μ
            for i in 1:S, j in 1:S
                if i != j
                    I_plus_mu[i,j] = mu
                end
            end
        end
    end

    if check_condition
        # Check condition number and if too large, reduce asymmetry if not symmetric
        max_attempts = 5
        attempts = 0
        while attempts < max_attempts
            cnum = cond(I_plus_mu)
            println("Condition number of (I+μ): ", cnum)
            if cnum <= condition_limit_number
                break
            end

            if symmetric
                # If symmetric and condition is large, there's not much we can do without changing parameters.
                println("Warning: Condition number still high in symmetric scenario. Consider changing μ or S.")
                break
            else
                # Reduce asymmetry by bringing off-diagonal values closer to μ
                # e.g., blend each off-diag value with μ: new_val = (val + μ)/2
                println("Condition > $condition_limit_number, reducing asymmetry by blending towards μ...")
                for i in 1:S, j in 1:S
                    if i != j
                        I_plus_mu[i,j] = (I_plus_mu[i,j] + mu)/2
                    end
                end
            end

            attempts += 1
        end
    end

    println("Final condition number of (I+μ): ", cond(I_plus_mu))
    if cond(I_plus_mu) > condition_limit_number && !symmetric
        println("Warning!!!!!!!!!: Condition number still high in asymmetric scenario. Consider changing μ or S.")
    end

    # Invert I+μ
    # Compute V = (I+μ)^{-1}
    V = inv(I_plus_mu)

    # Compute μ_matrix = (I+μ)-I
    μ_matrix = Matrix{Float64}(undef, S, S)
    for i in 1:S, j in 1:S
        if i == j
            μ_matrix[i,j] = I_plus_mu[i,j] - 1.0
        else
            μ_matrix[i,j] = I_plus_mu[i,j]
        end
    end

    return V, μ_matrix
end

V, mu_matrix = generate_competition_matrix(S, mu, asymmetry_competition)

#############################
# 4. Generate predator-predator matrix A
#    Similar approach: A_{αβ}. Predators also have self-regulation.
#    We'll assume diagonal = -1.0 for strong self-regulation, and off-diagonal mean = mu_pred
#############################
function generate_interaction_matrix_pred(R, mu_pred, asymmetry, connectivity)
    # We'll generate a baseline with diagonal = -1.0 (self-regulation)
    # and off-diagonal random with a certain mean. Only some entries present due to connectivity.
    A_matrix = fill(-1.0, R, R)
    # Fill off-diag
    # We only assign an interaction if rand() < connectivity
    # If connected, we assign a random positive/negative value scaled to achieve mean mu_pred
    # We'll do a two-step: generate random matrix, scale to mean mu_pred, then interpolate.

    # Symmetric baseline: all off-diag = mu_pred for simplicity
    A_symmetric = fill(mu_pred, R, R)
    for α in 1:R
        A_symmetric[α, α] = -1.0
    end

    A_random = Matrix{Float64}(undef, R, R)
    for α in 1:R
        for β in 1:R
            if α == β
                A_random[α, β] = -1.0
            else
                if rand() < connectivity
                    # Random positive or negative around mu_pred
                    A_random[α, β] = rand() # 0 to 1
                else
                    # No feeding link
                    A_random[α, β] = 0.0
                end
            end
        end
    end

    # Scale off-diagonal A_random to have mean off-diag = mu_pred
    off_vals = [A_random[α,β] for α in 1:R, β in 1:R if α != β]
    if !isempty(off_vals)
        current_mean = mean(off_vals)
        if current_mean != 0.0
            scaling_factor = mu_pred / current_mean
            for α in 1:R, β in 1:R
                if α != β
                    A_random[α, β] *= scaling_factor
                end
            end
        else
            # If no connections, just leave them as zero
        end
    end

    A = asymmetry * A_symmetric .+ (1.0 - asymmetry)*A_random
    return A
end

A = generate_interaction_matrix_pred(R, mu_pred, asymmetry_predators, connectivity_pp)

#############################
# 5. Herbivore-Predator feeding matrix a_{iα}
#    We'll create a SxR matrix with connectivity_hp chance of feeding link.
#    Mean interaction = mu_predation
#############################
function generate_predator_herbivore_matrix(rows::Int, cols::Int, connectivity::Float64, mu::Float64, asymmetry::Float64)
    @assert 0.0 ≤ connectivity ≤ 1.0 "connectivity must be between 0 and 1"
    @assert 0.0 ≤ asymmetry ≤ 1.0 "asymmetry must be between 0 and 1"
    @assert mu ≥ 0 "mu should be non-negative for this setup"

    mat = zeros(rows, cols)
    dist = Uniform(0, 2*mu)  # distribution with mean mu
    for i in 1:rows
        for j in 1:cols
            if rand() < connectivity
                # Realize this interaction
                # If asymmetry=1, value=mu. If 0, value~Uniform(0,2*mu).
                # Otherwise a combination:
                random_part = rand(dist) # mean mu
                value = asymmetry*mu + (1-asymmetry)*random_part
                mat[i, j] = value
            end
        end
    end
    return mat
end

a_matrix = generate_predator_herbivore_matrix(S, R, connectivity_hp, mu_predation, asymmetry_predation)

#############################
# 6. Compute C_ij and G_i from doc:
#    C_ij = Σ_{α,β} ε a_{iα} A^{-1}_{αβ} a_{jβ}
#    G_i = Σ_{α,β} a_{iα} A^{-1}_{αβ} m_β
#
# First invert A
#############################
A_inv = inv(A)

C = zeros(S,S)
G = zeros(S)

m_alpha = [M_mean for _ in 1:R] # Predator mortality rates for G_i calculation.
# Actually, doc states m_α for predators. Let's assign from A or user choice:
# We'll assume predator mortalities = M_mean for simplicity:
# If you want more complexity, draw them similarly as herbivores.

for i in 1:S
    for j in 1:S
        val = 0.0
        for α in 1:R
            for β in 1:R
                val += epsilon * a_matrix[i, α]*A_inv[α, β]*a_matrix[j, β]
            end
        end
        C[i,j] = val
    end
end

for i in 1:S
    val = 0.0
    for β in 1:R
        for α in 1:R
            val += a_matrix[i, α]*A_inv[α, β]*m_alpha[β]
        end
    end
    G[i] = val
end

#############################
# 7. Construct modified competition matrix:
#    μ_ij from V?
# The doc originally uses μ_ij and V differently. μ_ij = (1 - V_ij) maybe.
# In the doc: V is related to (I+μ)^{-1}, and μ_ij is difference from identity.
# Let's define μ_ij = 1 - V_ij (since V_ij=1 on diag, and V_ij <1 otherwise).
#
# Actually from original doc:
# If V = (I+μ)^{-1}, then μ_ij = ... But we have a chosen approach:
# Let's just derive μ_ij from V:
# For a stable system with no predators: (I+μ) = V^{-1}
# μ = V^{-1} - I
# μ_ij = (V^{-1})_ij - δ_ij
#
#############################
# WE DON'T DO THIS ANYMORE SINCE WE FIXED THE COMPETITION MATRIX BUILDER ABOVE
# V_inv = inv(V)
# mu_matrix = V_inv - I
# for i in axes(mu_matrix, 1), j in axes(mu_matrix, 2)
#     if i == j
#         mu_matrix[i, j] = 0.0
#     end
# end
#############################
# 8. Modified interaction due to predators:
#    μ_ij -> μ_ij + (C_ij * H_i^0)/m_i
#############################
M_modified = copy(mu_matrix)
for i in 1:S, j in 1:S
    M_modified[i,j] = mu_matrix[i,j] + (C[i,j]*H_i0[i]/m_i[i])
end

#############################
# 9. Compute hat_p_i:
#    hat_p_i = [ hat_H_i + Σ_j M_modified[i,j]*hat_H_j ] / h_i[i]
#############################
hat_p = similar(hat_H)
for i in 1:S
    interaction_sum = 0.0
    for j in 1:S
        interaction_sum += M_modified[i,j]*hat_H[j]
    end
    hat_p[i] = (hat_H[i] + interaction_sum)/h_i[i]
end

p = hat_p ./ sum(hat_p)

#############################
# 10. Solve for x using NPP condition:
# NPP = Σ_i (g_i + G_i)H_i
#
# g_i = x p_i m_i (from doc: g_i/m_i = x p_i => g_i = x p_i m_i)
#
# Also, we must find H_i at equilibrium:
# From doc with predators:
# (1/H_i) dH_i/dt=0 => (g_i+G_i)/m_i -1 = (H_i + Σ_j M_modified[i,j] H_j)/H_i^0
#
# In matrix form:
# Let M = M_modified for simplicity.
# Define vector R_i = (g_i+G_i)/m_i -1 = x p_i + G_i/m_i -1
# Then (I+M)(H/H^0) = R
# H/H^0 = (I+M)^(-1) R
#
# R_i = x p_i + G_i/m_i -1
# Then H_i = H_i^0 * Σ_j (I+M)^(-1)[i,j] (x p_j + G_j/m_j -1)
#
# Then NPP = Σ_i (x p_i m_i + G_i)*H_i must be solved for x.
#############################
IplusM = I + M_modified
IM_inv = inv(IplusM)

R_vector = [ (x -> x*p[j] + G[j]/m_i[j] -1 ) for j=1:S ] # function form

# We get a quadratic in x, let's precompute terms:
# H_i(x) = H_i^0 [ x A_i + B_i ], where
# A_i = Σ_j IM_inv[i,j]*p[j]
# B_i = Σ_j IM_inv[i,j]*(G[j]/m_i[j]-1)

A_vec = zeros(S)
B_vec = zeros(S)
for i in 1:S
    # Compute A_i and B_i
    A_val = 0.0
    B_val = 0.0
    for j in 1:S
        A_val += IM_inv[i,j]*p[j]
        B_val += IM_inv[i,j]*(G[j]/m_i[j] -1)
    end
    A_vec[i] = A_val
    B_vec[i] = B_val
end

# NPP(x) = Σ_i (x p_i m_i + G_i)*H_i
#        = Σ_i (x p_i m_i + G_i)*H_i^0 ( x A_i + B_i )
# Expand:
# = Σ_i [ x^2 p_i m_i H_i^0 A_i + x p_i m_i H_i^0 B_i + x G_i H_i^0 A_i + G_i H_i^0 B_i ]

# Collect terms in powers of x:
# x^2: Σ_i p_i m_i H_i^0 A_i
# x^1: Σ_i [p_i m_i H_i^0 B_i + G_i H_i^0 A_i]
# x^0: Σ_i G_i H_i^0 B_i

A_coef = sum(p[i]*m_i[i]*H_i0[i]*A_vec[i] for i in 1:S)
B_coef = sum(p[i]*m_i[i]*H_i0[i]*B_vec[i] + G[i]*H_i0[i]*A_vec[i] for i in 1:S)
C_coef = sum(G[i]*H_i0[i]*B_vec[i] for i in 1:S) - NPP

# Solve quadratic A_coef x^2 + B_coef x + C_coef = 0
# for x, choose positive root
discriminant = B_coef^2 - 4*A_coef*C_coef
if discriminant < 0
    error("No real solution for x found!")
end

x_candidates = [(-B_coef + sqrt(discriminant))/(2*A_coef), (-B_coef - sqrt(discriminant))/(2*A_coef)]
x_sol = maximum(filter(x->x>0, x_candidates))
if isempty(x_sol)
    error("No positive solution for x!")
end
x_final = x_sol

#############################
# 11. Compute g_i:
# g_i = x p_i m_i
#############################
g_i = [x_final * p[i] * m_i[i] for i in 1:S]

#############################
# Print results
#############################
# println("p: ", p)
# println("Sum p: ", sum(p))
# println("x: ", x_final)
# println("g_i: ", g_i)
# println("G_i: ", G)
println("g_i/m_i: ", g_i./m_i)
println("g_i/m_i: ", (g_i.+G)./m_i)
# println("C_ij: ", C)

# Initial conditions (example):
H_init = [H_i0[i] for i in 1:S] # start at characteristic density
P_init = fill(10.0, R)          # initial predator densities
u0 = vcat(H_init, P_init)

# Time span
tspan = (0.0, 500.0)

# Extinction threshold
EXTINCTION_THRESHOLD = 1e-6

#############################
# Define the ODE system
#############################
function ecosystem_dynamics!(du, u, p, t)
    # p holds all parameters, we destructure them:
    S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = p

    H = @view u[1:S]
    P = @view u[S+1:S+R]

    duH = zeros(S)
    duP = zeros(R)

    # Herbivore dynamics
    # dH_i/dt = H_i m_i [ ((g_i+G_i)/m_i -1 ) - (H_i + Σ_j M_ij H_j)/H_i^0 ]
    for i in 1:S
        if H[i] > 0.0
            m_ii = m_i[i]
            numerator = (g_i[i] + G[i])/m_ii - 1.0
            interaction_sum = H[i]
            for j in 1:S
                interaction_sum += M_modified[i,j]*H[j]
            end
            duH[i] = H[i]*m_ii*(numerator - interaction_sum/H_i0[i])
        else
            duH[i] = 0.0
        end
    end

    # Predator dynamics
    # dP_α/dt = P_α [ ε Σ_j a_{jα} H_j - m_α + Σ_β A_{αβ} P_β ]
    # Compute predation terms
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
            duP[α] = P[α]*(epsilon*predation_sum - m_alpha[α] + predator_interactions)
        else
            duP[α] = 0.0
        end
    end

    du[1:S] = duH
    du[S+1:S+R] = duP
end

#############################
# Set up ODE problem
#############################
params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)

#############################
# Define callbacks for extinction and positivity
#############################
callbacks = []
push!(callbacks, PositiveDomain())

# Herbivores extinction callback
for i in 1:S
    condition(u, t, integrator) = u[i] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[i] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

# Predators extinction callback
offset = S
for α in 1:R
    idx = offset + α
    condition(u, t, integrator) = u[idx] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[idx] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

cb = CallbackSet(callbacks...)

#############################
# Solve the system
#############################
sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

#############################
# Extract and plot results
#############################
times = sol.t
H_data = sol[1:S, :]
P_data = sol[S+1:S+R, :]

fig = Figure()
ax = Axis(fig[1,1], xlabel="Time", ylabel="Density", title="Dynamics")
for i in 1:S
    lines!(ax, times, H_data[i, :], label="H$i")
end
for α in 1:R
    lines!(ax, times, P_data[α, :], label="P$α", linestyle=:dash)
end
if legend 
    axislegend(ax, position=:rt)
end
display(fig)

#############################
# Summary statistics
#############################
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