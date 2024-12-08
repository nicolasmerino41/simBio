using Statistics, Distributions, LinearAlgebra

#############################
# Parameters
#############################
S = 5   # Number of herbivore species
R = 3   # Number of predator species

exponent_abundance = 0.5  # exponent for SAD (power law)
H0_mean = 100.0
M_mean = 0.1
H0_sd = H0_mean/10
M_sd = M_mean/10
NPP = 1000.0

# For dimensionless scaling:
# We will compute bar(H) and bar(M) as the averages of H_i^0 and m_i after we draw them.

# Interaction strengths
mu = 0.1                 # average herbivore-herbivore interaction strength
mu_pred = -0.05           # average predator-predator interaction strength (negative for stable)
mu_predation = 0.001      # average herbivore-predator interaction strength
asymmetry_competition = 0.8
asymmetry_predators = 0.7
asymmetry_predation = 0.5
epsilon = 0.1             # assimilation efficiency of predators

connectivity_hp = 0.7     # Herbivore-Predator connectivity
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
function generate_competition_matrix(S, mu, asymmetry)
    V_symmetric = ones(S, S)
    for i in 1:S, j in 1:S
        if i == j
            V_symmetric[i, j] = 1.0
        else
            V_symmetric[i, j] = 1/(1+mu)
        end
    end

    V_random = ones(S, S)
    for i in 1:S, j in 1:S
        if i != j
            V_random[i, j] = rand()
        else
            V_random[i, j] = 1.0
        end
    end

    off_diag_values = [V_random[i,j] for i in 1:S, j in 1:S if i != j]
    current_mean = mean(off_diag_values)
    desired_mean = 1/(1+mu)
    scaling_factor = desired_mean / current_mean

    for i in 1:S, j in 1:S
        if i != j
            V_random[i,j] *= scaling_factor
        end
    end

    V = asymmetry * V_symmetric + (1.0 - asymmetry)*V_random
    return V
end

V = generate_competition_matrix(S, mu, asymmetry_competition)

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
function generate_predation_matrix(S, R, mu_predation, asymmetry, connectivity)
    # Similar approach: we want diagonal doesn't matter here, just off-diag meaning i-th herb feeds predator α
    # Actually, there's no diagonal concept here since it's SxR.
    # We'll create a baseline with mean mu_predation and random variation.

    # Baseline: everything = mu_predation
    A_symmetric = fill(mu_predation, S, R)

    A_random = Matrix{Float64}(undef, S, R)
    for i in 1:S, α in 1:R
        if rand() < connectivity
            A_random[i, α] = rand() # random 0 to 1
        else
            A_random[i, α] = 0.0
        end
    end

    # Scale A_random off elements to have mean = mu_predation
    vals = [A_random[i, α] for i in 1:S, α in 1:R if A_random[i, α] != 0.0]
    if !isempty(vals)
        current_mean = mean(vals)
        if current_mean != 0.0
            scaling_factor = mu_predation / current_mean
            for i in 1:S, α in 1:R
                A_random[i, α] *= scaling_factor
            end
        end
    end

    a_matrix = asymmetry * A_symmetric .+ (1.0 - asymmetry)*A_random
    return a_matrix
end

a_matrix = generate_predation_matrix(S, R, mu_predation, asymmetry_predation, connectivity_hp)

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
V_inv = inv(V)
mu_matrix = V_inv .- I

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
println("p: ", p)
println("Sum p: ", sum(p))
println("x: ", x_final)
println("g_i: ", g_i)
println("G_i: ", G)
println("C_ij: ", C)
