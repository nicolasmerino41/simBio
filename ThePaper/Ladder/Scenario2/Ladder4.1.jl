using Pkg
cd(pwd())

dir = pwd()
# Packages
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Dates, Distributions, Serialization, StatsBase, Random
using DifferentialEquations, DiffEqCallbacks, LinearAlgebra, Logging, ForwardDiff
using GLM, Graphs

const DF = DataFrames
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
const EXTINCTION_THRESHOLD = 1e-6
################ INFORMATION #######################
# a) Ladder1.jl is for all the functions
# b) Ladder2.jl is for running the simulations
# c) Ladder3.jl is for post-processing and plotting
#############################################################################
#############################################################################
################## TROPHIC DYNAMICS ODE MODEL ###############################
#############################################################################
#############################################################################
function trophic_ode!(du, u, p, t)
    # Unpack parameters from tuple
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    
    # u is arranged as:
    # u[1:R]          ? resources
    # u[R+1 : R+C]  ? consumers
    
    # Resources dynamics (indices 1:R)
    for i in 1:R
        # Sum over predators of resource i (these are consumers, indices: R+1 to R+C)
        pred_sum = 0.0
        for j in R+1:R+C
            pred_sum += A[i, j] * u[j]
        end
        # Resource i dynamics
        du[i] = u[i] * d_res[i] * ( (r_res[i] / d_res[i]) - u[i] + pred_sum )
    end

    # Consumers dynamics (indices R+1 to R+C)
    for k in 1:C
        i = R + k  # global index for consumer k
        
        # Gains from prey (only A_used[i,j]>0)
        # println("A[i,j] = ", typeof(A))
        # println("epsilon[i,j] = ", typeof(epsilon))
        sum_prey = sum(
            epsilon[i,j] * A[i,j] * u[j]
            for j in 1:R+C if A[i,j] > 0.0;
            init = 0.0
        )

        # Losses to predators (only A_used[i,j]<0)
        sum_pred = sum(
            A[i,j] * u[j]
            for j in 1:R+C if A[i,j] < 0.0;
            init = 0.0
        )
        
        du[i] = u[i] * (m_cons[k] / xi_cons[k]) * ( -xi_cons[k] - u[i] + sum_prey + sum_pred )
    end
end

#############################################################################
#############################################################################
################## CALIBRATE PARAMETRISATION ################################
#############################################################################
#############################################################################
# ────────────────────────────────────────────────────────────────────────────────
# 1) Compute “safe” ξ and K so that M is strictly diagonally dominant
# ────────────────────────────────────────────────────────────────────────────────
function compute_default_thresholds(A::AbstractMatrix, ε::AbstractMatrix, R::Int; margin::Float64=0.05)
    S = size(A,1)
    C = S - R

    # (a) consumer thresholds ξ_k
    xi_cons = zeros(Float64, C)
    for k in 1:C
        i = R + k
        gain = sum(ε[i,j]*A[i,j] for j in 1:R if A[i,j] > 0.0; init=0.0)
        loss = sum(abs(A[j,i])    for j in 1:R if A[j,i] < 0.0; init=0.0)
        xi_cons[k] = (gain + loss)*(1 + margin)
    end

    # (b) resource capacities K_i
    K_res = zeros(Float64, R)
    for i in 1:R
        drain = sum(A[R+j,i] for j in 1:C if A[R+j,i] > 0.0; init=0.0)
        K_res[i] = drain*(1 + margin)
    end

    return xi_cons, K_res
end

# ────────────────────────────────────────────────────────────────────────────────
# 2) Given ξ and K build & solve M·B_eq = b, check positivity
# ────────────────────────────────────────────────────────────────────────────────
function calibrate_from_K_xi(
    xi_cons::AbstractVector,    # length C
    K_res::AbstractVector,      # length R
    ε::AbstractMatrix,          # S×S
    A::AbstractMatrix           # S×S
)
    R = length(K_res)
    S = size(A,1)
    C = S - R

    M = zeros(Float64, S, S)
    b = zeros(Float64, S)

    # ─── Resources (i=1:R):
    #    B_i - ∑_j A[i,j] * B_j =  K_res[i]
    for i in 1:R
        M[i,i] = 1.0
        for j in 1:S
            # move ∑ A[i,j]B_j to LHS
            M[i,j] -= A[i,j]
        end
        b[i] = K_res[i]
    end

    # ─── Consumers (i=R+1:S), indexed k=i−R:
    #    B_i 
    #   - ∑_{A[i,j]>0} ε[i,j] A[i,j] B_j 
    #   - ∑_{A[i,j]<0} (−A[i,j]) B_j 
    #   = -ξ_cons[k]
    for k in 1:C
        i = R + k
        M[i,i] = 1.0
        for j in 1:S
            if A[i,j] > 0
                # feeding gain term
                M[i,j] -=  ε[i,j] * A[i,j]
            elseif A[i,j] < 0
                # predation loss term (A[i,j]<0 => -A[i,j]>0)
                M[i,j] -= -A[i,j]
            end
        end
        # note the minus on the RHS
        b[i] = -xi_cons[k]
    end

    B_eq = M \ b
    return any(B_eq .<= 0.0) ? nothing :
           (B_eq[1:R], B_eq[R+1:end])
end

# ────────────────────────────────────────────────────────────────────────────────
# 3) Wrapper replacing your old cal_param
# ────────────────────────────────────────────────────────────────────────────────
"""
    cal_param(R::Int, C::Int, ε, A; margin)

Auto-computes ξ_cons, K_res, then solves for (R_eq, C_eq).  
Returns `(R_eq, C_eq, xi_cons, K_res)` or `nothing` if infeasible.
"""
function cal_param(R::Int, C::Int, ε::AbstractMatrix, A::AbstractMatrix; margin::Float64=0.05)
    xi_cons, K_res = compute_default_thresholds(A, ε, R; margin=margin)
    eq = calibrate_from_K_xi(xi_cons, K_res, ε, A)
    if eq === nothing
        return nothing
    else
        R_eq, C_eq = eq
        return (R_eq, C_eq, xi_cons, K_res)
    end
end

"""
    calibrate_from_K_xi(
      K_res::Vector{<:Real},    # length R: resource carrying capacities
      xi_cons::Vector{<:Real},  # length C: consumer thresholds
      ε::AbstractMatrix,        # S×S conversion efficiencies
      A::AbstractMatrix         # S×S interaction signs/magnitudes
    ) -> (R_eq, C_eq) or nothing

Solve the steady‐state abundances B = [R_eq; C_eq] from:

  For i=1:R (resources):
    0 = −K_i − B_i  − ∑_{j=R+1}^S A_{j,i} * B_j

  For i=R+1:S (consumers, indexed k=i−R):
    0 = −ξ_k − B_i + ∑_{j=1}^R ε_{i,j}*A_{i,j}*B_j − ∑_{j=1}^R A_{j,i}*B_j

Returns nothing if no positive solution exists.
"""
# ────────────────────────────────────────────────────────────────────────────────
# 2) Given ξ and K build & solve M·B_eq = b, check positivity
# ────────────────────────────────────────────────────────────────────────────────
function calibrate_from_K_xi(xi_cons::AbstractVector, K_res::AbstractVector,
                             ε::AbstractMatrix, A::AbstractMatrix)
    R = length(K_res)
    C = length(xi_cons)
    S = R + C

    M = zeros(S,S)
    b = zeros(S)

    # — resource rows: B_i + sum_j A[j,i]*B_j = K_res[i]
    for i in 1:R
        M[i,i] = 1.0
        for j in R+1:S
            M[i,j] = A[j,i]
        end
        b[i] =  K_res[i]    # <--- NO minus
    end

    # — consumer rows: B_i - sum εA_{i,j}B_j + sum |A_{j,i}|B_j = xi_cons[k]
    for k in 1:C
        i = R + k
        M[i,i] = 1.0
        for j in 1:R
            if A[i,j] > 0.0
                M[i,j] -= ε[i,j]*A[i,j]
            end
        end
        for j in 1:S
            if A[j,i] < 0.0
                M[i,j] += -A[j,i]
            end
        end
        b[i] = xi_cons[k]    # correct sign here
    end

    B_eq = M \ b
    return any(B_eq .<= 0.0) ? nothing : (B_eq[1:R], B_eq[R+1:end])
end

#############################################################################
#############################################################################
############################# make_A ########################################
#############################################################################
#############################################################################
function make_A(A::AbstractMatrix, R::Int, conn::Float64, scenario::Symbol;
                pareto_exponent::Float64=1.75, mod_gamma::Float64=5.0)
    S = size(A,1)
    C = S - R
    fill!(A, 0.0)

    prey_indices(i) = 1:R

    if scenario == :ER
        for i in (R+1):S, j in 1:R
            if rand() < conn
                A[i,j] = abs(randn()); A[j,i] = -abs(randn())
            end
        end

    elseif scenario == :PL
        raw = rand(Pareto(1.0, pareto_exponent), C)
        ks  = clamp.(floor.(Int,raw), 1, R)
        for (idx,k) in enumerate(ks)
            i = R+idx
            for j in sample(1:R, min(k,R); replace=false)
                A[i,j] = abs(randn()); A[j,i] = -abs(randn())
            end
        end

    elseif scenario == :MOD
        halfR,halfC = fld(R,2),fld(C,2)
        res1, res2 = 1:halfR, (halfR+1):R
        con1, con2 = (R+1):(R+halfC), (R+halfC+1):S

        for i in (R+1):S, j in 1:R
            same = (i in con1 && j in res1) || (i in con2 && j in res2)
            p = same ? conn*mod_gamma : conn/mod_gamma
            if rand() < clamp(p,0,1)
                A[i,j] = abs(randn()); A[j,i] = -abs(randn())
            end
        end
    else
        error("Unknown scenario $scenario")
    end

    return A
end

#############################################################################
#############################################################################
################## JACOBIAN, RESILIENCE AND RESISTANCE ######################
#############################################################################
#############################################################################
function compute_jacobian(B, p)
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    S = R + C

    # 1) Build D
    D = zeros(S,S)
    for i in 1:R
        D[i,i] = d_res[i] * B[i]
    end
    for k in 1:C
        i = R + k
        alpha = m_cons[k] / xi_cons[k]
        D[i,i] = alpha * B[i]
    end

    # 2) Build Mstar = -I + A*
    Mstar = Matrix{Float64}(-I(S))
    # resources → anything
    for i in 1:R, j in 1:S
        if i != j && A[i,j] != 0.0
            Mstar[i,j] += A[i,j]
        end
    end
    # consumers → anything
    for i in R+1:S
        for j in 1:S
            if A[i,j] > 0.0
                Mstar[i,j] += epsilon[i,j] * A[i,j]
            elseif A[j,i] < 0.0
                Mstar[i,j] += A[i,j]
            end
        end
    end

    return D, Mstar
end

# Resilience: negative of the largest real part of the Jacobian eigenvalues.
function compute_resilience(B, p)
    D, Mstar = compute_jacobian(B, p)
    J = D * Mstar
    ev = eigvals(J)
    return maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of the symmetric part of the Jacobian.
function compute_reactivity(B, p)
    D, Mstar = compute_jacobian(B, p)
    J = D * Mstar
    J_sym = (J + J') / 2
    ev_sym = eigvals(J_sym)
    return maximum(real.(ev_sym))
end

#############################################################################
#############################################################################
################## SIMULATE PRESS PERTURBATIONS #############################
#############################################################################
#############################################################################
function simulate_press_perturbation(
    u0, p, tspan, t_perturb, delta;
    solver=Tsit5(),
    plot=false,
    show_warnings=true,
    full_or_simple=true,
    cb = cb,
    species_specific_perturbation=false
)
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p

    # --- Phase 1: run up to t_perturb ---
    prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], t_perturb), p)
    sol1 = show_warnings ? solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8) :
                           with_logger(logger) do
                               solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8)
                           end
    pre_state = sol1.u[end]
    before_persistence = count(x -> x > EXTINCTION_THRESHOLD, pre_state) / length(pre_state)
    # pre_state[pre_state .> EXTINCTION_THRESHOLD] .= 0.0

    ########## FIRST PART: WHOLE COMMUNITY PERTURBATION ##########
    # --- Phase 2: apply press (reduce thresholds by delta) ---
    xi_press = xi_cons .+ delta
    # r_press = r_res .- delta
    p_press  = (R, C, m_cons, xi_press, r_press, d_res, epsilon, A)
    prob2    = ODEProblem(trophic_ode!, pre_state, (t_perturb, tspan[2]), p_press)
    sol2     = solve(prob2, solver; callback = cb, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)
    after_persistence = count(x -> x > EXTINCTION_THRESHOLD, new_equil) / length(new_equil)

    # --- Return Times ---
    return_times = fill(NaN, n)
    for i in 1:n
        target = new_equil[i]
        for (t, state) in zip(sol2.t, sol2.u)
            if abs(state[i] - target) / (abs(target) + 1e-8) < 0.10
                return_times[i] = t - t_perturb
                break
            end
        end
    end

    # --- Maximum Overshoot & Integrated Recovery Error ---
    overshoot = zeros(n)
    ire       = zeros(n)
    for i in 1:n
        # compute relative errors at each timepoint
        rel_errors = [abs(state[i] - new_equil[i]) / (abs(new_equil[i]) + 1e-8)
                      for state in sol2.u]
        overshoot[i] = maximum(rel_errors)
        ire[i]       = mean(rel_errors)
    end

    ########## SECOND PART: SPECIES SENSITIVITY TO PERTURBATION ##########
    species_specific_return_times = fill(NaN, C)

    if species_specific_perturbation
        for i in 1:C
            xi_press = copy(xi_cons)          # <-- copy, not alias
            xi_press[i] *= (1 - delta)         # perturb only species i
    
            p_press  = (R, C, m_cons, xi_press, r_res, d_res, epsilon, A)
            prob2    = ODEProblem(trophic_ode!, pre_state, (t_perturb, tspan[2]), p_press)
            sol2     = solve(prob2, solver; callback=cb, reltol=1e-8, abstol=1e-8)
    
            new_equil = sol2.u[end]
    
            # --- Return Times ---
            target = new_equil[R+i]
            for (t, state) in zip(sol2.t, sol2.u)
                if abs(state[i] - target) / (abs(target) + 1e-8) < 0.10
                    species_specific_return_times[i] = t - t_perturb
                    break
                end
            end
        end
    end

    # --- Optional plotting ---
    if plot && sol1.t[end] == t_perturb && all(isfinite, pre_state)
        fig = Figure(; size=(1600,800))
        ax1 = Axis(fig[1,1];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "Before Press (Full Model)" :
                                           "Before Press (Simplified)")
        ax2 = Axis(fig[1,2];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "After Press (Full Model)" :
                                           "After Press (Simplified)")
        for i in 1:R
            lines!(ax1, sol1.t[1:end-1], sol1[i, 1:end-1], color=:blue)
            lines!(ax2, sol2.t, sol2[i, :], color=:blue)
        end
        for i in R+1:R+C
            lines!(ax1, sol1.t[1:end-1], sol1[i, 1:end-1], color=:red)
            lines!(ax2, sol2.t, sol2[i, :], color=:red)
        end
        lines!(ax1, sol1.t[1:end-1], fill(0.0, length(sol1[1, 1:end-1])), color = :black)
        lines!(ax2, sol2.t, fill(0.0, length(sol2[1, :])), color = :black)
        display(fig)
    end

    return return_times, overshoot, ire, before_persistence, after_persistence, new_equil, species_specific_return_times
end

# -------------------------------
# Example parameter definitions
if false
    # Define target equilibrium abundances
    R_eq = fill(5.0, R)  # Equilibrium biomass for each resource
    C_eq = fill(1.0, C)  # Equilibrium biomass for each consumer

    # Compute the new parameters
    new_xi_cons, new_r_res = calibrate_parametrization(R_eq, C_eq, p_tuple)

    println("New consumer thresholds (xi_cons): ", new_xi_cons)
    println("New resource growth rates (r_res): ", new_r_res)

    # Numbers of resources and consumers
    R = 2      # e.g., two resources
    C = 3      # e.g., three consumers

    # Consumer parameters (vectors of length C)
    m_cons  = ones(C)             # Mortality rates (example: all ones)
    xi_cons = ones(C)             # Overpopulation thresholds

    # Resource parameters (vectors of length R)
    r_res = ones(R)               # Intrinsic growth rates
    d_res = ones(R)               # Scaling factors for resources

    # Conversion efficiency
    epsilon = 0.5

    # Global interaction matrix A (size: (R+C) x (R+C))
    # Here we create a random matrix for demonstration.
    A = randn(R+C, R+C)

    # Bundle parameters into a tuple in the specified order
    p_tuple = (R, C, m_cons, new_xi_cons, new_r_res, d_res, epsilon, A)
    # -------------------------------
    # Initial conditions:
    # Resources start at 5.0, consumers at 1.0
    u0 = [ones(R) .* 5.0; ones(C) .* 1.0]
    tspan = (0.0, 50.0)

    # Set up and solve the ODE problem
    prob = ODEProblem(trophic_ode!, u0, tspan, p_tuple)
    sol = solve(prob, Tsit5())

    # Plot the solution
    # Labels: first R are resources, next C are consumers.
    begin
        fig = Figure(; size=(800, 600), title="Trophic Dynamics")
        ax = Axis(fig[1, 1], xlabel="Time", ylabel="Population Size", title="Trophic Dynamics")
        for i in 1:(R)
            MK.lines!(ax, sol.t, sol[i, :], label="Resource $i", color=:blue, linewidth=1)
        end
        for i in 1:(C)
            MK.lines!(ax, sol.t, sol[R+i, :], label="Consumer $i", color=:red, linewidth=1)
        end
        display(fig)
    end
end
#############################################################################
#############################################################################
################## Build Jacobian & Analytical V ############################
#############################################################################
#############################################################################
function build_jacobian(B_eq, p)
    R, C, m_cons, xi_cons, r_res, d_res, epsilon, A = p
    total = R + C
    J = zeros(total, total)

    # Resource block
    for i in 1:R
        # diff(du_i)/diffB_i = - d_res[i] * B_eq[i]
        J[i,i] = -d_res[i]*B_eq[i]
        # diff(du_i)/diffB_j (j a consumer)
        for j in (R+1):total
            J[i,j] =  d_res[i]*B_eq[i]*A[i,j]
        end
    end

    # Consumer block
    for k in 1:C
        i = R + k
        psi = B_eq[i] / xi_cons[k]
        prefactor = m_cons[k] * psi
        for j in 1:total
            # -delta_{ij} + epsilon[i,j]*A[i,j]  - A[j,i]
            J[i,j] = prefactor * ( (i==j ? -1 : 0) + epsilon[i,j]*A[i,j] - A[j,i] )
        end
    end

    return J
end

"""
    compute_analytical_V(J, R, C, m_cons, xi_cons)

Return V = J^{-1} * D, with D_{ii}=m_i/xi_i for consumers (i>R), D_{ii}=0 for resources.
"""
function compute_analytical_V(J, R, C, m_cons, xi_cons)
    total = size(J,1)
    D = zeros(total, total)
    for k in 1:C
        D[R+k, R+k] = m_cons[k] / xi_cons[k]
    end
    # Solve J * V = D  ?  V = J \ D
    return J \ D
end

##########################################################################
##########################################################################
################### Ladder-of-simplification transforms ##################
##########################################################################
##########################################################################
function transform_for_ladder_step(step, A_adj, epsilon_full)
    total = size(A_adj,1)

    if step == 1 # Full
        return A_adj, epsilon_full
    elseif step == 2 # S2 epsilon mean per row
        epsilon2 = similar(epsilon_full)
        for i in 1:total
            epsilon2[i,:] .= mean(epsilon_full[i,:])
        end
        return A_adj, epsilon2
    elseif step == 3 # S3 epsilon global mean
        return A_adj, fill(mean(epsilon_full), total, total)
    elseif step == 4 # S4 epsilon re-randomised
        epsilon4 = similar(epsilon_full)
        for i in 1:total, j in 1:total
            epsilon4[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
        end
        return A_adj, epsilon4
    elseif 5 <= step <= 8 
        
        pos, neg = mean(A_adj[A_adj.>0]), mean(A_adj[A_adj.<0])
        
        A5 = map(x-> x>0 ? pos : x<0 ? neg : 0, A_adj)
        if step==5 # epsilon is full
            return A5, epsilon_full
        elseif step==6 # epsilon is mean per row
            epsilon6 = similar(epsilon_full)
            for i in 1:total 
                epsilon6[i,:] .= mean(epsilon_full[i,:])
            end
            return A5, epsilon6
        elseif step==7 # epsilon is global mean
            epsilon7 = fill(mean(epsilon_full), total, total)
            return A5, epsilon7
        elseif step==8 # epsilon is re-randomised
            epsilon8 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon8[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A5, epsilon8
        end
    elseif 9 <= step <= 12 # A is averaged across all non-zero values
        
        m = mean(abs.(A_adj[A_adj .!= 0]))
        A6 = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)
        
        if step==9 # epsilon is full
            return A6, epsilon_full
        elseif step==10 # epsilon is mean per row
            epsilon10 = similar(epsilon_full)
            for i in 1:total
                epsilon10[i,:] .= mean(epsilon_full[i,:])
            end
            return A6, epsilon10
        elseif step==11 # epsilon is global mean
            return A6, fill(mean(epsilon_full), total, total)
        elseif step==12 # epsilon is re-randomised
            epsilon12 = similar(epsilon_full)
            for i in 1:total*total
                epsilon12[i] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A6, epsilon12
        end
    elseif 13 <= step <= 16 # A is fully randomised with same sparsity pattern as A_adj
        # Fully randomized A matrix using Normal(mean, std) of abs non-zero A entries
        A_vals = abs.(A_adj[A_adj .!= 0])
        A_mean, A_std = mean(A_vals), std(A_vals)
    
        A_rand = similar(A_adj)
        for i in 1:total, j in 1:total
            if A_adj[i, j] != 0
                A_rand[i, j] = rand(Normal(A_mean, A_std)) * sign(A_adj[i, j])
            else
                A_rand[i, j] = 0.0
            end
        end
    
        if step == 13 # epsilon is full
            return A_rand, epsilon_full
        elseif step == 14 # epsilon is mean per row
            epsilon14 = similar(epsilon_full)
            for i in 1:total
                epsilon14[i, :] .= mean(epsilon_full[i, :])
            end
            return A_rand, epsilon14
        elseif step == 15 # epsilon is global mean
            return A_rand, fill(mean(epsilon_full), total, total)
        elseif step == 16 # epsilon is re-randomised
            epsilon16 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon16[i, j] = clamp(rand(Normal(mean(epsilon_full), std(epsilon_full))), 0, 1)
            end
            return A_rand, epsilon16
        end    
    else
        error("Unknown ladder step $step")
    end
end

#########################################################################
#########################################################################
######################### BUILD CALLBACKS ###############################
#########################################################################
#########################################################################
function build_callbacks(S::Int, EXTINCTION_THRESHOLD::Float64)
    # A) Continuous threshold-based extinctions
    callbacks = []

    # 1) Always ensure positivity
    push!(callbacks, PositiveDomain())

    # 2) Herbivores: set to zero if below EXTINCTION_THRESHOLD
    for x in 1:S
        function threshold_condition(u, t, integrator)
            # event if u[x] < EXTINCTION_THRESHOLD
            return u[x] - EXTINCTION_THRESHOLD
        end
        function threshold_affect!(integrator)
            integrator.u[x] = 0.0
        end
        push!(callbacks, ContinuousCallback(threshold_condition, threshold_affect!))
    end

    # Build callback set => no forced extinction
    cb_no_trigger = CallbackSet(callbacks...)

    return cb_no_trigger
end