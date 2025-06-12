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
    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p
    S = R + C

    ### Resources (1:R) ###
    for i in 1:R
        # predation term ∑ A[i,j] * u[j]  (j = consumers)
        pred_sum = zero(eltype(u))
        for j in R+1:S
            pred_sum += A[i,j] * u[j]
        end

        # logistic + predation:
        # du[i] = r_i * u[i] * (1 - u[i]/K_i + pred_sum/K_i)
        du[i] = r_res[i]/K_res[i] * u[i] * (K_res[i] - u[i] + pred_sum)
    end

    ### Consumers (R+1 : R+C) ###
    for k in 1:C
        i = R + k

        sum_prey = zero(eltype(u))
        sum_pred = zero(eltype(u))
        for j in 1:S
            if A[i,j] > 0
                sum_prey += epsilon[i,j] * A[i,j] * u[j]
            elseif A[i,j] < 0
                sum_pred += A[i,j] * u[j]
            end
        end

        # du[i] = u[i]*(m_i/xi_i)*(-xi_i - u[i] + gains + losses)
        du[i] = u[i] * (m_cons[k] / xi_cons[k]) * (-xi_cons[k] - u[i] + sum_prey + sum_pred)
    end
end

#############################################################################
#############################################################################
################## CALIBRATE PARAMETRISATION ################################
#############################################################################
#############################################################################
# ────────────────────────────────────────────────────────────────────────────────
# 1) Compute safe Xi and K so that M is strictly diagonally dominant
# ────────────────────────────────────────────────────────────────────────────────
# function compute_default_thresholds(
#     A::AbstractMatrix,
#     epsilon::AbstractMatrix,
#     R::Int;
#     margin::Float64 = 1.0
# )
#     S = size(A,1)
#     C = S - R

#     # 1) compute the raw gains+losses per consumer
#     gl = similar(zeros(C))
#     for k in 1:C
#         i = R + k
#         gain = sum(epsilon[i,j]*A[i,j] for j in 1:R if A[i,j] > 0; init =0.0)
#         loss = sum(-A[i,j] for j in R+1:S if A[i,j] < 0; init =0.0)
#         gl[k] = gain + loss
#     end

#     # 2) community‐wide floor for xi:
#     xi_floor = minimum(gl) * (1 + margin)

#     # 3) now randomise every xiᵢ uniformly *above* that floor:
#     xi_cons = xi_floor .+ rand(C) .* xi_floor

#     # 4) same for resources: raw drains
#     dr = similar(zeros(R))
#     for i in 1:R
#         dr[i] = sum(-A[i,j] for j in R+1:S if A[i,j] < 0; init=0.0)
#     end

#     # 5) community‐wide floor for K:
#     K_floor = minimum(dr) * (1 + margin)

#     # 6) randomise each Kᵢ above that floor:
#     K_res = K_floor .+ rand(R) .* K_floor

#     return xi_cons, K_res
# end

# function compute_default_thresholds(A::AbstractMatrix, epsilon::AbstractMatrix, R::Int; margin::Float64=0.1)
#     S = size(A,1)
#     C = S - R

#     xi_cons = zeros(Float64, C)
#     for k in 1:C
#         i = R + k
#         gain = sum(epsilon[i,j]*A[i,j] for j in 1:R if A[i,j] > 0; init =0.0)
#         loss = sum(-A[i,j] for j in R+1:S if A[i,j] < 0; init =0.0)
#         xi_cons[k] = (gain + loss)*(1 + abs(rand(Normal(3.0, 1.0))))
#     end

#     K_res = zeros(Float64, R)
#     for i in 1:R
#         drain = sum(-A[i,j] for j in R+1:S if A[i,j] < 0; init = 0.0)
#         K_res[i] = drain*(1 + abs(rand(Normal(30.0, 3.0))))
#     end

#     return xi_cons, K_res
# end

function compute_default_thresholds(A::AbstractMatrix, epsilon::AbstractMatrix, R::Int; margin::Float64=1.0)
    S = size(A,1)
    C = S - R

    xi_cons = zeros(Float64, C)
    for k in 1:C
        i = R + k
        gain = sum(epsilon[i,j]*A[i,j] for j in 1:R if A[i,j] > 0; init =0.0)
        loss = sum(-A[i,j]    for j in R+1:S if A[i,j] < 0; init =0.0)
        xi_cons[k] = (gain + loss)*(1 + margin)
    end

    K_res = zeros(Float64, R)
    for i in 1:R
        drain = sum(-A[i,j] for j in R+1:S if A[i,j] < 0; init = 0.0)
        K_res[i] = drain*(1 + margin)
    end

    return xi_cons, K_res
end


# ────────────────────────────────────────────────────────────────────────────────
# 2) Given Xi and K build & solve M B_eq = b, check positivity
# ────────────────────────────────────────────────────────────────────────────────
"""
    calibrate_from_K_xi(
      K_res::Vector{<:Real},    # length R: resource carrying capacities
      xi_cons::Vector{<:Real},  # length C: consumer thresholds
      epsilon::AbstractMatrix,        # SS conversion efficiencies
      A::AbstractMatrix         # SS interaction signs/magnitudes
    ) -> (R_eq, C_eq) or nothing

Solve the steady‐state abundances B = [R_eq; C_eq] from:

  For i=1:R (resources):
    0 = -K_i - B_i  - ∑_{j=R+1}^S A_{j,i} * B_j

  For i=R+1:S (consumers, indexed k=i-R):
    0 = -Xi_k - B_i + ∑_{j=1}^R epsilon_{i,j}*A_{i,j}*B_j - ∑_{j=1}^R A_{j,i}*B_j

Returns nothing if no positive solution exists.
"""
function calibrate_from_K_xi(
    xi_cons::AbstractVector, K_res::AbstractVector,
    epsilon::AbstractMatrix, A::AbstractMatrix
)
    R = length(K_res)
    S = size(A,1)
    C = S - R

    M = zeros(Float64, S, S)
    b = zeros(Float64, S)

    # resource rows
    for i in 1:R
        for j in 1:S
            M[i,j] = (i==j ? 1.0 : 0.0) - A[i,j]
        end
        b[i] = K_res[i]
    end

    # consumer rows
    for k in 1:C
        i = R + k
        for j in 1:S
            if j == i
                M[i,j] = 1.0
            elseif A[i,j] > 0
                M[i,j] = -epsilon[i,j]*A[i,j]   # prey term
            elseif A[i,j] < 0
                M[i,j] = -A[i,j]          # predation loss
            else
                M[i,j] = 0.0
            end
        end
        b[i] = -xi_cons[k]
    end

    return M \ b
end

function sample_Xi_K_from_A(A::AbstractMatrix, epsilon::AbstractMatrix, R::Int;
                            scale_factor::Float64=1.0)
    S = size(A,1)
    C = S - R

    # for each resource i, base K_i = total drain by consumers
    K = zeros(R)
    for i in 1:R
        # drain = sum(abs(A[i,j]) for j in 1:R; init = 0.0)
        K[i]  = abs(rand(Normal(2.0, 2.0))) # drain * scale_factor * rand(Uniform(0.5, 1.5))
    end

    # for each consumer k, base Xi_k = total potential gain+loss
    Xi = zeros(C)
    for k in 1:C
        i = R + k
        # raw_gain  = sum(epsilon[i,j]*abs(A[i,j]) for j in R+1:S if A[i,j] > 0; init = 0.0)
        # raw_loss  = sum(abs(A[i,j]) for j in R+1:S if A[i,j] < 0; init = 0.0)
        # Xi[k] = (raw_gain + raw_loss) * scale_factor * rand(Uniform(0.5, 1.5))
        Xi[k] = abs(rand(Normal(1.0, 1.0)))
    end

    return K, Xi
end

# ────────────────────────────────────────────────────────────────────────────────
# 3) generate_feasible_thresholds
# ────────────────────────────────────────────────────────────────────────────────
"""
For each margin delt in `margins`, form

    Xi_cons, K_res = compute_default_thresholds(A, epsilon, R; margin=delt)

and keep those where

1. `M / b` is finite & positive,
2. the Jacobian at that equilibrium is locally stable.

Returns a vector of
  (xi_cons=..., K_res=..., margin=delt, R_eq=..., C_eq=...)
for every delt that passes these checks.
"""
function generate_feasible_thresholds(
    A::AbstractMatrix, epsilon::AbstractMatrix, R::Int;
    margins = 1.0
)
    S = size(A,1)
    C = S - R
    out = NamedTuple[]
    for delt in margins
        # 1) propose thresholds
        xi_cons, K_res = compute_default_thresholds(A, epsilon, R; margin=delt)

        # 2) try solve equilibrium
        eq = try
            calibrate_from_K_xi(xi_cons, K_res, epsilon, A)
        catch
            continue
        end

        # 3) positivity & finiteness
        if any(!isfinite, eq) || any(x->x<=0, eq)
            continue
        end

        # 4) stability
        R_eq = view(eq, 1:R)
        C_eq = view(eq, R+1:S)
        p    = (R, C, xi_cons, xi_cons, K_res, ones(R), epsilon, A)
        D, M = compute_jacobian(vcat(R_eq, C_eq), p)
        J    = D*M
        if any(x -> !isfinite(x), eq) || any(x -> x <= 0, eq)
            continue
        end

        # 5) accept
        push!(out, (
            xi_cons = xi_cons,
            K_res    = K_res,
            margin   = delt,
            R_eq     = copy(R_eq),
            C_eq     = copy(C_eq),
        ))
    end

    return out
end
#############################################################################
#############################################################################
############################# make_A ########################################
#############################################################################
#############################################################################
function make_A(
    A::AbstractMatrix, R::Int, conn::Float64, scenario::Symbol;
    IS = 1.0,
    pareto_exponent::Float64=1.75,
    pareto_minimum_degree=1.0,
    mod_gamma::Float64=5.0,
    B_term::Bool=false,
    B_term_IS::Float64=0.1
)
    S = size(A,1)
    C = S - R
    fill!(A, 0.0)

    if scenario == :ER
        for i in (R+1):S, j in 1:R
            if rand() < conn && i != j
                A[i,j] = abs(rand(Normal(0.0, IS)))
                A[j,i] = -abs(rand(Normal(0.0, IS)))
            end
        end
        
    elseif scenario == :PL
        raw = rand(Pareto(pareto_minimum_degree, pareto_exponent), C)
        ks  = clamp.(floor.(Int, raw), 1, R)
        for (idx,k) in enumerate(ks)
            i = R + idx
            for j in sample(1:R, min(k,R); replace=false)
                if i != j
                    A[i,j] = abs(rand(Normal(0.0, IS)))
                    A[j,i] = -abs(rand(Normal(0.0, IS)))
                end
            end
        end

    elseif scenario == :MOD
        halfR, halfC = fld(R,2), fld(C,2)
        res1, res2   = 1:halfR, (halfR+1):R
        con1, con2   = (R+1):(R+halfC), (R+halfC+1):S

        for i in (R+1):S, j in 1:R
            same = (i in con1 && j in res1) || (i in con2 && j in res2)
            p    = same ? conn*mod_gamma : conn/mod_gamma
            if rand() < clamp(p,0,1) && i != j
                A[i,j] = abs(rand(Normal(0.0, IS)))
                A[j,i] = -abs(rand(Normal(0.0, IS)))
            end
        end

    else
        error("Unknown scenario $scenario")
    end

    # optionally sprinkle in consumer-consumer predation
    if B_term
        for i in (R+1):S, j in (R+1):S
            if i != j && rand() < conn
                A[i,j] = abs(rand(Normal(0.0, B_term_IS)))
                A[j,i] = -abs(rand(Normal(0.0, B_term_IS)))
            end
        end
    end

    return A
end
#############################################################################
#############################################################################
################## JACOBIAN, RESILIENCE AND RESISTANCE ######################
#############################################################################
#############################################################################
function compute_jacobian(B, p)
    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p
    S = R + C

    # 1) Build D
    D = zeros(S,S)
    for i in 1:R
        D[i,i] = r_res[i]/K_res[i] * B[i]
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

# ────────────────────────────────────────────────────────────────────────────────
# 0) Helper: compute the collectivity parameter phi for one (A, epsilon)
# ────────────────────────────────────────────────────────────────────────────────
"""
    compute_collectivity(A, epsilon)

Given the S S interaction matrix A and efficiency matrix epsilon, build the 
non‐dimensional direct‐interaction matrix A_prime with

  A_prime[i,j] =  epsilon[i,j]*A[i,j]   if A[i,j] > 0  (consumer gains)
           A[i,j]          if A[i,j] < 0  (losses)
           0.0             if A[i,j] == 0

and return phi = max(abs, eigvals(A_prime)).
"""
function compute_collectivity(A::AbstractMatrix, epsilon::AbstractMatrix)
    S = size(A,1)
    A_prime = zeros(eltype(A), S, S)
    for i in 1:S, j in 1:S
        if     A[i,j] >  0
            A_prime[i,j] =  epsilon[i,j] * A[i,j]
        elseif A[i,j] <  0
            A_prime[i,j] =  A[i,j]
        else
            A_prime[i,j] =  0.0
        end
    end
    # spectral radius
    vals = eigvals(A_prime)
    return maximum(abs, vals)
end

#############################################################################
#############################################################################
######################## COMPUTE COLLECTIVITY ###############################
#############################################################################
#############################################################################
# ────────────────────────────────────────────────────────────────────────────────
# 0) Helper: compute the collectivity parameter phi for one (A, epsilon)
# ────────────────────────────────────────────────────────────────────────────────
"""
    compute_collectivity(A, epsilon)

Given the S S interaction matrix A and efficiency matrix epsilon, build the 
non‐dimensional direct‐interaction matrix A_prime with

  A_prime[i,j] =  epsilon[i,j]*A[i,j]   if A[i,j] > 0  (consumer gains)
           A[i,j]          if A[i,j] < 0  (losses)
           0.0             if A[i,j] == 0

and return phi = max(abs, eigvals(A_prime)).
"""
function compute_collectivity(A::AbstractMatrix, epsilon::AbstractMatrix)
    S = size(A,1)
    A_prime = zeros(eltype(A), S, S)
    for i in 1:S, j in 1:S
        if     A[i,j] >  0
            A_prime[i,j] =  epsilon[i,j] * A[i,j]
        elseif A[i,j] <  0
            A_prime[i,j] =  A[i,j]
        else
            A_prime[i,j] =  0.0
        end
    end
    # spectral radius
    vals = eigvals(A_prime)
    return maximum(abs, vals)
end

#############################################################################
#############################################################################
################## SIMULATE PRESS PERTURBATIONS #############################
#############################################################################
#############################################################################
function simulate_press_perturbation(
    u0, p, tspan, t_perturb, delta;
    solver=AutoTsit5(Rosenbrock23()),
    plot=false,
    show_warnings=true,
    full_or_simple=true,
    cb = cb,
    species_specific_perturbation=false
)
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p

    # --- Phase 1: run up to t_perturb ---
    prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], t_perturb), p)
    sol1 = show_warnings ? solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8) :
                           with_logger(logger) do
                               solve(prob1, solver; callback = cb, reltol=1e-8, abstol=1e-8)
                           end
    pre_state = sol1.u[end]
    before_persistence = count(x -> x > EXTINCTION_THRESHOLD, pre_state) / length(pre_state)
    # pre_state[pre_state .> EXTINCTION_THRESHOLD] .= 0.0
    
    # @warn "Sol1 is $(sol1[:, :])"
    ########## FIRST PART: WHOLE COMMUNITY PERTURBATION ##########
    # --- Phase 2: apply press (reduce thresholds by delta) ---
    xi_press = xi_cons .* (1 + delta)  # perturb all species
    K_press = K_res .* (1 - delta)  # perturb all species
    
    # xi_press = copy(xi_cons)     
    # K_press = copy(K_res)
    # K_press[1] = K_press[1] * (1 - MAX_DELTA)  # perturb only resource 1
    
    p_press  = (R, C, m_cons, xi_press, r_res, K_press, epsilon, A)
    prob2    = ODEProblem(trophic_ode!, pre_state, (t_perturb, tspan[2]), p_press)
    sol2     = solve(prob2, solver; callback = cb, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    n = length(new_equil)
    after_persistence = count(x -> x > EXTINCTION_THRESHOLD, new_equil) / length(new_equil)
    # if new_equil[1] > EXTINCTION_THRESHOLD
    #     @warn "Resource 1 did not go extinct after perturbation, but should have!. We're in simulate_press_perturbation, it's whole process is:"
    #     @warn "$(sol2[:, :])"
    #     @warn "And delta was $(delta)"
    #     error()
    # end

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
    
            p_press  = (R, C, m_cons, xi_press, r_res, K_res, epsilon, A)
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
        fig = Figure(; size=(1600,650))
        ax1 = Axis(fig[1,1];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "Before Press (Full Model)" :
                                           "Before Press (Simplified)")
        ax2 = Axis(fig[1,2];
                   xlabel="Time", ylabel="Biomass",
                   title= full_or_simple ? "After Press (Full Model) delta=$(delta)" :
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
#############################################################################
#############################################################################
################## SIMULATE PULSE PERTURBATIONS #############################
#############################################################################
#############################################################################
"""
    simulate_pulse_perturbation(
        u0, p, tspan, t_pulse, delt;
        solver=AutoTsit5(Rosenbrock23()),
        plot=false,
        cb=nothing,
        species_specific_perturbation=false
    ) -> (return_times, before_persist, after_persist, eq_state, species_rt)

- `u0` : initial abundances (length S)
- `p`  : parameters tuple `(R, C, m_cons, xi_cons, r_res, K_res, epsilon, A)`
- `tspan` : `(t0, tfinal)`
- `t_pulse` : time at which to apply the pulse
- `delt`  : fractional pulse (e.g. 0.2 knocks everyone up by +20%; use `-0.5` to reduce by 50%)
- `cb` : a `CallbackSet` for extinction, defaults to none
- `species_specific_perturbation` : if `true`, also gives per-consumer pulse return-times

Returns:
1. `return_times::Vector{Float64}` - time for each of the S species to return within 10% of equilibrium  
2. `before_persist::Float64` - fraction of species alive just before the pulse  
3. `after_persist::Float64`  - fraction of species alive at the end  
4. `eq_state::Vector{Float64}` - the asymptotic state after the pulse  
5. `species_rt::Vector{Float64}` (length C) - if `species_specific_perturbation`, return-times for each consumer pulse
"""
function simulate_pulse_perturbation(
    u0, p, tspan, t_pulse, delt;
    solver=AutoTsit5(Rosenbrock23()),
    plot=false,
    cb=nothing,
    species_specific_perturbation=false
)
    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p
    S = R + C

    # 1) Integrate up to the pulse
    prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], t_pulse), p)
    sol1 = solve(prob1, solver; callback=cb, reltol=1e-8, abstol=1e-8)
    pre_state = sol1.u[end]
    before_persist = count(x->x>EXTINCTION_THRESHOLD, pre_state) / S

    # 2) Apply the pulse: fractional change delt
    pulsed = pre_state .* (1 + delt)

    # 3) Integrate from t_pulse to tfinal with same parameters
    prob2 = ODEProblem(trophic_ode!, pulsed, (t_pulse, tspan[2]), p)
    sol2 = solve(prob2, solver; callback=cb, reltol=1e-8, abstol=1e-8)

    # equilibrium after pulse (assumed to be sol2.u[end])
    eq_state = sol2.u[end]
    after_persist = count(x->x>EXTINCTION_THRESHOLD, eq_state) / S

    # 4) Return times: when each species returns within 10% of eq_state
    return_times = fill(NaN, S)
    for i in 1:S
        target = eq_state[i]
        for (t, u) in zip(sol2.t, sol2.u)
            if abs(u[i] - target) / (abs(target) + 1e-8) < 0.10
                return_times[i] = t - t_pulse
                break
            end
        end
    end

    # 5) Species‐specific pulse (only consumers) if requested
    species_rt = fill(NaN, C)
    if species_specific_perturbation
        for k in 1:C
            i = R + k
            # perturb only consumer k
            u_k = copy(pre_state)
            u_k[i] *= (1 + delt)
            prob_k = ODEProblem(trophic_ode!, u_k, (t_pulse, tspan[2]), p)
            sol_k  = solve(prob_k, solver; callback=cb, reltol=1e-8, abstol=1e-8)
            final = sol_k.u[end][i]
            # find return time for species i
            for (t, u) in zip(sol_k.t, sol_k.u)
                if abs(u[i] - final) / (abs(final) + 1e-8) < 0.10
                    species_rt[k] = t - t_pulse
                    break
                end
            end
        end
    end

    # 6) Optional plotting
    if plot
        fig = Figure(; size=(1200,600))
        ax1 = Axis(fig[1,1], title="Pre-Pulse Dynamics", xlabel="Time", ylabel="Biomass")
        ax2 = Axis(fig[1,2], title="Post-Pulse Dynamics", xlabel="Time", ylabel="Biomass")
        for i in 1:R
            lines!(ax1, sol1.t, [u[i] for u in sol1.u], color=:blue)
            lines!(ax2, sol2.t, [u[i] for u in sol2.u], color=:blue)
        end
        for i in R+1:S
            lines!(ax1, sol1.t, [u[i] for u in sol1.u], color=:red)
            lines!(ax2, sol2.t, [u[i] for u in sol2.u], color=:red)
        end
        display(fig)
    end

    return return_times, before_persist, after_persist, eq_state, species_rt
end

##########################################################################
##########################################################################
##########################################################################
################## TARGETTED PERTURBATIONS ###############################
##########################################################################
##########################################################################
const TOL = 1e-2    # precision on δK/K or δxi/xi
const MAX_DELTA = 0.99 # search up to 100% change

"""
    find_min_K_reduction(i, u0, p, tspan, tpert; cb)

For resource **i** (1 ≤ i ≤ R), finds the smallest δ∈[0,MAX_DELTA] such that reducing K_res[i] by
(i.e. K_res[i]→K_res[i]*(1-δ)) causes species i to go extinct (below EXTINCTION_THRESHOLD) by tspan[2].
Returns δ_crit (or `nothing` if even δ=MAX_DELTA does not kill it).
"""
function find_min_K_reduction(
    i::Int,
    u0::AbstractVector,
    p,
    tspan::Tuple{<:Real,<:Real},
    tpert::Real;
    cb
)

    ##################################################################################
    # Unpack parameters
    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p

    # --- Phase 1: run up to t_perturb ---
    prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], tpert), p)
    sol1 = solve(prob1, AutoTsit5(Rosenbrock23()); callback = cb, reltol=1e-6, abstol=1e-6)
    pre_state = sol1.u[end]
    
    ########## CHECKING MAXIMUM PERTURBATION ##########
    lo, hi = 0.0, MAX_DELTA

    # --- Phase 2: apply press (reduce thresholds by delta) ---
    xi_press = copy(xi_cons)     
    K_press = copy(K_res)
    K_press[i] = K_press[i] * (1 - hi)  # perturb only resource 1
    
    p_press  = (R, C, m_cons, xi_press, r_res, K_press, epsilon, A)
    prob2    = ODEProblem(trophic_ode!, pre_state, (tpert, tspan[2]), p_press)
    sol2     = solve(prob2, AutoTsit5(Rosenbrock23()); callback = cb, reltol=1e-8, abstol=1e-8)
    new_equil = sol2.u[end]
    if new_equil[i] > EXTINCTION_THRESHOLD
        # @warn "Resource $i does not go extinct even with maximum K reduction. Its final vector is:"
        # @warn "$(sol2[i, :])"
        # error()
        return 0.0  # already extinct
    end
    
    # bisection:
    while hi - lo > TOL
        mid = (lo+hi)/2
        # --- Phase 2: apply press (reduce thresholds by delta) ---
        xi_press = copy(xi_cons)     
        K_press = copy(K_res)
        K_press[i] = K_press[i] * (1 - mid)  # perturb only resource 1
        p_press  = (R, C, m_cons, xi_press, r_res, K_press, epsilon, A)
        # integrate:
        # --- Phase 1: run up to t_perturb ---
        prob1 = ODEProblem(trophic_ode!, u0, (tspan[1], tpert), p)
        sol1 = solve(prob1, AutoTsit5(Rosenbrock23()); callback = cb, reltol=1e-6, abstol=1e-6)
        pre_state = sol1.u[end]
        prob2    = ODEProblem(trophic_ode!, pre_state, (tpert, tspan[2]), p_press)
        sol2     = solve(prob2, AutoTsit5(Rosenbrock23()); callback = cb, reltol=1e-8, abstol=1e-8)
        new_equil = sol2.u[end]
        # --- Phase 2: apply press (reduce thresholds by delta) ---
        if new_equil[i] <= EXTINCTION_THRESHOLD
            hi = mid
        else
            lo = mid
        end
    end

    return hi
end

function find_min_K_reduction(i::Int,
                              u0::AbstractVector,
                              p::NamedTuple{(:R,:C,:m_cons,:xi_cons,:r_res,:K_res,:epsilon,:A),T},
                              tspan::Tuple{<:Real,<:Real},
                              tpert::Real;
                              cb) where T

    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p
       
    lo, hi = 0.0, MAX_DELTA*2

    K0 = copy(p.K_res)
    K0[i] = K0[i] * (1 - hi)  # perturb only resource i

    # check if even max does not kill:
    p_hi = (R, C, m_cons, xi_cons, r_res, K0, epsilon, A)
    sol = solve(ODEProblem(trophic_ode!, u0, (tspan[1],tpert), p), AutoTsit5(Rosenbrock23()); callback=cb)
    post_hi = solve(ODEProblem(trophic_ode!, sol.u[end], (tpert,tspan[2]), p_hi), AutoTsit5(Rosenbrock23()); callback=cb).u[end]
    if post_hi[i] > 0.01 #EXTINCTION_THRESHOLD
        error("Resource $i does not go extinct even with maximum K reduction. Its final vector is: $(post_hi[i])")
        return nothing
    end

    # bisection:
    while hi - lo > TOL
        mid = (lo+hi)/2
        p_mid = merge(p, (; K_res = K0 .* (1 .- mid) ))
        # integrate:
        sol1 = solve(ODEProblem(trophic_ode!, u0, (tspan[1],tpert), p), AutoTsit5(Rosenbrock23()); callback=cb)
        post = solve(ODEProblem(trophic_ode!, sol1.u[end], (tpert,tspan[2]), p_mid),
                     AutoTsit5(Rosenbrock23()); callback=cb).u[end]
        if post[i] <= 0.01 #EXTINCTION_THRESHOLD
            hi = mid
        else
            lo = mid
        end
    end

    return hi
end

"""
    find_min_xi_increase(j, u0, p, tspan, tpert; cb)

For **consumer** j (1 ≤ j ≤ C), finds the smallest δ∈[0,MAX_DELTA] such that increasing xi_cons[j] by
(i.e. xi_cons[j]→xi_cons[j]*(1+δ)) drives species R+j extinct by tspan[2].  Returns δ_crit or `nothing`.
"""
function find_min_xi_increase(j::Int,
                              u0::AbstractVector,
                              p,
                              tspan::Tuple{<:Real,<:Real},
                              tpert::Real;
                              cb)

    R, C, m_cons, xi_cons, r_res, K_res, epsilon, A = p
    
    lo, hi = 0.0, MAX_DELTA

    xi0 = copy(xi_cons)
    xi0[j] = xi0[j] * (1 + hi)  # perturb only consumer j
    # check if even max does not kill:
    p_hi = (R, C, m_cons, xi0, r_res, K_res, epsilon, A)
    sol = solve(ODEProblem(trophic_ode!, u0, (tspan[1],tpert), p), AutoTsit5(Rosenbrock23()); callback=cb)
    post_hi = solve(ODEProblem(trophic_ode!, sol.u[end], (tpert,tspan[2]), p_hi), AutoTsit5(Rosenbrock23()); callback=cb).u[end]
    if post_hi[R+j] > EXTINCTION_THRESHOLD
        # @warn "Consumer $j does not go extinct even with maximum xi increase. Its final vector is: $(post_hi[R+j])"
        return 0.0
    end

    # bisection:
    while hi - lo > TOL
        mid = (lo+hi)/2
        xi_mid = copy(xi0)
        xi_mid[j] *= 1 + mid
        p_mid = (R, C, m_cons, xi_mid, r_res, K_res, epsilon, A)
        sol1 = solve(ODEProblem(trophic_ode!, u0, (tspan[1],tpert), p), AutoTsit5(Rosenbrock23()); callback=cb)
        post = solve(ODEProblem(trophic_ode!, sol1.u[end], (tpert,tspan[2]), p_mid),
                     AutoTsit5(Rosenbrock23()); callback=cb).u[end]
        if post[R+j] <= EXTINCTION_THRESHOLD
            hi = mid
        else
            lo = mid
        end
    end

    return hi
end

##########################################################################
##########################################################################
################### Ladder-of-simplification transforms ##################
##########################################################################
##########################################################################
function transform_for_ladder_step(step, A_adj, epsilon_full)
    total = size(A_adj,1)

    if step == 1
        # - no change -
        return A_adj, epsilon_full
    elseif step == 2
        # S2: per‐row mean epsilon
        epsilon2 = similar(epsilon_full)
        for i in 1:total
            epsilon2[i, :] .= mean(epsilon_full[i, :])
        end
        return A_adj, epsilon2

    elseif step == 3
        # S3: global mean epsilon
        mu = mean(epsilon_full)
        return A_adj, fill(mu, total, total)

    elseif step == 4
        # S4: re-randomize epsilon from Normal(mean,std)
        mu, sigma = mean(epsilon_full), std(epsilon_full)
        epsilon4 = similar(epsilon_full)
        for i in 1:total, j in 1:total
            epsilon4[i,j] = clamp(rand(Normal(mu,sigma)), 0, 1)
        end
        return A_adj, epsilon4

    elseif 5 <= step <= 8
        # S5-S8: replace A by its average positive & negative values
        pos = mean(A_adj[A_adj .> 0])
        neg = mean(A_adj[A_adj .< 0])
        A5  = map(x -> x>0 ? pos : x<0 ? neg : 0, A_adj)

        if step == 5
            return A5, epsilon_full
        elseif step == 6
            # per‐row mean epsilon
            epsilon6 = similar(epsilon_full)
            for i in 1:total
                epsilon6[i, :] .= mean(epsilon_full[i, :])
            end
            return A5, epsilon6
        elseif step == 7
            mu = mean(epsilon_full)
            return A5, fill(mu, total, total)
        else  # step == 8
            mu, sigma = mean(epsilon_full), std(epsilon_full)
            epsilon8 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon8[i,j] = clamp(rand(Normal(mu,sigma)), 0, 1)
            end
            return A5, epsilon8
        end

    elseif 9 <= step <= 12
        # S9-S12: set every nonzero A entry to the mean magnitude m
        m   = mean(abs.(A_adj[A_adj .!= 0]))
        A6  = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)

        if step == 9
            return A6, epsilon_full
        elseif step == 10
            epsilon10 = similar(epsilon_full)
            for i in 1:total
                epsilon10[i, :] .= mean(epsilon_full[i, :])
            end
            return A6, epsilon10
        elseif step == 11
            mu = mean(epsilon_full)
            return A6, fill(mu, total, total)
        else  # step == 12
            mu, sigma = mean(epsilon_full), std(epsilon_full)
            epsilon12 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon12[i,j] = clamp(rand(Normal(mu,sigma)), 0, 1)
            end
            return A6, epsilon12
        end

    elseif 13 <= step <= 16
        # S13-S16: fully re-randomize A over the same sparsity pattern
        vals   = abs.(A_adj[A_adj .!= 0])
        mu, sigma   = mean(vals), std(vals)
        A_rand = zeros(eltype(A_adj), total, total)
        for i in 1:total, j in 1:total
            if A_adj[i,j] != 0
                A_rand[i,j] = abs(rand(Normal(mu,sigma))) * sign(A_adj[i,j])
            end
        end

        if step == 13
            return A_rand, epsilon_full
        elseif step == 14
            epsilon14 = similar(epsilon_full)
            for i in 1:total
                epsilon14[i, :] .= mean(epsilon_full[i, :])
            end
            return A_rand, epsilon14
        elseif step == 15
            mu = mean(epsilon_full)
            return A_rand, fill(mu, total, total)
        elseif step == 16
            mu, sigma = mean(epsilon_full), std(epsilon_full)
            epsilon16 = similar(epsilon_full)
            for i in 1:total, j in 1:total
                epsilon16[i,j] = clamp(rand(Normal(mu,sigma)), 0, 1)
            end
            return A_rand, epsilon16
        end
    elseif step == 17
        # S17: global mean of epsilon*A
        pos_vals = [
            epsilon_full[i,j]*A_adj[i,j] 
            for i in 1:total, j in 1:total if A_adj[i,j] > 0 
        ]
        neg_vals = abs.([A_adj[i,j] for i in 1:total, j in 1:total if A_adj[i,j] < 0])
        mu = mean(vcat(pos_vals, neg_vals))

        # build A17: mu in the same sparsity pattern, zeros preserved
        A17 = map(x -> x>0 ? mu : x<0 ? -mu : 0.0, A_adj)

        # epsilon17 = all ones
        epsilon17 = ones(eltype(epsilon_full), total, total)

        return A17, epsilon17
    else
        error("Unknown ladder step $step")
    end
end

function short_transform_for_ladder_step(step, A_adj, epsilon_full)
    total = size(A_adj,1)

    if step == 1
        # - no change -
        return A_adj, epsilon_full
    elseif step == 2
        # S9-S12: set every nonzero A entry to the mean magnitude m
        m   = mean(abs.(A_adj[A_adj .!= 0]))
        A2  = ifelse.(A_adj .!= 0, m*sign.(A_adj), 0.0)

        mu = mean(epsilon_full)
        return A2, fill(mu, total, total)
    elseif step == 3
        # S17: global mean of epsilon*A
        pos_vals = [
            epsilon_full[i,j]*A_adj[i,j] 
            for i in 1:total, j in 1:total if A_adj[i,j] > 0 
        ]
        neg_vals = abs.([A_adj[i,j] for i in 1:total, j in 1:total if A_adj[i,j] < 0])
        mu = mean(vcat(pos_vals, neg_vals))

        # build A17: mu in the same sparsity pattern, zeros preserved
        A3 = map(x -> x>0 ? mu : x<0 ? -mu : 0.0, A_adj)

        # epsilon17 = all ones
        epsilon3 = ones(eltype(epsilon_full), total, total)

        return A3, epsilon3
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

# """
#     median_return_rate(J::AbstractMatrix, B::AbstractVector;
#                        t::Real, n::Int=1000, rng::AbstractRNG = Random.GLOBAL_RNG)
    
# For Jacobian J and equilibrium B (length S), draw `n` random pulses
# u[k] with entries ∝ B.^2, normalize each to unit norm, then compute
# alpha_k = |exp(t*J)*u[k]|.  Returns
# R_med(t) = -(1/t) * log(median(alpha_k)).
# """
function median_return_rate(J, B; t=0.01, n=1000, rng=Random.GLOBAL_RNG)
    S = length(B)
    # precompute matrix exponential
    E = exp(t*J)
    alpha = Vector{Float64}(undef, n)

    # draw pulses
    for k in 1:n
        u = randn(rng, S)# .* (B .^ 2)         # biomass-weighted
        u ./= norm(u)                         # unit length
        alpha[k] = norm(E * u)                    # growth factor
    end

    return -log(median(alpha)) / t
end

# function species_return_rates(J::AbstractMatrix, B::AbstractVector; t=0.01, n=1000, rng=Random.GLOBAL_RNG)
#     S = length(B)
#     if any(!isfinite, J)
#         return fill(NaN, S)
#     end

#     E = exp(t*J)
    
#     absV = Matrix{Float64}(undef, S, n)

#     for k in 1:n
#         u = randn(rng, S)# .* (B .^ 2)
#         u ./= norm(u)
#         v = E * u
#         @inbounds absV[:, k] = abs.(v)
#     end

#     # now median over the 2nd dimension, and convert to rates
#     med = mapslices(median, absV; dims=2)[:, 1]
#     return @. -(1/t) * log(med)
# end

function species_return_rates(J::AbstractMatrix, B::AbstractVector;
                              t::Real=0.01, n::Int=1000,
                              rng::AbstractRNG=Random.GLOBAL_RNG)
    S = length(B)
    if any(!isfinite, J)
        return fill(NaN, S)
    end

    # precompute matrix exponential
    E = exp(t*J)

    # store per‐species alpha values
    alpha = Matrix{Float64}(undef, S, n)
    for k in 1:n
        # draw a random perturbation
        u = randn(rng, S)# .* (B .^ 2)
        u ./= norm(u)
        v = E * u
        @inbounds for i in 1:S
            alpha[i,k] = abs(v[i] / u[i])
        end
    end

    # median over trials, then convert to a rate
    med = mapslices(median, alpha; dims=2)[:,1]
    return -log.(med) ./ t
end

"""
    analytical_species_return_rates(J::AbstractMatrix; t::Real=0.01)

Compute the exact, species-level return rate vector
r_i(t) = [exp(t * J) * exp(t * J)']_{i,i} / (2 * t),
for i = 1:S, where S = size(J,1).  Returns an S-element Vector{Float64}.

In the short-time limit (t -> 0), this should approach -diag(J).
"""
function analytical_species_return_rates(J::AbstractMatrix; t::Real = 0.01)
    S = size(J, 1)
    # If J contains any non-finite values, return NaNs
    if any(!isfinite, J)
        return fill(NaN, S)
    end

    # Compute exp(t * J)
    E = exp(t * J)

    # Compute the product exp(tJ) * (exp(tJ))'
    M = E * transpose(E)

    # Extract the diagonal and divide by 2*t
    rates = Vector{Float64}(undef, S)
    for i in 1:S
        rates[i] = (M[i,i] - 1.0) / (2*t)
    end
    return -rates # flip the sign to make a positive return rate?? Is this correct?
end

"""
    analytical_median_return_rate(J; t=0.01)

Computes the median across species of your fully-analytical
species-level return-rates, i.e. the direct analogue of
median_return_rate but without Monte Carlo.
"""
function analytical_median_return_rate(J::AbstractMatrix; t::Real=0.01)
    rates = analytical_species_return_rates(J; t=t)
    return median(rates)
end

"""
    degree_scaled_A(A, R, scen)

If scen == :PL, rescales each consumer's outgoing links and each resource's incoming
links by 1/degree, where "degree" is the number of non-zero links.  Otherwise returns A unchanged.
"""
function degree_scaled_A(A::Matrix{Float64}, R::Int, scen::Symbol)
    S = size(A,1)
    if scen != :PL
        return A
    end

    A2 = copy(A)

    # 1) Resources 1:R: count how many predators each has, then scale its column
    for res in 1:R
        deg = count(i -> A2[i,res] != 0.0, R+1:S)
        if deg > 0
            α = 1/deg
            for i in R+1:S
                if A2[i,res] != 0.0
                    A2[i,res] *= α      # consumer resource gain
                    A2[res,i] *= α      # resource consumer loss
                end
            end
        end
    end

    # 2) Consumers R+1:S: count how many preys each has, then scale its row
    for con in (R+1):S
        deg = count(j -> A2[con,j] != 0.0, 1:R)
        if deg > 0
            α = 1/deg
            for j in 1:R
                if A2[con,j] != 0.0
                    A2[con,j] *= α      # consumer resource gain
                    A2[j,con] *= α      # resource consumer loss
                end
            end
        end
    end

    return A2
end