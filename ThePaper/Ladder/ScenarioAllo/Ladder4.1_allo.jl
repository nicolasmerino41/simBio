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
function allo_ode!(du, u, p, t)
    S       = p.S
    r       = p.r
    K       = p.K
    x       = p.x
    d       = p.d
    attack  = p.attack
    epsilon = p.epsilon
    B0_ref  = p.B0_ref
    i0      = p.i0
    h       = p.h

    # zero out
    @inbounds for i in 1:S
        du[i] = 0
    end

    Bh = u .^ h

    @inbounds for i in 1:S
        # intrinsic + metabolism + mortality
        Gi = r[i]>0 ? (1 - u[i]/K[i]) : 0
        du[i] = r[i]*u[i]*Gi - x[i]*u[i] - d[i]*u[i]

        gain = 0.0
        loss = 0.0

        for j in 1:S
            if attack[i,j] > 0
                # functional response denominator
                denom = B0_ref[i]^h + i0[i]*Bh[i] + sum(attack[i,k]*Bh[k] for k in 1:S if attack[i,k]>0)
                Fij   = attack[i,j]*Bh[j] / denom
                gain += x[i]*u[i]*Fij / epsilon[i,j]
            end
            if attack[j,i] > 0
                denom = B0_ref[j]^h + i0[j]*Bh[j] + sum(attack[j,k]*Bh[k] for k in 1:S if attack[j,k]>0)
                Fji   = attack[j,i]*Bh[i] / denom
                loss += x[j]*u[j]*Fji
            end
        end

        du[i] += gain - loss
    end
end

#############################################################################
#############################################################################
########################### is_locally_stable and survives ##################
#############################################################################
function is_locally_stable(J::AbstractMatrix)
    if any(!isfinite, J)
        return false
    end
    lambda = eigvals(J)
    maximum(real.(lambda)) < 0 
end
function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(allo_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false,sol.u[end]) : (all(sol.u[end] .> EXTINCTION_THRESHOLD),sol.u[end])
end
#############################################################################
#############################################################################
################## CALIBRATE PARAMETRISATION ################################
#############################################################################
#############################################################################
function compute_default_thresholds(A::AbstractMatrix, epsilon::AbstractMatrix, R::Int; margin::Float64=10.0)
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
    mod_gamma::Float64=5.0,
    B_term::Bool=false
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
        raw = rand(Pareto(1.0, pareto_exponent), C)
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

    # optionally sprinkle in consumer→consumer predation
    if B_term
        for i in (R+1):S, j in (R+1):S
            if i != j && rand() < conn
                A[i,j] = abs(rand(Normal(0.0, IS*0.1)))
                A[j,i] = -abs(rand(Normal(0.0, IS*0.1)))
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
using ForwardDiff

function compute_jacobian_allo(Bstar::AbstractVector, p)
    # g(u) returns f(u)---the time‐derivative at state u
    g(u) = begin
        du = similar(u)
        allo_ode!(du, u, p, 0.0)
        return du
    end

    # allocate output
    S = length(Bstar)
    J = Matrix{Float64}(undef, S, S)

    # fill J in‐place
    ForwardDiff.jacobian!(J, g, Bstar)

    return J
end

# Resilience: negative of the largest real part of the Jacobian eigenvalues.
function compute_resilience(B, p)
    D, Mstar = compute_jacobian_allo(B, p)
    J = D * Mstar
    ev = eigvals(J)
    return maximum(real.(ev))
end

# Reactivity: maximum eigenvalue of the symmetric part of the Jacobian.
function compute_reactivity(B, p)
    D, Mstar = compute_jacobian_allo(B, p)
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
function simulate_press_perturbation_allo!(
    u0, p, tspan, tpert, delta;
    solver = Tsit5(),
    cb
)
    S       = p.S
    r       = p.r
    K       = p.K
    x       = p.x
    d       = p.d
    attack  = p.attack
    epsilon       = p.epsilon
    B0_ref  = p.B0_ref
    i0      = p.i0
    h       = p.h

    # Phase 1
    prob1 = ODEProblem(allo_ode!, u0, (tspan[1], tpert), p)
    sol1  = solve(prob1, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    pre   = sol1.u[end]
    before = count(>(EXTINCTION_THRESHOLD), pre)/S

    # Phase 2: press on K
    K2 = K .* (1 .- delta)
    p2 = merge(p, (K = K2,))
    prob2 = ODEProblem(allo_ode!, pre, (tpert, tspan[2]), p2)
    sol2  = solve(prob2, solver; callback=cb, abstol=1e-8, reltol=1e-8)
    post  = sol2.u[end]
    after = count(>(EXTINCTION_THRESHOLD), post)/S

    rtn = fill(NaN, S)
    for i in 1:S
        tgt = post[i]
        for (t, st) in zip(sol2.t, sol2.u)
            if abs(st[i] - tgt)/(abs(tgt)+1e-8) < 0.10
                rtn[i] = t - tpert
                break
            end
        end
    end

    return rtn, before, after, post
end

#############################################################################
#############################################################################
################## SIMULATE PULSE PERTURBATIONS #############################
#############################################################################
#############################################################################
"""
    simulate_pulse_perturbation(
        u0, p, tspan, t_pulse, delt;
        solver=Tsit5(),
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
    solver=Tsit5(),
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
        fig = Figure(resolution=(1200,600))
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
# u[k] with entries B.^2, normalize each to unit norm, then compute
# alpha_k = |exp(t*J)*u[k]|.  Returns
# R_med(t) = -(1/t) * log(median(alpha_k)).
# """
function median_return_rate(J, B; t=1.0, n=1000, rng=Random.GLOBAL_RNG)
    # if J has any non-finite entry, skip
    if any(!isfinite, J)
        return NaN
    end

    E = exp(t*J)
    S = length(B)
    alpha = Vector{Float64}(undef, n)

    for k in 1:n
        u = randn(rng, S) .* (B .^ 2)
        u ./= norm(u)
        alpha[k] = norm(E * u)
    end

    return -log(median(alpha)) / t
end
