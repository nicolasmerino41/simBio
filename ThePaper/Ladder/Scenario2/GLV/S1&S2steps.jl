# — analytic resilience/reactivity (if not already in scope) —
measure_resilience(J) = maximum(real(eigvals(J)))
measure_reactivity(J) = maximum(real(eigvals((J + J')/2)))

# — simulated persistence —
"""
    measure_persistence(r, K, A, Bstar; δ=0.1, Tpost=100.0, tol=1e-3)

Starting from equilibrium Bstar, apply a downward perturbation
Bₚ = max.(Bstar .- δ .* K, 0.0), then simulate GLV dB/dt = ... from Bₚ
for time Tpost.  Returns:

    fraction_survived = (# species with B(Tpost) > tol) / (# species with Bstar > tol)

Use δ in (0,1), Tpost sufficiently long (e.g. 50–200), and tol your extinction threshold.
"""
function measure_persistence(r, K, A, Bstar; δ=0.1, Tpost=100.0, tol=1e-3)
    S = length(Bstar)
    # 1) how many were alive just before?
    pre = count(x -> x > tol, Bstar)
    if pre == 0
        return 0.0
    end

    # 2) build perturbed initial condition
    Bp = max.(Bstar .- δ .* K, 0.0)
    # simulate GLV from Bp
    p = (r=r, K=K, A=A)
    prob = ODEProblem(glv!, Bp, (0.0, Tpost), p)
    sol = solve(prob; abstol=1e-8, reltol=1e-8)
    Bend = sol.u[end]

    # 3) how many survive after?
    post = count(x -> x > tol, Bend)

    return post / pre
end

# assume df is your DataFrame of size N
resilience   = Float64[]
reactivity   = Float64[]
persistence  = Float64[]
for row in eachrow(df)
    # analytic
    push!(resilience,   measure_resilience(row.Jstar))
    push!(reactivity,   measure_reactivity(row.Jstar))

    # simulated
    push!(persistence,
        measure_persistence(
            row.r, row.K, row.A_mat, row.Bstar;
            δ    = 0.2,         # e.g. 20% perturbation
            Tpost=100.0,
            tol  = 1e-3
        ))
end

df.resilience  = resilience
df.reactivity  = reactivity
df.persistence = persistence

###############################################################################
###############################################################################
###############################################################################
# — helper to reshuffle **within** each trophic layer (not between them) —
function reshuffle_within_groups(A::AbstractMatrix, Rb::Int, Ri::Int, Rt::Int)
    S = size(A,1)
    levels = vcat(fill(1,Rb), fill(2,Ri), fill(3,Rt))
    B = copy(A)
    for lvl in (1,2,3)
        # collect all off‐diagonal entries (i≠j) within this level
        idx = findall(x->x==lvl, levels)
        pairs = [(i,j) for i in idx, j in idx if i!=j]
        vals = [A[i,j] for (i,j) in pairs] .* -1.0
        p   = randperm(length(vals))
        for k in 1:length(pairs)
            i,j = pairs[k]
            B[i,j] = vals[p[k]]
        end
    end
    return B
end

"""
    reshuffle_between_levels(A, levels; ordered_levels)

Permute **only** the links from each level ℓ₁ → ℓ₂, for all (ℓ₁,ℓ₂) in `ordered_levels`,
while keeping the values themselves fixed.  `levels` is a vector of length S giving
the level‐index of each species (e.g. 1=basal, 2=inter, 3=top).  
`ordered_levels` is a collection of pairs, e.g. `[(1,2),(2,3)]`.
"""
function reshuffle_between_levels(A::AbstractMatrix, levels::AbstractVector{<:Integer};
                                   ordered_levels = [(1,2), (2,3)])
    S = size(A,1)
    B = copy(A)
    for (ℓ1,ℓ2) in ordered_levels
        # collect all (i,j) with i in ℓ1 and j in ℓ2
        idx1 = findall(x->x==ℓ1, levels)
        idx2 = findall(x->x==ℓ2, levels)
        pairs = [(i,j) for i in idx1 for j in idx2]
        vals   = [ A[i,j] for (i,j) in pairs ]
        perm   = randperm(length(vals))
        for (k,(i,j)) in enumerate(pairs)
            B[i,j] = vals[perm[k]]
        end
    end
    return B
end

function add_S1!(df::DataFrame)
    n = nrow(df)

    # 1) Pre‐allocate the new columns
    df.persistence_S1 = Vector{Float64}(undef, n)
    df.resilience_S1  = Vector{Float64}(undef, n)
    df.reactivity_S1  = Vector{Float64}(undef, n)
    df.A_mat_S1       = Vector{Matrix{Float64}}(undef, n)
    df.Jstar_S1       = Vector{Matrix{Float64}}(undef, n)
    cb = build_callbacks(df.Rb[1]+df.Ri[1]+df.Rt[1], EXTINCTION_THRESHOLD)
    # 2) Loop by row‐index
    for i in 1:n
        Rb, Ri, Rt = df.Rb[i], df.Ri[i], df.Rt[i]
        A0 = df.A_mat[i]
        Bstar = df.Bstar[i]
        r0 = df.r[i]
        K0 = df.K[i]

        # a) reshuffle within-groups
        # A1 = reshuffle_within_groups(A0, Rb, Ri, Rt)
        levels = vcat(fill(1,Rb), fill(2,Ri), fill(3,Rt))
        A1 = reshuffle_between_levels(A0, levels; ordered_levels=[(1,2),(2,3),(1,3)])

        # b) recompute Jacobian at the same equilibrium
        p1 = (r = r0, K = K0, A = A1)
        
        prob = ODEProblem(glv!, Bstar, (0.0,200), p1)
        sol = solve(prob; callback=cb, abstol=1e-8, reltol=1e-8)
        Bstar = sol.u[end]
        J1 = jacobian_glv(Bstar, p1)
        # c) measure metrics
        pers1 = measure_persistence(r0, K0, A1, Bstar)
        res1  = measure_resilience(J1)
        rea1  = measure_reactivity(J1)

        # d) write into the dataframe
        df.persistence_S1[i] = pers1
        df.resilience_S1[i]  = res1
        df.reactivity_S1[i]  = rea1
        df.A_mat_S1[i]       = A1
        df.Jstar_S1[i]       = J1
    end

    return df
end

# suppose `df` already has all your original columns:
df2 = add_S1!(df)
df2.Δres = df2.resilience - df2.resilience_S1
df2.Δrea = df2.reactivity - df2.reactivity_S1
df2.Δper = df2.persistence - df2.persistence_S1

"""
    reshuffle_between_groups(A0, r0, Rb, Ri, Rt, change)

Move `change` species out of the basal layer into the intermediate/top layers:
  newRb = Rb - change
  newRi = Ri + change÷2
  newRt = Rt + change÷2

Returns:
  A1   = A₀ with all off‐diagonal links reshuffled *within* each of the new three layers,
  r1   = new growth‐rate vector (first newRb entries copied from r₀, rest zero),
  newRb, newRi, newRt
"""
function reshuffle_between_groups(
    A0::AbstractMatrix{<:Real},
    r0::AbstractVector{<:Real},
    Rb::Int, Ri::Int, Rt::Int,
    change::Int
)
    # compute new layer sizes
    newRb = Rb - change
    newRi = Ri + change ÷ 2
    newRt = Rt + change ÷ 2
    S = size(A0,1)
    @assert newRb ≥ 0 && newRi ≥ 0 && newRt ≥ 0 && newRb+newRi+newRt==S
        "Invalid change: layers must still sum to S"

    # rebuild r:
    r1 = zeros(eltype(r0), S)
    # carry over the first newRb basal rates
    r1[1:newRb] = r0[1:newRb]

    # reshuffle A within each new layer
    levels = vcat(fill(1,newRb), fill(2,newRi), fill(3,newRt))
    A1 = zeros(eltype(A0), S, S)
    # copy diagonal from A0 (if any) or zero
    for i in 1:S
        A1[i,i] = A0[i,i]
    end
    # for each layer, reshuffle its off‐diagonal entries
    for lvl in (1,2,3)
        idx = findall(x->x==lvl, levels)
        # collect 
        pairs = [(i,j) for i in idx for j in idx if i≠j]
        vals = [A0[i,j] for (i,j) in pairs]
        perm = randperm(length(vals))
        for k in 1:length(pairs)
            i,j = pairs[k]
            A1[i,j] = vals[perm[k]]
        end
    end

    return A1, r1, newRb, newRi, newRt
end

function add_S2!(df::DataFrame)
    n = nrow(df)

    # 1) Pre‐allocate the new columns
    df.persistence_S2 = Vector{Float64}(undef, n)
    df.resilience_S2  = Vector{Float64}(undef, n)
    df.reactivity_S2  = Vector{Float64}(undef, n)
    df.A_mat_S2       = Vector{Matrix{Float64}}(undef, n)
    df.Jstar_S2       = Vector{Matrix{Float64}}(undef, n)
    cb = build_callbacks(df.Rb[1]+df.Ri[1]+df.Rt[1], EXTINCTION_THRESHOLD)
    # 2) Loop by row‐index
    for i in 1:n
        Rb, Ri, Rt = df.Rb[i], df.Ri[i], df.Rt[i]
        A0 = df.A_mat[i]
        Bstar = df.Bstar[i]
        r0 = df.r[i]
        K0 = df.K[i]

        # a) reshuffle within-groups
        A1, r1, Rb, Ri, Rt = reshuffle_between_groups(A0, r0, Rb, Ri, Rt, 4)

        ########################################################################
        # S = Rb+Ri+Rt
        # for _ in 1:max_tries
        #     # intrinsic rates & capacities
        #     r = r0
        #     r[Rb+1:Rb+4] .= 0.0
        #     K = K0
        #     Rb = Rb - 4
        #     Ri = Ri + 2
        #     Rt = Rt + 2
        #     # r = vcat(fill(1.0* r_scale, Rb), zeros(Ri+Rt))
        #     # K = fill(10.0* K_scale, S)

        #     # interaction matrix with strict three‐layer links
        #     A = zeros(S,S)
        #     ## basal→inter
        #     for b in 1:Rb, i in 1:Ri
        #         if rand()<conn_bi
        #             A[b,   Rb+i] =  rand()*df.IS
        #             A[Rb+i,b   ] = -rand()*df.IS
        #         end
        #     end
        #     ## inter→top
        #     for i in 1:Ri, t in 1:Rt
        #         if rand()<conn_it
        #             A[Rb+i, Rb+Ri+t] =  rand()*df.IS
        #             A[Rb+Ri+t, Rb+i] = -rand()*df.IS
        #         end
        #     end
        #     ## basal→top (omnivory)
        #     # for b in 1:Rb, t in 1:Rt
        #     #     if rand()<conn_bt
        #     #         A[b,   Rb+Ri+t] =  rand()*IS*0.5
        #     #         A[Rb+Ri+t, b   ] = -rand()*IS*0.5
        #     #     end
        #     # end

        #     # simulate GLV
        #     p = (r=r, K=K, A=A)
        #     prob = ODEProblem(glv!, fill(0.1,S), (0.0,T), p)
        #     sol = solve(prob; callback=cb, abstol=1e-8, reltol=1e-8)
        #     Bstar = sol.u[end]
        #     survivors = count(x-> x>1e-3, Bstar)
        #     if survivors >= surv_frac*S && sol.t[end]==T
        #         # local stability
        #         Jstar = jacobian_glv(Bstar, p)
        #         if all(real(eigvals(Jstar)) .< 0)
        #             return (r=r, K=K, A=A, Bstar=Bstar, Jstar=Jstar, sol=sol)
        #         end
        #     end
        # end
        ##########################################################################

        # b) recompute Jacobian at the same equilibrium
        p1 = (r = r1, K = K0, A = A1)
        
        prob = ODEProblem(glv!, Bstar, (0.0,200), p1)
        sol = solve(prob; callback=cb, abstol=1e-8, reltol=1e-8)
        Bstar = sol.u[end]
        J1 = jacobian_glv(Bstar, p1)
        # c) measure metrics
        pers1 = measure_persistence(r0, K0, A1, Bstar)
        res1  = measure_resilience(J1)
        rea1  = measure_reactivity(J1)

        # d) write into the dataframe
        df.persistence_S2[i] = pers1
        df.resilience_S2[i]  = res1
        df.reactivity_S2[i]  = rea1
        df.A_mat_S2[i]       = A1
        df.Jstar_S2[i]       = J1
    end

    return df
end

df3 = add_S2!(df2)

Δres2 = df3.resilience - df3.resilience_S2
Δrea2 = df3.reactivity - df3.reactivity_S2
Δper2 = df3.persistence - df3.persistence_S2

df4 = filter(
    row -> row.resilience_S1 .< 0 &&
    #  row.resilience_S2 .< 0 &&
      abs(row.resilience_S1) .< 5.0,
    #    abs(row.resilience_S2) .< 5.0,
    df3
)