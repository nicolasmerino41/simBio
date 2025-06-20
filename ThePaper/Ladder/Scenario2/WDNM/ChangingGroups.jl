# ─────────────────────────────────────────────────────────────────────────────
# 1) GLV RHS & Jacobian
# ─────────────────────────────────────────────────────────────────────────────
function glv!(dB,B,p,t)
    r, K, A = p.r, p.K, p.A
    S = length(B)
    @inbounds for i in 1:S
        dB[i] = r[i] * B[i] * (1 - B[i]/K[i])
        for j in 1:S
            if i != j
                dB[i] += A[i,j] * B[i] * B[j]
            end
        end
    end
end

function jacobian_glv(B,p)
    r, K, A = p.r, p.K, p.A
    S = length(B)
    J = zeros(S,S)
    @inbounds for i in 1:S
        # diagonal
        J[i,i] = r[i]*(1 - 2B[i]/K[i])
        for j in 1:S
            if i != j
                J[i,i] += A[i,j]*B[j]
                J[i,j]  = A[i,j]*B[i]
            end
        end
    end
    return J
end

measure_reactivity(J) = maximum(eigvals((J+J')/2))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Build & test feasibility
# ─────────────────────────────────────────────────────────────────────────────
"""
    build_and_test(Rb,Ri,Rt; conn_bi, conn_it, conn_bt, max_tries)

Attempts up to max_tries random draws of an interaction web connecting
basal→inter (p=conn_bi), inter→top (conn_it), basal→top (conn_bt).
Assigns positive r to basal, zero to others, K=10, A sampled ≈ Uniform(–1,1).
Keeps only those where all species survive a GLV run.
Returns (r,K,A,Bstar) on success, or nothing on failure.
"""
function build_and_test(Rb, Ri, Rt; conn_bi=0.3, conn_it=0.3, conn_bt=0.1, max_tries=50)
    S = Rb+Ri+Rt
    for _ in 1:max_tries
        # growth & capacity
        r = vcat(fill(1.0,Rb), zeros(Ri+Rt))
        K = fill(10.0, S)
        # interaction matrix
        A = zeros(S,S)
        # basal→inter
        for i in 1:Ri, j in 1:Rb
            if rand()<conn_bi
                A[j,Rb+i] =  randn()*0.5
                A[Rb+i,j] = -randn()*0.5
            end
        end
        # inter→top
        for i in 1:Rt, j in 1:Ri
            if rand()<conn_it
                A[Rb+j, Rb+Ri+i] =  randn()*0.5
                A[Rb+Ri+i, Rb+j] = -randn()*0.5
            end
        end
        # basal→top (omnivory)
        for i in 1:Rt, j in 1:Rb
            if rand()<conn_bt
                A[j, Rb+Ri+i] =  randn()*0.2
                A[Rb+Ri+i, j] = -randn()*0.2
            end
        end

        # simulate
        p = (r=r, K=K, A=A)
        prob = ODEProblem(glv!, fill(0.1,S), (0.0,200.0), p)
        sol = solve(prob, abstol=1e-6, reltol=1e-6)
        Bstar = sol.u[end]
        if all(x->x>1e-3, Bstar)
            return r, K, A, Bstar
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 3) Reshaping schemes
# ─────────────────────────────────────────────────────────────────────────────
function reshuffle_matrix(A, levels, scheme::Symbol)
    S = size(A,1)
    B = copy(A)
    idxs = Tuple{Int,Int}[]
    for i in 1:S, j in 1:S
        if i==j
            continue
        end
        same = levels[i]==levels[j]
        if scheme==:full ||
           (scheme==:inter && !same) ||
           (scheme==:intra && same)
            push!(idxs,(i,j))
        end
    end
    vals = shuffle(A[i] for i in idxs)
    for (k,(i,j)) in enumerate(idxs)
        B[i,j] = vals[k]
    end
    return B
end

# ─────────────────────────────────────────────────────────────────────────────
# 4) Experiment
# ─────────────────────────────────────────────────────────────────────────────
function experiment(Rb,Ri,Rt; nref=50, reps=20)
    results = Dict(s=>Float64[] for s in (:full,:inter,:intra))
    levels = vcat(fill(1,Rb), fill(2,Ri), fill(3,Rt))

    for _ in 1:nref
        got = build_and_test(Rb,Ri,Rt)
        if isnothing(got)
            continue
        end
        r,K,A,B0 = got
        p0 = (r=r,K=K,A=A)
        J0 = jacobian_glv(B0,p0)
        react0 = measure_reactivity(J0)
        pers0  = 1.0  # by construction

        for scheme in (:full,:inter,:intra)
            Δr, Δp = Float64[], Float64[]
            for _ in 1:reps
                A2 = reshuffle_matrix(A, levels, scheme)
                p2 = (r=r,K=K,A=A2)
                sol2 = solve(ODEProblem(glv!,B0,(0.0,200.0),p2),
                             abstol=1e-6, reltol=1e-6)
                B2 = sol2.u[end]
                push!(Δp, abs(count(>=(1e-3),B2)/length(B2) - pers0))

                J2 = jacobian_glv(B2,p2)
                push!(Δr, abs(measure_reactivity(J2) - react0))
            end
            # trimmed‐mean
            sort!(Δp); trimp = mean(Δp[ceil(Int,0.05*reps)+1 : floor(Int,0.95*reps)])
            sort!(Δr); trimr = mean(Δr[ceil(Int,0.05*reps)+1 : floor(Int,0.95*reps)])
            push!(results[scheme], (trim_persistence=trimp, trim_reactivity=trimr))
        end
    end

    return results
end

# ─────────────────────────────────────────────────────────────────────────────
# 5) Run & summarize
# ─────────────────────────────────────────────────────────────────────────────
res = experiment(20,15,10; nref=100, reps=30)
for s in (:full,:inter,:intra)
    arr = res[s]
    println(s,": Δ persistence = ",
      mean(getfield.(arr, :trim_persistence)),
      "   Δ reactivity = ",
      mean(getfield.(arr, :trim_reactivity)))
end
