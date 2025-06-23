# ─────────────────────────────────────────────────────────────────────────────
# 1) GLV RHS & Jacobian
# ─────────────────────────────────────────────────────────────────────────────
function glv!(dB,B,p,t)
    r, K, A = p.r, p.K, p.A
    S = length(B)
    @inbounds for i in 1:S
        dB[i] = r[i]*B[i]*(1 - B[i]/K[i])
        @inbounds for j in 1:S
            dB[i] += A[i,j]*B[i]*B[j]
        end
    end
end

function jacobian_glv(B,p)
    r, K, A = p.r, p.K, p.A; S = length(B)
    J = zeros(S,S)
    @inbounds for i in 1:S
        J[i,i] = r[i]*(1 - 2B[i]/K[i])
        @inbounds for j in 1:S
            if i != j
                J[i,i] += A[i,j]*B[j]
                J[i,j]  = A[i,j]*B[i]
            end
        end
    end
    return J
end

# ─────────────────────────────────────────────────────────────────────────────
# 2) Attempt to build one stable community
# ─────────────────────────────────────────────────────────────────────────────
function build_one(Rb,Ri,Rt;
    conn_bi, conn_it, conn_bt,
    IS, r_scale, K_scale,
    surv_frac=0.8, T=200.0, max_tries=50)

    S = Rb+Ri+Rt
    for _ in 1:max_tries
        # intrinsic rates & capacities
        r = vcat(fill(1.0* r_scale, Rb), zeros(Ri+Rt))
        K = fill(10.0* K_scale, S)

        # interaction matrix with strict three‐layer links
        A = zeros(S,S)
        ## basal→inter
        for b in 1:Rb, i in 1:Ri
            if rand()<conn_bi
                A[b,   Rb+i] =  randn()*IS
                A[Rb+i,b   ] = -randn()*IS
            end
        end
        ## inter→top
        for i in 1:Ri, t in 1:Rt
            if rand()<conn_it
                A[Rb+i, Rb+Ri+t] =  randn()*IS
                A[Rb+Ri+t, Rb+i] = -randn()*IS
            end
        end
        ## basal→top (omnivory)
        # for b in 1:Rb, t in 1:Rt
        #     if rand()<conn_bt
        #         A[b,   Rb+Ri+t] =  randn()*IS*0.5
        #         A[Rb+Ri+t, b   ] = -randn()*IS*0.5
        #     end
        # end

        # simulate GLV
        p = (r=r, K=K, A=A)
        prob = ODEProblem(glv!, fill(0.1,S), (0.0,T), p)
        sol = solve(prob, abstol=1e-8, reltol=1e-8)
        Bstar = sol.u[end]
        survivors = count(x-> x>1e-3, Bstar)
        if survivors >= surv_frac*S && sol.t[end]==T
            # local stability
            Jstar = jacobian_glv(Bstar, p)
            if all(real(eigvals(Jstar)) .< 0)
                return (r=r, K=K, A=A, Bstar=Bstar, Jstar=Jstar, sol=sol)
            end
        end
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# 3) Sweep over everything & record
# ─────────────────────────────────────────────────────────────────────────────
function explore_communities(; S_total=30,
    conn_bi_vals=0.1:0.1:0.5,
    conn_it_vals=0.1:0.1:0.5,
    conn_bt_vals=0.05:0.05:0.05,
    IS_vals=0.1:0.2:1.0,
    r_scales=[0.5,1.0,2.0],
    K_scales=[0.5,1.0,2.0],
    surv_frac=0.8,
    T=200.0,
    max_tries=30)

    rows = NamedTuple[]
    combo_id = 0
    locki = ReentrantLock()

    @threads for (Rb, Ri) in (5:5:15,5:5:10)
        Rt = S_total - Rb - Ri
        Rt<1 && continue
        for conn_bi in conn_bi_vals, conn_it in conn_it_vals, conn_bt in conn_bt_vals,
            IS in IS_vals, r_sc in r_scales, K_sc in K_scales

            combo_id += 1
            got = build_one(
                Rb,Ri,Rt;
                conn_bi=conn_bi, conn_it=conn_it, conn_bt=conn_bt,
                IS=IS, r_scale=r_sc, K_scale=K_sc,
                surv_frac=surv_frac, T=T, max_tries=max_tries
            )

            if got !== nothing
                lock(locki) do
                    push!(rows, (
                    combo_id=combo_id,
                    Rb=Rb, Ri=Ri, Rt=Rt,
                    conn_bi=conn_bi,
                    conn_it=conn_it,
                    conn_bt=conn_bt,
                    IS=IS,
                    r_scale=r_sc,
                    K_scale=K_sc,
                    Bstar=got.Bstar,
                    Jstar=got.Jstar
                    ))
                end   
            end
        end
    end

    return DataFrame(rows)
end

# ─────────────────────────────────────────────────────────────────────────────
# 4) Run it and save
# ─────────────────────────────────────────────────────────────────────────────
df = explore_communities(
  S_total=30,
  conn_bi_vals=0.1:0.2:0.5,
  conn_it_vals=0.1:0.2:0.5,
  conn_bt_vals=0.05:0.05:0.05,
  IS_vals=0.1:0.4:1.0,
  r_scales=[1.5, 1.0],
  K_scales=[0.5,1.0],#2.0],
  surv_frac=0.8,
  T=200.0,
  max_tries=1
)

println("Found $(nrow(df)) feasible communities.")
CSV.write("communities.csv", df)
