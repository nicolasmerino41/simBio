# ——— Core “random‐matrix” metrics ————————————————————————————————————————————

# Compute resilience = –max Re(eig(J))
function measure_resilience(J::Matrix{Float64})
    vals = eigvals(J)
    return -maximum(real(vals))
end

# Compute reactivity = max eigenvalue of the symmetrized Jacobian
function measure_reactivity(J::Matrix{Float64})
    Jsym = (J + J')/2
    return maximum(real(eigvals(Jsym)))
end

# Analytical median return rate
function analytical_Rmed(J; t=0.01)
    E = exp(t*J)
    M = E*E'
    rates = diag(M) .- 1 .|> x-> -x/(2t)
    return mean(rates)
end

# Pulse‐return‐time approx = 1/resilience
pulse_return_time(J) = 1/measure_resilience(J)

# Persistence ≈ fraction of eigenvalues with negative real part
measure_persistence(J) = count(λ->real(λ)<0, eigvals(J))/length(eigvals(J))

# Collectivity = spectral radius of G part (zero‐mean)
function measure_collectivity(d,S;σ=0.1)
    G = randn(S,S)*σ
    for i in 1:S
        G[i,i] = 0
    end
    return maximum(abs, eigvals(G))
end

# Henrici non‐normality
function nonnormality(J::AbstractMatrix)
    C = J*J' - J'*J
    return norm(C)
end

# ——— Build a “real” J from your full pipeline —————————————————————————————

"""
    build_real_J(R, C; conn, IS, scen, eps_scale)

Construct a full Jacobian J = D*M from your trophic framework:

- make_A → ε
- sample Xi,K → equilibrium B0
- compute_jacobian → (D,M)
Returns: (J, A, ε, B0, K_full)
"""
function build_real_J(R::Int, C::Int; conn, IS, scen, eps_scale=1.0)
    S = R + C

    # 1) network A + ε
    A = make_A(zeros(S,S), R, conn, scen; IS=IS)
    ε = clamp.(randn(S,S)*eps_scale, 0, 1)

    # 2) sample one feasible Xi,K and equilibrium
    K_res, xi_cons = sample_Xi_K_from_A(A, ε, R; scale_factor=1.0)
    B0 = calibrate_from_K_xi(xi_cons, K_res, ε, A)

    # 3) compute Jacobian at B0 with unit rates (we only care about structure)
    #    here we pass dummy m_cons = ones(C), r_res = ones(R)
    p = (R, C, fill(1.0,C), xi_cons, fill(1.0,R), K_res, ε, A)
    D, M = compute_jacobian(B0, p)

    # full carrying‐capacity vector:
    K_full = vcat(K_res, xi_cons)

    return D*M, A, ε, B0, K_full
end

# ——— Compute all Δ‐metrics on “real” J’s ——————————————————————————————————

"""
    Δmetrics_real(d, S; σs, reps, R, C, scen)

For each σ in σs, build two independent J’s (with same IS=σ),
compute all eight Δ‐metrics, and return
a tuple of eight vectors.
"""
# compute a trimmed mean by discarding the lowest p and highest p fraction
function trimmed_mean(v::Vector{<:Real}; p::Float64 = 0.05)
    @assert 0 ≤ p < 0.5 "trim fraction must be in [0,0.5)"
    n = length(v)
    if n == 0
        return NaN
    elseif n < 2
        return mean(v)
    end
    lo = floor(Int,  p * n)
    hi = ceil(Int,  (1-p) * n)
    vs = sort(v)
    return mean(vs[lo+1 : hi])
end

function Δmetrics_real(d::Vector{Float64}, S::Int; σs, reps, R, C, scen)
    Δρ    = Float64[]; Δrea  = Float64[]; ΔRmed = Float64[]
    Δprt  = Float64[]; Δpers = Float64[]; Δphi  = Float64[]
    Δnn   = Float64[]; Δsl   = Float64[]

    for σ in σs
        bufρ=Float64[]; bufrea=Float64[]; bufRmed=Float64[]
        bufprt=Float64[]; bufpers=Float64[]; bufphi=Float64[]
        bufnn=Float64[]; buffsl=Float64[]

        for _ in 1:reps
            J0,A0,ε0,B0,K0 = build_real_J(R, C; conn=0.2, IS=σ, scen=scen, eps_scale=1.0)
            J1,A1,ε1,B1,K1 = build_real_J(R, C; conn=0.2, IS=σ, scen=scen, eps_scale=1.0)

            push!(bufρ,    abs(measure_resilience(J1)   - measure_resilience(J0)))
            push!(bufrea,  abs(measure_reactivity(J1)    - measure_reactivity(J0)))
            push!(bufRmed, abs(analytical_Rmed(J1)       - analytical_Rmed(J0)))
            push!(bufprt,  abs(pulse_return_time(J1)     - pulse_return_time(J0)))
            push!(bufpers, abs(measure_persistence(J1)   - measure_persistence(J0)))
            push!(bufphi,  abs(measure_collectivity(d,S;σ=σ) - measure_collectivity(d,S;σ=σ)))
            push!(bufnn,   abs(nonnormality(J1)         - nonnormality(J0)))

            sl0 = compute_SL(A0, K0)
            sl1 = compute_SL(A1, K1)
            push!(buffsl,  mean(abs.(sl1 .- sl0)))
        end

        # use trimmed_mean instead of mean
        push!(Δρ,    trimmed_mean(bufρ))
        push!(Δrea,  trimmed_mean(bufrea))
        push!(ΔRmed, trimmed_mean(bufRmed))
        push!(Δprt,  trimmed_mean(bufprt))
        push!(Δpers, trimmed_mean(bufpers))
        push!(Δphi,  trimmed_mean(bufphi))
        push!(Δnn,   trimmed_mean(bufnn))
        push!(Δsl,   trimmed_mean(buffsl))
    end

    return (Δρ, Δrea, ΔRmed, Δprt, Δpers, Δphi, Δnn, Δsl)
end

function test_real(R=50, C=20; scen=:ER, reps=200)
    Random.seed!(123)

    Jref, _, _, _, _ = build_real_J(R,C; conn=0.2, IS=0.1, scen=scen)
    d     = abs.(-diag(Jref))
    min_d = minimum(d)

    # Sweep from 10⁻³ all the way up to 10²
    ratios     = 10 .^ range(-3, 2; length=25)

    σs_IS = range(0.01, 10.0; length=25)

    σs_sigmd  = ratios .* min_d

    out_IS   = Δmetrics_real(d, R+C; σs=σs_IS,    reps=reps, R=R, C=C, scen=scen)
    out_sig  = Δmetrics_real(d, R+C; σs=σs_sigmd, reps=reps, R=R, C=C, scen=scen)

    fig = Figure(; size=(1400,600))
    titles   = ["Δ Res","Δ Rea","Δ Rmed","Δ PulseRT","Δ Pers","Δ Coll","Δ Nonnorm","Δ SL"]
    datasets = (out_IS, out_sig)
    xaxes    = (σs_IS, σs_sigmd)
    xlabels  = ("IS (σ)", "σ/min(d)")

    for (k,(Δs,xs,xlab)) in enumerate(zip(datasets, xaxes, xlabels))
        for i in 1:8
            col, row = divrem(i-1,4)
            ax = Axis(fig[row+1, col*2 + k])
            lines!(ax, xs, Δs[i]; linewidth=2)
            scatter!(ax, xs, Δs[i]; markersize=6)
            ax.xlabel = xlab
            ax.ylabel = titles[i]
            if k == 2
                ax.xscale = log10
            end
        end
    end

    display(fig)
end

# run it!
test_real()

