using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")   # defines make_A, sample_Xi_K_from_A, calibrate_from_K_xi, compute_jacobian, trophic_ode!, etc.

# ── helper to build J and return its diagonal D, “mixing” M, and full J = D*M
function build_real_J(R::Int, C::Int; conn, IS, scen, eps_scale=1.0)
    S = R + C
    A = make_A(zeros(S,S), R, conn, scen; IS=IS)
    ε = clamp.(randn(S,S) * eps_scale, 0, 1)
    K_res, xi_cons = sample_Xi_K_from_A(A, ε, R; scale_factor=1.0)
    full_K = vcat(K_res, xi_cons)      # now length S
    B0     = calibrate_from_K_xi(xi_cons, K_res, ε, A)
    D, M   = compute_jacobian(B0, (R, C,
                                    fill(1.0,C), xi_cons,
                                    fill(1.0,R), K_res, ε, A))
    return D, M, D*M, A, ε, B0, full_K
end

# ── the eight Δ-metric functions
function measure_resilience(J::Matrix)
    return -maximum(real(eigvals(J)))
end

function measure_reactivity(J::Matrix)
    Js = (J + J')/2
    maximum(real(eigvals(Js)))
end

function analytical_Rmed(J; t=0.01)
    E = exp(t*J)
    M = E*E'
    rates = -((diag(M) .- 1.0) ./ (2t))
    mean(rates)
end

pulse_return_time(J) = 1/measure_resilience(J)

function measure_persistence(J::Matrix)
    ev = real.(eigvals(J))
    count(<(0), ev) / length(ev)
end

function measure_collectivity(d::Vector, S; σ=0.1)
    G = randn(S,S)*σ
    for i in 1:S; G[i,i]=0; end
    maximum(abs, eigvals(G))
end

function nonnormality(J::AbstractMatrix)
    C = J*J' - J'*J
    norm(C)
end

function compute_SL(A::Matrix, K::Vector)
    B = (I - A) \ K
    B ./ K
end

# ── compare_diag_effect
function compare_diag_effect(R=50, C=20; scen=:ER, σ=1.0, reps=200)
    Random.seed!(1234)

    # 1) Build reference once
    D0, M0, J0, A0, ε0, B0, K0 = build_real_J(R, C; conn=0.2, IS=σ, scen=scen)

    metrics = (:Res, :Rea, :Rmed) # :PulseRT, :Pers, :Coll, :Nonnorm, :SL)
    Δ_full = Dict(m=>Float64[] for m in metrics)
    Δ_diag = deepcopy(Δ_full)

    for _ in 1:reps
        # full rebuild
        D1, M1, J1, A1, ε1, B1, K1 = build_real_J(R, C; conn=0.2, IS=σ, scen=scen)

        # diag-fixed: rebuild then re-impose D0
        _, M2, _, A2, ε2, B2, _ = build_real_J(R, C; conn=0.2, IS=σ, scen=scen)
        J2 = D0 * M2

        # full vs ref
        push!(Δ_full[:Res],     abs(measure_resilience(J1)    - measure_resilience(J0)))
        push!(Δ_full[:Rea],     abs(measure_reactivity(J1)     - measure_reactivity(J0)))
        push!(Δ_full[:Rmed],    abs(analytical_Rmed(J1)        - analytical_Rmed(J0)))
        # push!(Δ_full[:PulseRT], abs(pulse_return_time(J1)      - pulse_return_time(J0)))
        # push!(Δ_full[:Pers],    abs(measure_persistence(J1)    - measure_persistence(J0)))
        # push!(Δ_full[:Coll],    abs(measure_collectivity(-diag(J1),R+C;σ=σ) -
        #                           measure_collectivity(-diag(J0),R+C;σ=σ)))
        # push!(Δ_full[:Nonnorm], abs(nonnormality(J1)           - nonnormality(J0)))

        # # SL full:
        # sl_ref0  = compute_SL(A0, K0)
        # sl_full0 = compute_SL(A1, K1)
        # push!(Δ_full[:SL], mean(abs.(sl_full0 .- sl_ref0)))

        # diag-fixed vs ref
        push!(Δ_diag[:Res],     abs(measure_resilience(J2)    - measure_resilience(J0)))
        push!(Δ_diag[:Rea],     abs(measure_reactivity(J2)     - measure_reactivity(J0)))
        push!(Δ_diag[:Rmed],    abs(analytical_Rmed(J2)        - analytical_Rmed(J0)))
        # push!(Δ_diag[:PulseRT], abs(pulse_return_time(J2)      - pulse_return_time(J0)))
        # push!(Δ_diag[:Pers],    abs(measure_persistence(J2)    - measure_persistence(J0)))
        # # collectivity identical by construction:
        # push!(Δ_diag[:Coll], Δ_full[:Coll][end])
        # push!(Δ_diag[:Nonnorm], abs(nonnormality(J2)         - nonnormality(J0)))

        # # SL diag-fixed (use K0):
        # sl_diag0 = compute_SL(A2, K0)
        # push!(Δ_diag[:SL], mean(abs.(sl_diag0 .- sl_ref0)))
    end

    # compute means
    mean_full = Dict(k=>mean(v) for (k,v) in Δ_full)
    mean_diag = Dict(k=>mean(v) for (k,v) in Δ_diag)
    # metrics = metrics[1:1]
    fig = Figure(; size=(1000,450))
    for (i, m) in enumerate(metrics)
        row, col = divrem(i-1,4)  # row = 0 or 1, col = 0…3
        ax = Axis(fig[row+1, col+1])
        xs = [1,2]
        ys = [mean_full[m], mean_diag[m]]
        barplot!(ax, xs, ys)
        # ax.title = titles[i]
        ax.xticks = (xs, ["full","diag-fix"])
        # optionally fix y-limits across all panels:
        # ylims!(ax, 0, maximum(ys)*1.1)
    end

    display(fig)

    return mean_full, mean_diag
end

# ── run it at several σ’s
for σ in (0.1, 1.0, 5.0)
    mf, md = compare_diag_effect(50,20; scen=:ER, σ=σ, reps=500)
    println("σ=$σ full vs diag-fixed:")
    @show mf; @show md
end
