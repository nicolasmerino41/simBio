# ── helper to build J and return its diagonal D, "mixing" M, and full J = D*M
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
    rates = -((diag(M) .- 1.0) ./ (2 * t))
    @assert all(isfinite, rates) "Found NaN or Inf in Rmed computation"
    return mean(rates)
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

function compare_diag_effect(
    R=50, C=20;
    scen=:ER, σ=1.0,
    reps=200, MAX_TRIES=20,
)
    Random.seed!(1234)

    # 1) Build a *stable* reference J0 (retry up to MAX_TRIES)
    D0 = nothing; M0 = nothing; J0 = nothing
    for attempt in 1:MAX_TRIES
        d0, m0, j0, _, _, _, _ = build_real_J(R, C; conn=0.5, IS=σ, scen=scen)
        if true
            D0, M0, J0 = d0, m0, j0
            break
        end
    end
    if J0 === nothing
        error("Could not generate a stable reference Jacobian J0 after $MAX_TRIES attempts")
    end

    metrics = (:Res, :Rea, :PulseRT, :Pers, :Coll, :Nonnorm)

    Δ_full = Dict(m=>Float64[] for m in metrics)
    Δ_diag = deepcopy(Δ_full)

    for _ in 1:reps
        # --- get a new 'full' J1 that is stable ---
        J1 = nothing
        tries = 0
        while tries < MAX_TRIES
            tries += 1
            _, M1, _, _, _, _, _ = build_real_J(R, C; conn=0.2, IS=σ, scen=scen)
            J1_candidate = D0 * M1
            if is_locally_stable(J1_candidate)
                J1 = J1_candidate
                break
            end
        end
        if J1 === nothing 
            @warn "J1 never stabilized after $MAX_TRIES tries" 
        end
        J1 = J1 === nothing ? D0*M0 : J1  # fallback

        # --- get a new 'diag‐fixed' J2 that is stable ---
        J2 = nothing
        tries = 0
        while tries < MAX_TRIES
            tries += 1
            _, M2, _, _, _, _, _ = build_real_J(R, C; conn=0.2, IS=σ, scen=scen)
            J2_candidate = D0 * M2
            if is_locally_stable(J2_candidate)
                J2 = J2_candidate
                break
            end
        end
        if J2 === nothing
            @warn "J2 never stabilized after $MAX_TRIES tries"
        end
        J2 = J2 === nothing ? D0*M0 : J2  # fallback

        # --- now compute your Δ‐metrics against J0 ---
        push!(Δ_full[:Res],     abs(measure_resilience(J1) - measure_resilience(J0)))
        push!(Δ_full[:Rea],     abs(measure_reactivity(J1) - measure_reactivity(J0)))
        push!(Δ_full[:PulseRT], abs(pulse_return_time(J1) - pulse_return_time(J0)))
        push!(Δ_full[:Pers],    abs(measure_persistence(J1) - measure_persistence(J0)))
        push!(Δ_full[:Coll],    abs(measure_collectivity(-diag(J1),R+C;σ=σ) - measure_collectivity(-diag(J0),R+C;σ=σ)))
        push!(Δ_full[:Nonnorm], abs(nonnormality(J1) - nonnormality(J0)))

        push!(Δ_diag[:Res],     abs(measure_resilience(J2) - measure_resilience(J0)))
        push!(Δ_diag[:Rea],     abs(measure_reactivity(J2) - measure_reactivity(J0)))
        push!(Δ_diag[:PulseRT], abs(pulse_return_time(J2) - pulse_return_time(J0)))
        push!(Δ_diag[:Pers],    abs(measure_persistence(J2) - measure_persistence(J0)))
        push!(Δ_diag[:Coll],    Δ_full[:Coll][end])
        push!(Δ_diag[:Nonnorm], abs(nonnormality(J2) - nonnormality(J0)))
    end

    # compute means
    mean_full = Dict(k=>mean(v) for (k,v) in Δ_full)
    mean_diag = Dict(k=>mean(v) for (k,v) in Δ_diag)

    # plot
    titles = ["Δ Res","Δ Rea", "Δ PulseRT","Δ Pers","Δ Coll","Δ Nonnorm"]
    fig = Figure(; size=(1000,450))
    Label(fig[0, 1:3], "σ = $σ"; fontsize=14)
    for (i, m) in enumerate(metrics)
        row, col = divrem(i-1,4)
        ax = Axis(fig[row+1, col+1])
        xs = [1,2]
        ys = [mean_full[m], mean_diag[m]]
        barplot!(ax, xs, ys)
        ax.title = titles[i]
        ax.xticks = (xs, ["full","diag-fix"])
    end
    display(fig)

    return mean_full, mean_diag
end

# run it at several σ’s
for σ in (0.1, 1.0, 5.0)
    mf, md = compare_diag_effect(50,20; scen=:ER, σ=σ, reps=20)
    println("σ=$σ full vs diag-fixed:")
    @show mf; @show md
end