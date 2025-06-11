using Random, LinearAlgebra, Statistics, DataFrames, CairoMakie

# --- Helper functions from pipeline ---
function compute_jacobian(B, p)
    R, C, m_cons, xi_cons, r_res, K_res, ε, A = p
    S = R + C
    # Build D
    D = zeros(S,S)
    for i in 1:R
        D[i,i] = r_res[i]/K_res[i] * B[i]
    end
    for k in 1:C
        i = R + k
        D[i,i] = (m_cons[k]/xi_cons[k]) * B[i]
    end
    # Build Mstar = -I + A*
    Mstar = Matrix{Float64}(I, S, S)
    for i in 1:R, j in 1:S
        if i!=j && A[i,j] != 0
            Mstar[i,j] += A[i,j]
        end
    end
    for i in R+1:S, j in 1:S
        if A[i,j] > 0
            Mstar[i,j] += ε[i,j]*A[i,j]
        elseif A[i,j] < 0
            Mstar[i,j] += A[i,j]
        end
    end
    return D*Mstar
end

# Extract sigma/min(d)
function sigma_over_min_d(A, J)
    d = -diag(J)             # from J = D(-I + A)
    min_d = minimum(d)
    offs = [A[i,j] for i in 1:size(A,1), j in 1:size(A,1) if i!=j]
    σ = std(offs)
    return σ/min_d
end

# Rewire A in place
function rewire_A!(A, change_connectance=false, new_conn=0.1, σ=1.0)
    S = size(A,1)
    # optionally change sparsity
    if change_connectance
        mask = rand(S,S) .< new_conn
        mask .= mask .& .!I(S)
    else
        mask = A .!= 0
    end
    # resample nonzero entries
    for i in 1:S, j in 1:S
        if mask[i,j] && i!=j
            A[i,j] = randn()*σ
        end
    end
    return A
end

# --- Main experiment ---
function test_rewiring_pipeline(; S=50, C=20, conn=0.1, σ=1.0, reps=100)
    R = S - C
    results = DataFrame(step=Int[], resilience=Float64[], reactivity=Float64[],
                        sl=Float64[], sigma_ratio=Float64[])
    Random.seed!(42)

    for rep in 1:reps
        # 1) Sample parameters
        A = zeros(S,S)
        mask = rand(S,S) .< conn
        mask .= mask .& .!I(S)
        A[mask] .= randn(sum(mask))*σ
        # sample epsilon arbitrarily
        ε = clamp.(randn(S,S)*0.1, 0, 1)
        # feasible thresholds
        xi_cons = abs.(randn(C)*0.5 .+ 1.0)
        K_res   = abs.(randn(R)*0.5 .+ 1.0)
        # equilibrium abundances (toy: B_i=K_i)
        B = vcat(K_res, xi_cons)
        # build J
        p = (R, C, fill(1.0,C), xi_cons, fill(1.0,R), K_res, ε, A)
        J = compute_jacobian(B, p)

        # Compute metrics at step 0 (original)
        push!(results, (0,
                        measure_resilience(J),
                        measure_reactivity(J),
                        mean(compute_SL(A, vcat(K_res, xi_cons))),
                        sigma_over_min_d(A, J)))

        # Step 1: simple rewiring
        A1 = copy(A)
        rewire_A!(A1, false, conn, σ)
        J1 = compute_jacobian(B, (R,C,fill(1.0,C),xi_cons,fill(1.0,R),K_res,ε,A1))
        push!(results, (1,
                        measure_resilience(J1),
                        measure_reactivity(J1),
                        mean(compute_SL(A1, vcat(K_res, xi_cons))),
                        sigma_over_min_d(A1, J1)))

        # Step 2: rewiring + new connectance & strength
        A2 = copy(A)
        new_conn = conn * 2.0  # e.g. double connectance
        new_σ    = σ * 2.0     # e.g. double strength
        rewire_A!(A2, true, new_conn, new_σ)
        J2 = compute_jacobian(B, (R,C,fill(1.0,C),xi_cons,fill(1.0,R),K_res,ε,A2))
        push!(results, (2,
                        measure_resilience(J2),
                        measure_reactivity(J2),
                        mean(compute_SL(A2, vcat(K_res, xi_cons))),
                        sigma_over_min_d(A2, J2)))
    end

    # Plot
    fig = Figure(resolution=(800, 300))
    # List of (column, label)
    metrics = [
        (:resilience, "Resilience"),
        (:reactivity, "Reactivity"),
        (:sl,         "SL"),
        (:sigma_ratio,"σ/min(d)")
    ]

    for (i, (col, ylabel)) in enumerate(metrics)
        ax = Axis(fig[1, i];
            title=ylabel,
            xlabel="Step",
            ylabel=ylabel
        )
        for step in 0:2
            vals = results[results.step .== step, col]
            scatter!(ax, fill(step, length(vals)), vals; alpha=0.3)
            lines!(ax, [step-0.2, step+0.2], [mean(vals), mean(vals)]; color=:black, linewidth=2)
        end
    end

    display(fig)
    return results
end

res = test_rewiring_pipeline()
