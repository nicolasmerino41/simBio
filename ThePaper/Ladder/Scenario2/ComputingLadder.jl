using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
include("Ladder4.1.jl")

# ────────────────────────────────────────────────────────────────
# 5) Stability & survival checks
# ────────────────────────────────────────────────────────────────────────────────
function is_locally_stable(J::AbstractMatrix)
    if any(!isfinite, J)
        return false
    end
    lambda = eigvals(J)
    maximum(real(lambda)) < 0
end

function survives!(fixed, p; tspan=(0.,500.), cb)
    prob = ODEProblem(trophic_ode!, fixed, tspan, p)
    sol  = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)
    return sol.t[end]<tspan[2] ? (false,sol.u[end]) : (all(sol.u[end] .> EXTINCTION_THRESHOLD),sol.u[end])
end

function calibrate_with_margins(
    A::AbstractMatrix, epsilon::AbstractMatrix, R::Int;
    margins=[0.05,0.2,0.5,1.0,2.0]
)
    S = size(A,1)
    C = S - R

    for delt in margins
        xi_cons, K_res = compute_default_thresholds(A, epsilon, R; margin=delt)
        B = try
            calibrate_from_K_xi(xi_cons, K_res, epsilon, A)
        catch
            continue  # singular
        end

        # positivity & finiteness
        if any(!isfinite, B) || any(x->x<=0, B)
            continue
        end

        # split out
        B_res = view(B, 1:R)
        B_cons = view(B, R+1:S)

        # build Jacobian and test local stability:
        # (this is a *placeholder*; replace with your real compute_jacobian!)
        D = Diagonal(1.0 ./ B)     # crude approx: D=diag(1/B_i)
        Mj = zeros(S,S)             # you'd fill this with the true partials
        J = D*Mj
        if any(!isfinite, J)
            continue
        end
        # Here we skip eigvals if Mj is zero; in your real code,
        # use `eigvals(J)` and test `maximum(real(lambda))<0`.

        # success
        return B_res, B_cons, xi_cons, K_res, delt
    end

    return nothing
end

# ────────────────────────────────────────────────────────────────────────────────
# Main sweep: fix S & C, vary structure, collect full & ladder pers.
# ────────────────────────────────────────────────────────────────────────────────
function ComputingLadder(
    S::Int=30, C::Int=6;
    conn_vals=[0.05,0.1,0.2],
    IS_vals=[0.1,1.0],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[0.1],
    delta_vals=[1.0,3.0],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100,
    threshold_steps=50
)
    R = S - C
    A = zeros(S,S)
    results = Vector{NamedTuple}()
    locki = ReentrantLock()
    cb = build_callbacks(S,EXTINCTION_THRESHOLD)

    combos = collect(Iterators.product(
        conn_vals,IS_vals,scenarios,delta_vals,eps_scales
    ))

    @threads for (conn,IS,scen,delta,eps) in 
        sample(combos, min(length(combos), number_of_combinations); replace=false)

        # 1) build A & epsilon
        A = make_A(A,R,conn,scen).*IS
        epsilon = clamp.(randn(S,S).*eps,0,1)

        # --- NEW: compute collectivity phi ---
        psi = compute_collectivity(A, epsilon)

        # 2) calibrate Xi,K → get (R_eq,C_eq,xi_cons,K_res)
        # 2) generate *all* feasible (Xi,K) sets for this A,epsilon
        thr_sets = generate_feasible_thresholds(A, epsilon, R;
                        margins = 10 .^(range(-2, 1, length=threshold_steps)),)
        isempty(thr_sets) && ( @warn "No feasible thresholds for $conn, $IS, $scen, $eps)" ; continue )

        # 2b) now loop over each feasible threshold‐set
        for t in thr_sets
            xi_cons   = t.xi_cons
            K_res     = t.K_res
            used_margin = t.margin
            R_eq        = t.R_eq
            C_eq        = t.C_eq

            # set up parameters & run the rest of your pipeline exactly as before
            m_cons = xi_cons
            d_res  = ones(R)
            r_res  = K_res

            p     = (R, C, m_cons, xi_cons, r_res, d_res, epsilon, A)
            fixed = vcat(R_eq, C_eq)

            # 3) stability & survival
            D, M = compute_jacobian(fixed, p)
            !is_locally_stable(D*M) && continue
            ok, B0 = survives!(fixed, p; cb=cb)
            !ok && continue

            # resilience_full is –max Re(eig(J)), but compute_resilience returns max Re
            resilience_full = compute_resilience(fixed, p)
            reactivity_full =  compute_reactivity(fixed, p)

            # 4) full‐model persistence
            rt_press, os, ire, before_full, after_persistence_full, new_equil, _ = simulate_press_perturbation(
                B0, p, tspan, tpert, delta;
                solver=Tsit5(),
                plot=false,
                show_warnings=true,
                full_or_simple=true,
                cb=cb,
                species_specific_perturbation=false
            )
            # 4a) full‐model pulse
            rt_pulse, _, _, _, _ = simulate_pulse_perturbation(
                B0, p, tspan, tpert, delta;
                solver=Tsit5(),
                plot=false,
                cb=cb,
                species_specific_perturbation=false
            )
            rt_press_full = mean(filter(!isnan, rt_press))
            rt_pulse_full = mean(filter(!isnan, rt_pulse))

            collectivity_full = psi 

            # 5) ladder persistence
            before_persistence_S = Dict(i => NaN for i in 1:16)
            after_persistence_S  = Dict(i => NaN for i in 1:16)
            rt_press_S   = Dict(i => NaN for i in 1:16)
            rt_pulse_S   = Dict(i => NaN for i in 1:16)
            collectivity_S = Dict(i => NaN for i in 1:16)
            resilience_S  = Dict(i=>NaN for i in 1:16)
            reactivity_S  = Dict(i=>NaN for i in 1:16)
            @info "Running ladder"

            for step in 1:16
                # build simplified A and epsilon
                A_s, epsilon_s = transform_for_ladder_step(step, A, epsilon)

                # original equilibrium abundances
                R_eq_full, C_eq_full = R_eq, C_eq

                psi_s = compute_collectivity(A_s, epsilon_s)

                # 5a) Recompute xi_hat
                xi_hat = zeros(C)
                for k in 1:C
                    i = R + k
                    # 1) feeding gains (A>0)
                    gain = sum(epsilon_s[i,j]*A_s[i,j]*R_eq_full[j] for j in 1:R if A_s[i,j] > 0; init=0.0 )
                    # 2) predation losses: A_s[i,j]<0, but we need "+(-A)B"
                    loss = sum(A_s[i,j]*C_eq_full[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                    # consumer eq: xi = B_i – gain – loss
                    xi_hat[k] = -C_eq_full[k] + gain + loss
                end

                # 5b) Recompute K_hat
                K_hat = zeros(R)
                for i in 1:R
                    # resource eq uses A[i,j] (j consumer) directly:
                    drain = sum(A_s[i,j]*C_eq_full[j-R] for j in R+1:S if A_s[i,j] < 0; init=0.0 )
                    # K_i = B_i + ∑ A[i,j] B_j
                    K_hat[i] = R_eq_full[i] + drain
                end

                # 5c) Solve for new equilibrium
                eq = calibrate_from_K_xi(xi_hat, K_hat, epsilon_s, A_s)
                if eq === nothing
                    @warn("Ladder step $step infeasible: xi_hat=$(xi_hat), K_hat=$(K_hat), step=$step")
                    continue
                end
                R_eq_s, C_eq_s = eq

                # 5d) simulate simplified model
                p_s = (R, C, xi_hat, xi_hat, K_hat, ones(R), epsilon_s, A_s)
                rt_press2, _, _, before_s, after_s, _, _ = simulate_press_perturbation(
                    B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    full_or_simple=false
                )
                rt_pulse3, _, _, _, _ = simulate_pulse_perturbation(
                    B0, p_s, tspan, tpert, delta;
                    solver=Tsit5(),
                    plot=false,
                    cb=cb,
                    species_specific_perturbation=false
                )

                before_persistence_S[step] = before_s
                after_persistence_S[step]  = after_s
                rt_press_S[step]   = mean(filter(!isnan, rt_press2))
                rt_pulse_S[step]   = mean(filter(!isnan, rt_pulse3))
                collectivity_S[step] = psi_s
                resilience_S[step] = compute_resilience(B0, p_s)
                reactivity_S[step] = compute_reactivity(B0, p_s)
            end

            step_pairs = collect(Iterators.flatten(
                ([
                    Symbol("before_persistence_S$i") => before_persistence_S[i],
                    Symbol("aft_S$i") => after_persistence_S[i],
                    Symbol("rt_press_S$i") => rt_press_S[i],
                    Symbol("rt_pulse_S$i") => rt_pulse_S[i],
                    Symbol("collectivity_S$i") => collectivity_S[i],
                    Symbol("resilience_S$i") => resilience_S[i],
                    Symbol("reactivity_S$i") => reactivity_S[i]
                ] for i in 1:16)
            ))

            rec = (
                conn=conn, IS=IS, scen=scen, margin=used_margin,
                delta =delta, eps=eps,
                before_persistence_full=before_full, after_persistence_full=after_persistence_full, rt_press_full=rt_press_full, rt_pulse_full=rt_pulse_full,
                collectivity_full=collectivity_full, resilience_full=resilience_full, reactivity_full=reactivity_full,
                step_pairs...,  # Properly flattened pairs
                p_final = p,
                R_eq = R_eq,
                C_eq = C_eq
            )

            lock(locki) do
                push!(results, rec)
            end
        end
    end

    return DataFrame(results)
end


# ────────────────────────────────────────────────────────────────────────────────
# 9) Run it
# ────────────────────────────────────────────────────────────────────────────────
A = ComputingLadder(
    50, 20;
    conn_vals=0.01:0.05:0.9,
    IS_vals=[0.01, 0.1, 1.0, 2.0],
    scenarios=[:ER,:PL,:MOD],
    delta_vals=[1.0,3.0],
    eps_scales=[1.0, 0.1],
    tspan=(0.,500.), tpert=250.0,
    number_of_combinations = 100,
    threshold_steps=2,
)


A = deserialize("ThePaper/Ladder/Outputs/Final_results.jls")

function plot_step_correlations(
    df::DataFrame,
    var::Symbol;             # :before_persistence, :after_persistence or :rt
    color_by::Symbol = :conn # e.g. :conn, :IS, :scen, :delta, :R_eq, :C_eq, etc.
)
    # 1) your keys & labels
    step_keys = ["S1","S2","S3","S4","S5","S6","S7","S8","S9","S10",
                 "S11","S12","S13","S14","S15"]
    step_names = [
        "Full Model",
        "Full A (Species-specific ϵ)", "Full A (Global ϵ)", "Full A (Re-randomised ϵ)",
        "Global -Aij & Aij (ϵ_full)", "Global -Aij & Aij (Species-specific ϵ)",
        "Global -Aij & Aij (Global ϵ)", "Global -Aij & Aij (Re-randomised ϵ)",
        "Global A (ϵ_full)", "Global A (Species-specific ϵ)",
        "Global A (Global ϵ)", "Global A (Re-randomised ϵ)",
        "Re-randomised A (ϵ_full)", "Re-randomised A (Species-specific ϵ)",
        "Re-randomised A (Global ϵ)", "Re-randomised A (Re-randomised ϵ)"
    ]
    @assert length(step_keys)+1 == length(step_names)

    # 2) assemble the panel‐column symbols
    full_col   = Symbol(string(var), "_full")
    panel_cols = [ full_col;
                   map(k -> Symbol(string(var), "_", k), step_keys) ]

    # 3) pull data vectors
    full_vals  = df[!, full_col]
    
    color_vals = df[!, color_by]
    cmin, cmax = minimum(color_vals), maximum(color_vals)

    # 4) grid layout (16 panels → 4×4)
    ns   = length(panel_cols)
    rows = ceil(Int, sqrt(ns))
    cols = ceil(Int, ns/rows)

    fig = Figure(; size = (900, 650))
    # optional big title:
    Label(fig[0, 1:cols], uppercase(string(var, " correlations")); fontsize=24)

    # 5) loop panels
    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1

        ax = Axis(
            fig[r, c];
            xlabel = string(full_col),
            ylabel = string(panel_cols[idx]),
            title  = step_names[idx],
            titlesize = 10,
            xlabelsize = 12,
            ylabelsize = 12,
        )

        step_vals = df[!, panel_cols[idx]]

        scatter!(
            ax, full_vals, step_vals;
            markersize   = 6,
            color        = color_vals,
            colormap     = :viridis,
            colorrange   = (cmin, cmax)
        )

        mn = min(minimum(full_vals), minimum(step_vals))
        mx = max(maximum(full_vals), maximum(step_vals))
        lines!(ax, [mn,mx], [mn,mx]; color=:black, linestyle=:dash)

        r_val = cor(full_vals, step_vals)
        text!(
            ax, "r=$(round(r_val, digits=2))";
            position = (mx, mn),
            align    = (:right, :bottom),
            fontsize = 9
        )
    end

    # # 6) single colorbar on the right
    # cbax = Axis(fig[1:rows, cols+1], showaxis=false)
    # Colorbar(cbax;
    #     colormap = :viridis,
    #     limits   = (cmin, cmax),
    #     label    = string(color_by)
    # )

    display(fig)
    return fig
end

# e.g. colour by connectance
plot_step_correlations(A, :rt_press;  color_by = :conn)

# or colour by number of surviving species after full press
plot_step_correlations(A, :before_persistence; color_by = :conn)

# or by scenario
plot_step_correlations(A, :after_persistence;  color_by = :conn)

# or by delta
plot_step_correlations(A, :rt_pulse;)

# or by equilibrium abundances
plot_step_correlations(A, :collectivity;)

# or by resilience
plot_step_correlations(A, :resilience;)

# or by reactivity
plot_step_correlations(A, :reactivity;)
