# =============================================================================
# 2. ODE Dynamics and Simulation (using your omnivore_dynamics!)
# =============================================================================
# Here we assume that "omnivore_dynamics!" is defined elsewhere and works.
# We wrap it for our simulation purposes.
function simulate_community(u0, params, time_end=500.0; initial_deviation=0.0)
    u0 = u0 .+ randn(length(u0)) .* initial_deviation
    prob = ODEProblem(omnivore_dynamics!, u0, (0.0, time_end), params)
    sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
    function f_wrapper(u)
        du = similar(u)
        omnivore_dynamics!(du, u, params, 0.0)
        return du
    end
    J = ForwardDiff.jacobian(f_wrapper, sol.u[end])
    return sol, J
end

# =============================================================================
# 3. Metrics Computation
# =============================================================================
function compute_resilience(J)
    # Resilience approximated as the spectral abscissa (max real part of eigenvalues)
    return maximum(real.(eigen(J).values))
end

function compute_reactivity(J)
    S_sym = (J + transpose(J)) / 2
    return maximum(real.(eigen(S_sym).values))
end

function compute_persistence(sol; threshold=1e-6)
    final_state = sol.u[end]
    n = length(final_state)
    survivors = count(x -> x > threshold, final_state)
    return survivors / n
end

# =============================================================================
# 4. Search for a Stable Configuration
# =============================================================================
"""
    find_stable_configuration(H, O, R; conn, mu_range, eps_range, m_alpha_range, time_end)

Loops over specified ranges for μ, ε, and mₐ, builds the abstract community and simulates it,
and returns the first configuration that yields a locally stable equilibrium (i.e. all eigenvalues of J have negative real parts).
"""
function find_stable_configuration(H::Int, O::Int, R::Int; conn=0.2,
        mu_range=0.4:0.1:0.6, eps_range=0.8:0.1:1.2, m_alpha_range=0.05:0.05:0.2,
        mean_m_alpha=0.1, time_end=500.0, cell_abundance_h=nothing, cell_abundance_p=nothing)
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = a_parametrise_the_community(H, O, R;
                    conn=conn, mu=mu_val, epsilon_val=eps_val, mean_m_alpha=m_alpha_val,
                    cell_abundance_h=cell_abundance_h, cell_abundance_p=cell_abundance_p)
                u0 = vcat(comm.H_eq, comm.P_eq)
                params = (comm.S, comm.R, comm.K_i, comm.r_i, mu_val, comm.nu,
                          comm.P_matrix, comm.O_matrix, comm.T_matrix, comm.epsilon,
                          comm.m_alpha, comm.K_alpha, comm.B_matrix, comm.D_matrix,
                          comm.nu*1.0, comm.nu*1.0)
                sol, J = simulate_community(u0, params, time_end)
                if all(real.(eigen(J).values) .< 0)
                    println("Stable configuration found: mu=$(mu_val), eps=$(eps_val), m_alpha=$(m_alpha_val)")
                    return (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val)
                end
            end
        end
    end
    error("No stable configuration found for H=$H, O=$O, R=$R, conn=$conn")
end

# =============================================================================
# 5. Experiment Pipeline: Loop over Community Configurations
# =============================================================================
"""
    run_experiment(; S_range, O_range, R_range, conn_range, ...)

For each combination of species numbers and connectance, uses find_stable_configuration
to obtain a locally stable community, then computes metrics (resilience, reactivity, persistence).
Returns a DataFrame with the results.
"""
function run_experiment(
    ; S_range=10:10:30, O_range=0:10:30, R_range=5:5:15, conn_range=0.1:0.1:1.0,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:0.2,
    time_end=500.0
)
    results = DataFrame(S=Int[], O=Int[], R=Int[], conn=Float64[],
                        mu=Float64[], eps=Float64[], m_alpha=Float64[],
                        resilience=Float64[], reactivity=Float64[], persistence=Float64[])
    for S in S_range
        for O in O_range
            for R in R_range
                for conn in conn_range
                    try
                        stable = find_stable_configuration(S, O, R; conn=conn, mu_range=mu_range,
                                                eps_range=eps_range, m_alpha_range=m_alpha_range,
                                                mean_m_alpha=mean_m_alpha, time_end=time_end)
                        sol, J = stable.sol, stable.J
                        resil = compute_resilience(J)
                        react = compute_reactivity(J)
                        persist = compute_persistence(sol)
                        push!(results, (S, O, R, conn, stable.mu, stable.eps, stable.m_alpha, resil, react, persist))
                    catch e
                        println("Skipping configuration S=$S, O=$O, R=$R, conn=$conn: ", e)
                    end
                end
            end
        end
    end
    return results
end

# =============================================================================
# 6. Plotting the Experiment Results
# =============================================================================
# Then, you can plot the results using a helper function:
function plot_experiment_results(results::DataFrame; altogether=true)
    if !altogether
        fig1 = Figure(resolution=(800,600))
        ax1 = Axis(fig1[1,1], xlabel="Connectance", ylabel="Resilience", title="Resilience vs. Connectance")
        MK.scatter!(ax1, results.conn, results.resilience, markersize=8, color=:blue)

        prop_omn = results.O ./ (results.S .+ results.O)
        fig2 = Figure(resolution=(800,600))
        ax2 = Axis(fig2[1,1], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
        MK.scatter!(ax2, prop_omn, results.reactivity, markersize=8, color=:green)

        total_species = results.S .+ results.O .+ results.R
        fig3 = Figure(resolution=(800,600))
        ax3 = Axis(fig3[1,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Species Richness")
        MK.scatter!(ax3, total_species, results.persistence, markersize=8, color=:red)

        display(fig1)
        display(fig2)
        display(fig3)
    else
        fig1 = Figure(resolution=(800,600))
        ax1 = Axis(fig1[1,1], xlabel="Connectance", ylabel="Resilience", title="Resilience vs. Connectance")
        MK.scatter!(ax1, results.conn, results.resilience, markersize=8, color=:blue)

        prop_omn = results.O ./ (results.S .+ results.O)
        ax2 = Axis(fig1[1,2], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
        MK.scatter!(ax2, prop_omn, results.reactivity, markersize=8, color=:green)

        total_species = results.S .+ results.O .+ results.R
        ax3 = Axis(fig1[2,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Species Richness")
        MK.scatter!(ax3, total_species, results.persistence, markersize=8, color=:red)

        ax4 = Axis(fig1[2,2], xlabel="Total Species", ylabel="Resilience", title="Resilience vs. Species Richness")
        MK.scatter!(ax4, total_species, results.resilience, markersize=8, color=:red)

        display(fig1)
    end
end

function plot_3D_results(results::DataFrame)

    total_species = results.S .+ results.O .+ results.R
    fig1 = Figure(resolution=(800,600))
    ax1 = Axis3(fig1[1,1], xlabel="Connectance", ylabel="Resilience", zlabel="Species Richness", title="Connectance vs. Resilience vs. Species Richness")
    MK.scatter!(ax1, results.conn, results.resilience, total_species, markersize=8, color=:blue)

    display(fig1)
end

# -------------------------------
# 7. Running the Pipeline
# -------------------------------
A_results = run_experiment()
println(A_results)
plot_experiment_results(A_results)
plot_3D_results(A_results)