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
    find_stable_configuration(S, O, R; conn, mu_range, eps_range, m_alpha_range, time_end)

Loops over specified ranges for μ, ε, and mₐ, builds the abstract community and simulates it,
and returns the first configuration that yields a locally stable equilibrium (i.e. all eigenvalues of J have negative real parts).
"""
function find_stable_configurations(S::Int, O::Int, R::Int; conn=0.2,
        mu_range=0.0:0.1:1.0, eps_range=0.0:0.1:1.0, m_alpha_range=0.05:0.05:0.2,
        time_end=500.0)
    stable_configs = []
    for mu_val in mu_range
        for eps_val in eps_range
            for m_alpha_val in m_alpha_range
                comm = a_parametrise_the_community(S, O, R;
                    conn=conn, mu=mu_val, epsilon_val=eps_val, mean_m_alpha=m_alpha_val,
                    cell_abundance_h = [i <= S ? rand()*10.0 : rand()*5.0 for i in 1:S+O],  # if provided, else defaults will be ones
                    cell_abundance_p = [rand() for i in 1:R]
                    )
                u0 = vcat(comm.H_eq, comm.P_eq)
                params = (comm.S, comm.R, comm.K_i, comm.r_i, mu_val, comm.nu,
                          comm.P_matrix, comm.O_matrix, comm.T_matrix, comm.epsilon,
                          comm.m_alpha, comm.K_alpha, comm.B_matrix, comm.D_matrix,
                          comm.nu*1.0, comm.nu*1.0)
                sol, J = simulate_community(u0, params, time_end)
                if sol[:, end] == u0 #all(real.(eigen(J).values) .< 0) #
                    push!(stable_configs,
                        (comm=comm, sol=sol, J=J, mu=mu_val, eps=eps_val, m_alpha=m_alpha_val))
                end
            end
        end
    end
    if isempty(stable_configs)
        error("No stable configuration found for H=$S, O=$O, R=$R, conn=$conn")
    else
        return stable_configs
    end
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
                        stable = find_stable_configurations(S, O, R; conn=conn, mu_range=mu_range,
                                    eps_range=eps_range, m_alpha_range=m_alpha_range,
                                    time_end=time_end)
                        for s in stable
                            sol, J = s.sol, s.J
                            resil = compute_resilience(J)
                            react = compute_reactivity(J)
                            persist = compute_persistence(sol)
                            push!(results, (S, O, R, conn, s.mu, s.eps, s.m_alpha, resil, react, persist))
                        end
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
function plot_experiment_results(results::DataFrame; altogether=true, averages=false)
    # Helper: if data are constant, force a small range for axis limits.
    function fix_data_range(data::Vector{<:Real})
        if minimum(data) == maximum(data)
            return data .+ (-0.5:0.1:0.5)[1:length(data)]
        else
            return data
        end
    end

    # Compute derived quantities.
    # Proportion Omnivores: handle O==0 case.
    prop_omn = results.O ./ (results.S .+ results.O)
    prop_omn = fix_data_range(prop_omn)
    # Total species
    total_species = results.S .+ results.O .+ results.R

    if !altogether
        # Plot each scatter separately.
        fig1 = Figure(resolution=(800,600))
        ax1 = Axis(fig1[1,1], xlabel="Connectance", ylabel="Resilience", title="Resilience vs. Connectance")
        scatter!(ax1, results.conn, results.resilience, markersize=8, color=:blue)

        fig2 = Figure(resolution=(800,600))
        ax2 = Axis(fig2[1,1], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
        scatter!(ax2, prop_omn, results.reactivity, markersize=8, color=:green)

        fig3 = Figure(resolution=(800,600))
        ax3 = Axis(fig3[1,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Species Richness")
        scatter!(ax3, total_species, results.persistence, markersize=8, color=:red)

        display(fig1)
        display(fig2)
        display(fig3)
    elseif altogether && !averages
        # Plot all points together in one figure with four subplots.
        fig1 = Figure(resolution=(800,600))
        ax1 = Axis(fig1[1,1], xlabel="Connectance", ylabel="Resilience", title="Resilience vs. Connectance")
        scatter!(ax1, results.conn, results.resilience, markersize=8, color=:blue)

        ax2 = Axis(fig1[1,2], xlabel="Proportion Omnivores", ylabel="Reactivity", title="Reactivity vs. Proportion Omnivores")
        scatter!(ax2, prop_omn, results.reactivity, markersize=8, color=:green)

        ax3 = Axis(fig1[2,1], xlabel="Total Species", ylabel="Persistence", title="Persistence vs. Species Richness")
        scatter!(ax3, total_species, results.persistence, markersize=8, color=:red)

        ax4 = Axis(fig1[2,2], xlabel="Total Species", ylabel="Resilience", title="Resilience vs. Species Richness")
        scatter!(ax4, total_species, results.resilience, markersize=8, color=:red)

        display(fig1)
    elseif altogether && averages
        
        # For resilience vs. connectance: group by connectance.
        grouped_conn = groupby(results, :conn)
        conn_vals = [first(g.conn) for g in grouped_conn]
        mean_resilience = [mean(g.resilience) for g in grouped_conn]
        sd_resilience   = [std(g.resilience) for g in grouped_conn]

        # For reactivity vs. proportion omnivores: group by O (assume S is constant so prop_omn = O/(S+O))
        grouped_O = groupby(results, :O)
        mean_prop_omn = [mean(g.O ./ (g.S .+ g.O)) for g in grouped_O]
        mean_reactivity = [mean(g.reactivity) for g in grouped_O]
        sd_reactivity   = [std(g.reactivity) for g in grouped_O]

        # For persistence and resilience vs. total species: add a column.
        results[!, :total_species] = results.S .+ results.O .+ results.R
        grouped_ts = groupby(results, :total_species)
        total_species_vals = [first(g.total_species) for g in grouped_ts]
        mean_persistence = [mean(g.persistence) for g in grouped_ts]
        sd_persistence   = [std(g.persistence) for g in grouped_ts]
        mean_resilience_ts = [mean(g.resilience) for g in grouped_ts]
        sd_resilience_ts   = [std(g.resilience) for g in grouped_ts]
        
        fig1 = Figure(resolution=(800,600))
        ax1 = Axis(fig1[1,1], xlabel="Connectance", ylabel="Resilience",
                   title="Resilience vs. Connectance")
        MK.errorbars!(ax1, conn_vals, mean_resilience, sd_resilience, color = :blue)
        scatter!(ax1, conn_vals, mean_resilience, color = :blue)

                # Compute and store proportion omnivores
        results[!, :prop_omn] = results.O ./ (results.S .+ results.O)

        # Group by rounded prop_omn (or bin them if you want)
        grouped_prop = groupby(results, :prop_omn)
        prop_vals = [first(g.prop_omn) for g in grouped_prop]
        mean_reactivity = [mean(g.reactivity) for g in grouped_prop]
        sd_reactivity   = [std(g.reactivity) for g in grouped_prop]

        # Sort by prop_omn for clean x-axis
        sorted_idx = sortperm(prop_vals)
        prop_vals_sorted = prop_vals[sorted_idx]
        mean_reactivity_sorted = mean_reactivity[sorted_idx]
        sd_reactivity_sorted = sd_reactivity[sorted_idx]

        ax2 = Axis(fig1[2,2], xlabel="Proportion Omnivores", ylabel="Reactivity",
        title="Reactivity vs. Proportion Omnivores")
        errorbars!(ax2, prop_vals_sorted, mean_reactivity_sorted, sd_reactivity_sorted, color = :green)
        scatter!(ax2, prop_vals_sorted, mean_reactivity_sorted, color = :green)

        ax3 = Axis(fig1[2,1], xlabel="Total Species", ylabel="Persistence",
                   title="Persistence vs. Species Richness")
        errorbars!(ax3, total_species_vals, mean_persistence, sd_persistence, color = :red)
        scatter!(ax3, total_species_vals, mean_persistence, color = :red)

        ax4 = Axis(fig1[1,2], xlabel="Total Species", ylabel="Resilience",
                   title="Resilience vs. Species Richness")
        errorbars!(ax4, total_species_vals, mean_resilience_ts, sd_resilience_ts, color = :purple)
        scatter!(ax4, total_species_vals, mean_resilience_ts, color = :purple)

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
# AA_results = run_experiment(
#     S_range=1:5:30, O_range=0:5:10, R_range=0:5:10, conn_range=0.1:0.2:1.0,
# )
# println(A_results)
plot_experiment_results(A_results; averages = true)
plot_3D_results(A_results)

begin
    results = deepcopy(A_results)
    results[!, :community] = string.("S", results.S, "_O", results.O, "_R", results.R)

    results[!, :feasible] .= true

    grouped = groupby(results, :community)
    community_names = [first(g.community) for g in grouped]
    n_total = [nrow(g) for g in grouped]
    n_feasible = [sum(g.feasible) for g in grouped]
    feasibility_ratio = n_feasible ./ n_total

    fig = Figure(resolution=(1000, 600))
    ax = Axis(fig[1, 1], xlabel="Community", ylabel="Feasible Proportion", 
            title="Feasible Parameter Space per Community")
    barplot!(ax, community_names, feasibility_ratio, color=:skyblue, rotation=π/4)
    fig
end

function plot_feasibility_heatmaps(results::DataFrame;
    S_vals=1:5:30, O_vals=0:5:10, R_vals=0:5:10,
    mu_vals=unique(results.mu), eps_vals=unique(results.eps))

    # Get unique combinations of communities
    communities = [(S, O, R) for S in S_vals, O in O_vals, R in R_vals if S + O + R > 0]

    # Build lookup set from actual results
    results[!, :key] = string.("S", results.S, "_O", results.O, "_R", results.R, "_mu", results.mu, "_eps", results.eps)
    existing_keys = Set(results.key)

    # Prepare grid for layout
    ncols = 4
    nrows = ceil(Int, length(communities) / ncols)
    fig = Figure(resolution=(ncols * 300, nrows * 300))

    for (i, (S, O, R)) in enumerate(communities)
        row = ceil(Int, i / ncols)
        col = i % ncols == 0 ? ncols : i % ncols

        # Build matrix of 1/0 for mu x eps
        z = zeros(Bool, length(mu_vals), length(eps_vals))

        for (i_mu, mu) in enumerate(mu_vals), (i_eps, eps) in enumerate(eps_vals)
            key = "S$(S)_O$(O)_R$(R)_mu$(mu)_eps$(eps)"
            if key in existing_keys
                z[i_mu, i_eps] = true
            end
        end

        ax = Axis(fig[row, col],
            title = "S=$S O=$O R=$R",
            xlabel = "eps",
            ylabel = "mu",
            xticks = (1:length(eps_vals), string.(eps_vals)),
            yticks = (1:length(mu_vals), string.(mu_vals))
        )
        heatmap!(ax, z; colormap = [:white, :blue])
    end

    fig
end


plot_feasibility_heatmaps(A_results)