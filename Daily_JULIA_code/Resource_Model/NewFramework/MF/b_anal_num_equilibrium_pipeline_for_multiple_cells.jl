using DifferentialEquations, ForwardDiff, LinearAlgebra, DataFrames, CSV

# --- Analytical equilibrium function ---
function analytical_equilibrium(
    cell, 
    mu_val, eps_val, mean_m_alpha;
    delta_nu = 0.05,
    d_alpha = 1.0, d_i = 1.0,
    include_predators = true,
    sp_removed_name = nothing,
    artificial_pi = false, pi_size = 1.0,
    H_init = nothing,
    P_init = nothing
)
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # Filter species if needed.
    if !include_predators
        sp_nm = [sp for sp in sp_nm if sp in herbivore_names]
    end
    if !isnothing(sp_removed_name)
        sp_nm = [sp for sp in sp_nm if sp != sp_removed_name]
    end
    
    # Get community parameters.
    params_setup = b_attempt_setup_community(
        local_i, local_j,
        mu_val,
        eps_val,
        mean_m_alpha;
        species_names = sp_nm,
        artificial_pi = artificial_pi, pi_size = pi_size,
        delta_nu = delta_nu,
        d_alpha = d_alpha,
        d_i = d_i
    )
    if isnothing(params_setup)
        @error "Error: params_setup is nothing for cell $cell"
        return nothing
    end
    
    # Destructure parameters.
    S          = params_setup.S
    R          = params_setup.R
    H_eq       = params_setup.H_eq
    P_eq       = params_setup.P_eq
    r_i        = params_setup.r_i
    K_i        = params_setup.K_i
    nu_val     = params_setup.nu
    P_matrix   = params_setup.P_matrix
    epsilon    = params_setup.epsilon
    m_alpha    = params_setup.m_alpha
    K_alpha    = params_setup.K_alpha
    H_star     = params_setup.H_star
    P_star     = params_setup.P_star
    
    # Set initial conditions.
    if isnothing(H_init)
        H_init = H_star
    end
    if R > 0
        if isnothing(P_init)
            P_init = P_star
        end
        u0 = vcat(H_init, P_init)
    else
        u0 = H_init
    end
    # Uncomment for debugging:
    # println("Cell $cell: typeof(u0) = ", typeof(u0))
    
    # Build parameter tuple for the ODE.
    p = (S, R, K_i, r_i, mu_val, nu_val, P_matrix, epsilon, m_alpha, K_alpha)
    
    # f_wrapper for Jacobian computation.
    function f_wrapper(u)
        du = similar(u)
        bipartite_dynamics!(du, u, p, 0.0)
        return du
    end
    
    # Compute the Jacobian at equilibrium using ForwardDiff.
    J = ForwardDiff.jacobian(f_wrapper, u0)
    # Uncomment for debugging:
    # println("Cell $cell: Jacobian = ", J)
    
    return (
        equilibrium = (H_star = H_star, P_star = P_star, u0 = u0),
        parameters  = (
            S = S,
            R = R,
            r_i = r_i,
            K_i = K_i,
            mu = mu_val,
            nu = nu_val,
            P_matrix = P_matrix,
            epsilon = epsilon,
            m_alpha = m_alpha,
            K_alpha = K_alpha
        ),
        Jacobian = J
    )
end

# --- Sensitivity metrics function ---
function compute_sensitivity_metrics(A_eq, A_p; perturbation=0.01, tspan=(0.0, 50.0))
    # Extract equilibrium state and total equilibrium biomass.
    u0 = A_eq.u0
    total_eq = sum(u0)
    n = length(u0)
    sensitivity = zeros(n)
    
    # For each species, perturb its equilibrium value by a small fraction,
    # simulate the system, and compute the maximum absolute deviation
    # in total biomass relative to the equilibrium total.
    for i in 1:n
         u0_perturbed = copy(u0)
         u0_perturbed[i] += perturbation * u0[i]
         prob = ODEProblem(bipartite_dynamics!, u0_perturbed, tspan, A_p)
         sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
         deviations = [abs(sum(sol(u)) - total_eq) for u in sol.t]
         sensitivity[i] = maximum(deviations)
    end
    
    # Compute species traits.
    biomass = copy(u0)
    
    # "Degree": For herbivores (indices 1:S) the degree is the number of predators feeding on them,
    # and for predators (indices S+1:n) it is the number of prey they feed on.
    S = A_p.S
    R = A_p.R
    P_matrix = A_p.P_matrix
    degree = zeros(n)
    for i in 1:S
        degree[i] = sum(P_matrix[i, :])
    end
    for alpha in 1:R
        degree[S+alpha] = sum(P_matrix[:, alpha])
    end
    
    return (sensitivity = sensitivity, biomass = biomass, degree = degree)
end

# --- Function to plot relationships (using GLMakie) ---
function plot_sensitivity_traits(metrics)
    fig = Figure(resolution = (900, 400))
    
    ax1 = Axis(fig[1,1], xlabel = "Equilibrium Biomass", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Equilibrium Biomass")
    scatter!(ax1, metrics.biomass, metrics.sensitivity, markersize = 10)
    
    ax2 = Axis(fig[1,2], xlabel = "Degree (Number of Interactions)", ylabel = "Sensitivity (max deviation)",
               title = "Sensitivity vs Degree")
    scatter!(ax2, metrics.degree, metrics.sensitivity, markersize = 10)
    
    return fig
end

# --- Main Pipeline: Process multiple cells and collect results ---
function process_cells(cell_ids; mu_val=0.5, eps_val=1.0, mean_m_alpha=0.1,
                       delta_nu=0.05, d_alpha=1.0, d_i=1.0,
                       include_predators=true, artificial_pi=false, pi_size=1.0,
                       perturbation=0.01, tspan=(0.0,50.0))
    
    # DataFrame to store species-level results across cells.
    df = DataFrame(cell = Int[], species = Int[], species_type = String[],
                   sensitivity = Float64[], biomass = Float64[], degree = Float64[],
                   max_eigen = Float64[], stable = Bool[])
    
    for cell in cell_ids
        result = analytical_equilibrium(cell, mu_val, eps_val, mean_m_alpha;
                                        delta_nu=delta_nu, d_alpha=d_alpha, d_i=d_i,
                                        include_predators=include_predators,
                                        artificial_pi=artificial_pi, pi_size=pi_size,
                                        H_init=nothing, P_init=nothing)
        if result === nothing
            @warn "Skipping cell $cell due to setup error."
            continue
        end
        # Unpack equilibrium and parameters.
        A_eq = result.equilibrium
        A_p  = result.parameters
        J    = result.Jacobian
        
        # Compute eigenvalues and stability.
        eigvals = eigen(J).values
        max_eig = maximum(real.(eigvals))
        is_stable = all(real.(eigvals) .< 0)
        
        # Compute sensitivity metrics and species traits.
        metrics = compute_sensitivity_metrics(A_eq, A_p; perturbation=perturbation, tspan=tspan)
        n = length(A_eq.u0)
        
        # For each species, determine its type and record the data.
        for i in 1:n
            species_type = (i <= A_p.S) ? "herbivore" : "predator"
            push!(df, (cell, i, species_type, metrics.sensitivity[i],
                       metrics.biomass[i], metrics.degree[i], max_eig, is_stable))
        end
    end
    return df
end

# --- Run pipeline over multiple cells ---
# Define a range of cell indices (adjust as needed for your dataset)
cells_to_process = 1:1  # for example, cells 1 to 10

df_results = process_cells(cells_to_process; mu_val=0.5, eps_val=1.0, mean_m_alpha=0.1,
                             delta_nu=0.05, d_alpha=1.0, d_i=1.0,
                             include_predators=true, artificial_pi=false, pi_size=1.0,
                             perturbation=0.01, tspan=(0.0,50.0))

println("Collected results across cells:")
println(df_results)

# --- Optionally, save the DataFrame to CSV ---
# CSV.write("stability_sensitivity_results.csv", df_results)

# --- Evaluate relationships across cells ---
# Example: Plot sensitivity vs equilibrium biomass for all species.
begin
    fig_overall = Figure(resolution = (600, 500))
    ax_overall = Axis(fig_overall[1,1], xlabel = "Equilibrium Biomass", ylabel = "Sensitivity",
                    title = "Sensitivity vs Equilibrium Biomass across cells")
    MK.scatter!(ax_overall, df_results.degree, df_results.sensitivity, markersize = 8, 
            #  color = df_results.species_type .== "herbivore" ? :blue : :red, 
            label = "Species")
    # axislegend(ax_overall)
    display(fig_overall)
end
# You can similarly plot sensitivity vs degree or run regressions/correlations on df_results.
