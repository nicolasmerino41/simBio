# Function to compute the full removal impact Δϕᵢ
function compute_ecosystem_function_impact(A_eq, A_p; removal_fraction=1.0, tspan=(0.0,50.0), callbacks=true)
    # Get the equilibrium state vector and its total biomass (ϕ)
    u0 = A_eq.u0
    baseline_phi = sum(u0)
    n = length(u0)
    delta_phi = zeros(n)
    
    # For each species, remove it (or reduce its biomass by removal_fraction)
    for i in 1:n
        u0_removed = copy(u0)
        u0_removed[i] *= (1.0 - removal_fraction)  # if removal_fraction = 1.0 then species i is removed
        
        # Solve the system with the modified initial condition
        prob = ODEProblem(omnivore_dynamics!, u0_removed, tspan, A_p)
        sol = if callbacks
                  solve(prob, Tsit5(); callback = cb_no_trigger, abstol=1e-8, reltol=1e-6)
              else
                  solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
              end
        final_state = sol(sol.t[end])
        phi_removed = sum(final_state)
        
        # Full impact: how much the ecosystem function (total biomass) dropped
        delta_phi[i] = baseline_phi - phi_removed
    end
    return delta_phi
end

# Function to partition the impact into direct and indirect contributions.
# Here we assume that the ecosystem function is total biomass so that the direct contribution is just the equilibrium biomass.
function compute_direct_indirect_impacts(A_eq, A_p; removal_fraction=1.0, tspan=(0.0,50.0), callbacks=true)
    u0 = A_eq.u0
    # The direct contribution of each species (its equilibrium biomass)
    direct = copy(u0)
    # Full impact from removal (computed as before)
    full_impact = compute_ecosystem_function_impact(A_eq, A_p; removal_fraction=removal_fraction, tspan=tspan, callbacks=callbacks)
    # The indirect contribution is the additional change beyond the direct loss.
    indirect = full_impact .- direct
    return (full=full_impact, direct=direct, indirect=indirect)
end

# Function to plot the relationship between sensitivity and the partitioned impacts.
function plot_sensitivity_vs_partitioned_impact(A_eq, A_p; removal_fraction=1.0, tspan=(0.0,50.0), callbacks=true)
    # Get sensitivity metrics using your existing function (max deviation in total biomass)
    metrics = compute_sensitivity_metrics(A_eq, A_p; perturbation=removal_fraction, tspan=tspan, callbacks=callbacks)
    
    # Compute the full, direct, and indirect impacts
    impacts = compute_direct_indirect_impacts(A_eq, A_p; removal_fraction=removal_fraction, tspan=tspan, callbacks=callbacks)
    
    # Get species names (assumed to be the concatenation of herbivore and predator names)
    sp_names = vcat(A_eq.herbivore_list, A_eq.predator_list)
    
    fig = Figure(resolution = (1000, 700))
    ax1 = Axis(fig[1,1],
               xlabel = "Sensitivity (Max Deviation in Total Biomass)",
               ylabel = "Direct Contribution",
               title = "Sensitivity vs. Direct Contribution", 
               xticklabelrotation = π/2.5)
    ax2 = Axis(fig[1,2],
               xlabel = "Sensitivity (Max Deviation in Total Biomass)",
               ylabel = "Indirect Contribution",
               title = "Sensitivity vs. Indirect Contribution",
               xticklabelrotation = π/2.5)
    ax3 = Axis(fig[2,1],
               xlabel = "Sensitivity (Max Deviation in Total Biomass)",
               ylabel = "Full Impact (Δϕ)",
               title = "Sensitivity vs. Full Impact",
               xticklabelrotation = π/2.5)
    
    # Color points by guild
    colors = [metrics.guild[i] == "Herbivore" ? :blue :
              metrics.guild[i] == "Omnivore"   ? :green :
              metrics.guild[i] == "Predator"   ? :red : :gray
              for i in 1:length(metrics.guild)]
    
    MK.scatter!(ax1, metrics.sensitivity, impacts.direct, markersize = 10, color = colors)
    MK.scatter!(ax2, metrics.sensitivity, impacts.indirect, markersize = 10, color = colors)
    MK.scatter!(ax3, metrics.sensitivity, impacts.full, markersize = 10, color = colors)
    
    display(fig)
    return fig
end

# --- Example Usage ---
# Assuming A_eq and A_p were obtained from your analytical_equilibrium function
fig = plot_sensitivity_vs_partitioned_impact(A_eq, A_p; removal_fraction=0.5, tspan=(0.0, 50.0), callbacks=false)

# --- Merged Function ---
function run_stability_and_sensitivity_analysis(cell;
    removal_fraction=1.0, tspan=(0.0,50.0), callbacks=true,
    mu_range=0.0:0.1:1.0, eps_range=0.0:0.01:1.0, m_alpha_range=0.0:0.1:1.0)
    
    stable_found = false
    stable_result = nothing

    # Search for a stable equilibrium configuration over the parameter ranges.
    Threads.@threads for eps_val in eps_range
        for mu_val in mu_range
            for m_alpha in m_alpha_range
                result = analytical_equilibrium(cell, mu_val, eps_val, m_alpha;
                    delta_nu = 0.05,
                    d_alpha = 1.0, d_i = 1.0,
                    include_predators = true,
                    include_omnivores = true,
                    sp_removed_name = nothing,
                    artificial_pi = true, pi_size = 10.0,
                    H_init = nothing, P_init = nothing,
                    nu_omni_proportion = 1.0, nu_b_proportion = 1.0,
                    r_omni_proportion = 1.0,
                    callbacks = callbacks, plot = false)
                
                eigvals = eigen(result.Jacobian).values
                if all(real.(eigvals) .< 0)
                    println("Stable equilibrium found for cell $cell with μ = $mu_val, ε = $eps_val, mₐ = $m_alpha")
                    stable_found = true
                    stable_result = result
                    break
                end
            end
            if stable_found
                break
            end
        end
        if stable_found
            break
        end
    end
    
    if !stable_found
        println("No stable equilibrium found for cell $cell")
        return nothing
    end
    
    # Extract equilibrium and parameter data from the stable result.
    A_eq = stable_result.equilibrium    # Contains H_star, P_star, u0, herbivore_list, predator_list
    A_p  = stable_result.parameters     # Contains S, R, r_i, K_i, mu, nu, interaction matrices, etc.
    
    # Call the plotting function that combines sensitivity and ecosystem function impact.
    fig = plot_sensitivity_vs_partitioned_impact(A_eq, A_p; removal_fraction=removal_fraction, tspan=tspan, callbacks=callbacks)
    
    return (stable_result = stable_result, fig = fig)
end

for cell in 1:1
    run_stability_and_sensitivity_analysis(
        cell;
        removal_fraction=0.1, tspan=(0.0,1000.0), callbacks=false,
        mu_range=0.0:0.1:1.0, eps_range=0.0:0.01:1.0, m_alpha_range=0.0:0.1:1.0
    )
end