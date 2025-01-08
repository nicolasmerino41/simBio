include("DA_birmmals_with_pi.jl")
include("Functions/generate_competition_matrix.jl")
include("Functions/species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
##############################################################################################
# We apply the dynamics using the species, variables and parameters from a chosen cell.
##############################################################################################
total_biomass_df = DataFrame(herbivore_biomass = Float64[], predator_biomass = Float64[], total_biomass = Float64[], i_star = Float64[])
i_star_idx = 4 
begin

    ## NOTE, COMMENT THE include(callbacks.jl) IF YOU'VE ALREADY RUN IT ONCE FOR THAT CELL
    extinction_trigger = false
    # Choose a cell
    aco = 2500
    i, j = idx[aco][1], idx[aco][2]

    # Extract the NPP from the cell
    NPP = 1000.0 #Float64(npp_DA[i, j])

    H0_vector = Vector{Float64}(H0_DA[i, j].a)
    
    S, R, species_names, herbivore_list, predator_list, H_i0, m_i, p_vec, x_final, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha = setup_community_from_cell(
        i, j;
        NPP = NPP,
        M_mean = 0.1,
        mu = 0.5,
        symmetrical_competition = false,
        mean_m_alpha = 0.1,
        epsilon_val = 0.01,
        mu_predation = 0.01,
        iberian_interact_NA=iberian_interact_NA,
        species_dict=species_dict,
        m_standard_deviation = 0.0,
        h_standard_deviation = 0.0,
        artificial_pi = false,
        real_H0 = false,
        H0_vector = H0_vector
    )

    # include("Callbacks.jl")
    
    params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
    H_init = H_i0
    P_init = H_i0[1:R] ./ 10 # fill(10.0, R)
    u0 = vcat(H_init, P_init)
    tspan = (0.0, 500.0)
    EXTINCTION_THRESHOLD = 1e-6
    T_ext = 100.0
    i_star = i_star_idx

    # Decide which callbacks to use:
    if extinction_trigger
        cb = cb_trigger
    else
        # Just the threshold-based extinction
        cb = cb_no_trigger
    end
    
    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    times = sol.t
    H_data = sol[1:S, :]
    P_data = sol[S+1:S+R, :]

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time", ylabel="Density", title="Dynamics")
    for x in 1:S
        lines!(ax, times, H_data[x, :], label="H$x")
    end
    for α in 1:R
        lines!(ax, times, P_data[α, :], label="P$α", linestyle=:dash)
    end
    # axislegend(ax; position=:rt)
    display(fig)

    # Summary stats
    herb_survivors = count(H_data[:, end] .> EXTINCTION_THRESHOLD)
    pred_survivors = count(P_data[:, end] .> EXTINCTION_THRESHOLD)

    println("Herbivores survived: $herb_survivors/$S")
    println("Predators survived: $pred_survivors/$R")

    H_biomass = sum(H_data[:, end][H_data[:, end] .> EXTINCTION_THRESHOLD])
    P_biomass = sum(P_data[:, end][P_data[:, end] .> EXTINCTION_THRESHOLD])
    println("Herbivore biomass at end: ", H_biomass)
    println("Predator biomass at end: ", P_biomass)
    println("Total biomass: ", H_biomass + P_biomass)
    println("herb_pred_ratio: ", P_biomass / H_biomass)
    println("NPP == ∑g_iH_i? ", isapprox(NPP, sum(g_i .* H_data[:, end]), atol=50.0))
    println("∑g_iH_i = ", sum(g_i .* H_data[:, end]))
    println("NPP = ", NPP)

    total_biomass_df = vcat(total_biomass_df, DataFrame(herbivore_biomass = H_biomass, predator_biomass = P_biomass, total_biomass = H_biomass + P_biomass, i_star = i_star))
end

###############################################################################
# 1) Build a DataFrame of final biomasses, and convert to abundances
###############################################################################
# We'll collect final data in a DataFrame with columns:
#   :species, :type, :final_biomass, :bodyMass, :final_abundance

df_final = DataFrame(
    species         = String[],
    type            = String[],   # "herb" or "pred"
    final_biomass   = Float64[],
    bodyMass        = Float64[],
    final_abundance = Float64[]
)

# A) For herbivores, final biomass is H_data[i, end].
for (idx, sp_name) in enumerate(herbivore_list)
    final_bio = H_data[idx, end]
    # If final_bio <= EXTINCTION_THRESHOLD, species effectively extinct => skip or keep as zero
    if final_bio <= EXTINCTION_THRESHOLD
        continue
    end

    # Find row in birmmals_biomass_fixed that matches species name
    row_idx = findfirst(row -> row[:species] == sp_name, eachrow(birmmals_biomass_fixed))
    if row_idx === nothing
        @warn "Species $sp_name not found in birmmals_biomass_fixed; skipping."
        continue
    end

    # Extract bodyMass
    bm_val = birmmals_biomass_fixed[row_idx, :bodyMass]
    # Convert biomass -> abundance
    final_abund = final_bio / bm_val

    push!(df_final, (
        species         = sp_name,
        type            = "herb",
        final_biomass   = final_bio,
        bodyMass        = bm_val,
        final_abundance = final_abund
    ))
end

# B) For predators, final biomass is P_data[α, end].
for (pidx, sp_name) in enumerate(predator_list)
    final_bio = P_data[pidx, end]
    if final_bio <= EXTINCTION_THRESHOLD
        continue
    end

    row_idx = findfirst(row -> row[:species] == sp_name, eachrow(birmmals_biomass_fixed))
    if row_idx === nothing
        @warn "Species $sp_name not found in birmmals_biomass_fixed; skipping."
        continue
    end

    bm_val = birmmals_biomass_fixed[row_idx, :bodyMass]
    final_abund = final_bio / bm_val

    push!(df_final, (
        species         = sp_name,
        type            = "pred",
        final_biomass   = final_bio,
        bodyMass        = bm_val,
        final_abundance = final_abund
    ))
end

###############################################################################
# 2) Sort in descending order of final_abundance
###############################################################################
df_final_sorted = sort(df_final, :final_abundance, rev=true)
df_final_sorted_bio = sort(df_final, :final_biomass, rev=true)
###############################################################################
# 3) Plot a bar chart of final_abundance by species
###############################################################################
begin
    log = false

    # Sort species by abundance in descending order
    sorted_indices = sortperm(df_final_sorted.final_abundance, rev=true)
    sorted_species = df_final_sorted.species[sorted_indices]
    sorted_abundance = df_final_sorted.final_abundance[sorted_indices]
    sorted_biomass = df_final_sorted.final_biomass[sorted_indices]

    fig_bar = Figure(resolution=(1500, 800))

    # Plot abundance
    ax_abundance = Axis(fig_bar[1, 1], xlabel="Species", ylabel="Abundance",
                        title="Species Abundance (last time step)",
                        yscale=log ? log10 : identity)
    species_indices = 1:length(sorted_species)  # Numeric indices for sorted species
    barplot!(ax_abundance, species_indices, sorted_abundance)
    ax_abundance.xticks = (species_indices, sorted_species)
    ax_abundance.xticklabelrotation = π / 2  # Rotate x-tick labels 90 degrees

    # Plot biomass
    ax_biomass = Axis(fig_bar[1, 2], xlabel="Species", ylabel="Biomass",
                      title="Species Biomass (last time step)",
                      yscale=log ? log10 : identity)
    barplot!(ax_biomass, species_indices, sorted_biomass)
    ax_biomass.xticks = (species_indices, sorted_species)
    ax_biomass.xticklabelrotation = π / 2  # Rotate x-tick labels 90 degrees

    display(fig_bar)
end



