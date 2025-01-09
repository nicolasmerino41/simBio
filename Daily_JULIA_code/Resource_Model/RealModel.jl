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

# for i_star_idx in 1:10 
begin

    ## NOTE, COMMENT THE include(callbacks.jl) IF YOU'VE ALREADY RUN IT ONCE FOR THAT CELL
    extinction_trigger = true
    # Choose a cell
    aco = 417
    i, j = idx[aco][1], idx[aco][2]

    # Extract the NPP from the cell
    NPP = 1000.0 #Float64(npp_DA[i, j])

    H0_vector = Vector{Float64}(H0_DA[i, j].a)
    
    S, R, species_names, herbivore_list, predator_list, H_i0, m_i, p_vec, x_final, g_i, localHatH, G, M_modified, a_matrix, A, epsilon, m_alpha = setup_community_from_cell(
        i, j;
        NPP = NPP,
        M_mean = 0.1,
        mu = 0.5,
        symmetrical_competition = true,
        mean_m_alpha = 0.1,
        epsilon_val = 1.0,
        mu_predation = 0.01,
        iberian_interact_NA=iberian_interact_NA,
        species_dict=species_dict,
        m_standard_deviation = 0.0,
        h_standard_deviation = 0.0,
        artificial_pi = false,
        real_H0 = true,
        H0_vector = H0_vector
    )

    # include("Callbacks.jl")
    
    params = (S, R, H_i0, m_i, g_i, G, M_modified, a_matrix, A, epsilon, m_alpha)
    H_init = H_i0
    P_init = H_i0[1:R] ./ 10 # fill(10.0, R)
    u0 = vcat(H_init, P_init)
    tspan = (0.0, 500.0)
    EXTINCTION_THRESHOLD = 1e-6
    T_ext = 250.0
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
# end

begin
    # -------------------------------------------------------------------------
    # 0) Set parameters
    # -------------------------------------------------------------------------
    include_predators   = false  
    sort_by_abundance   = true  # Sort by abundance true or sort by biomass false
    log                 = false  # Choose log scale

    # -------------------------------------------------------------------------
    # 1) Build a DataFrame of final biomasses, and convert to abundances
    # -------------------------------------------------------------------------
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
            type            = "herb",
            final_biomass   = final_bio,
            bodyMass        = bm_val,
            final_abundance = final_abund
        ))
    end

    # B) For predators (only if include_predators is true)
    if include_predators
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
    end

    # -------------------------------------------------------------------------
    # 2) Split df_final into herbivores and predators
    # -------------------------------------------------------------------------
    herb_mask = df_final.type .== "herb"
    pred_mask = df_final.type .== "pred"  # will be empty if include_predators=false

    # -- Sort herbivores
    if sort_by_abundance
        sorted_herb_inds = sortperm(df_final[herb_mask, :final_abundance], rev=true)
    else
        sorted_herb_inds = sortperm(df_final[herb_mask, :final_biomass], rev=true)
    end

    herb_species_sorted    = df_final[herb_mask, :species][sorted_herb_inds]
    herb_abundance_sorted  = df_final[herb_mask, :final_abundance][sorted_herb_inds]
    herb_biomass_sorted    = df_final[herb_mask, :final_biomass][sorted_herb_inds]

    # -- Sort predators (only if include_predators)
    pred_species_sorted    = String[]
    pred_abundance_sorted  = Float64[]
    pred_biomass_sorted    = Float64[]
    if include_predators && !isempty(df_final[pred_mask, :species])
        if sort_by_abundance
            sorted_pred_inds = sortperm(df_final[pred_mask, :final_abundance], rev=true)
        else
            sorted_pred_inds = sortperm(df_final[pred_mask, :final_biomass], rev=true)
        end

        pred_species_sorted    = df_final[pred_mask, :species][sorted_pred_inds]
        pred_abundance_sorted  = df_final[pred_mask, :final_abundance][sorted_pred_inds]
        pred_biomass_sorted    = df_final[pred_mask, :final_biomass][sorted_pred_inds]
    end

    # -------------------------------------------------------------------------
    # 3) Create figure with either a 1×2 or 2×2 layout
    # -------------------------------------------------------------------------
    # If you prefer always a 2×2 figure (even if predators are not included),
    # you can omit this if-else and just do a 2×2. For demonstration, we make
    # the figure "shorter" when not including predators.
    # -------------------------------------------------------------------------
    if include_predators
        fig_bar = Figure(resolution=(1500, 1200))
    else
        fig_bar = Figure(resolution=(1500, 900))
    end

    # -------------------------------------------------------------------------
    # 3A) Top row: Herbivore Abundance and Biomass
    # -------------------------------------------------------------------------
    # Herbivore biomass
    ax_biomass_herb = Axis(
        fig_bar[1, 1],
        # xlabel="Herbivore Species",
        ylabel="Biomass",
        title="Herbivores: Biomass",
        yscale=log ? log10 : identity
    )
    herb_indices = 1:length(herb_species_sorted)
    barplot!(ax_biomass_herb, herb_indices, herb_biomass_sorted)
    ax_biomass_herb.xticks = (herb_indices, herb_species_sorted)
    ax_biomass_herb.xticklabelrotation = π/2.5
    
    # Herbivore abundance
    ax_abundance_herb = Axis(
        fig_bar[1, 2],
        # xlabel="Herbivore Species",
        ylabel="Abundance",
        title="Herbivores: Abundance",
        yscale=log ? log10 : identity
    )
    herb_indices = 1:length(herb_species_sorted)
    barplot!(ax_abundance_herb, herb_indices, herb_abundance_sorted)
    ax_abundance_herb.xticks = (herb_indices, herb_species_sorted)
    ax_abundance_herb.xticklabelrotation = π/2.5

    # -------------------------------------------------------------------------
    # 3B) Bottom row: Predator Abundance and Biomass (only if include_predators)
    # -------------------------------------------------------------------------
    if include_predators && !isempty(pred_species_sorted)
        ax_biomass_pred = Axis(
            fig_bar[2, 1],
            # xlabel="Predator Species",
            ylabel="Biomass",
            title="Predators: Biomass",
            yscale=log ? log10 : identity
        )
        barplot!(ax_biomass_pred, pred_indices, pred_biomass_sorted)
        ax_biomass_pred.xticks = (pred_indices, pred_species_sorted)
        ax_biomass_pred.xticklabelrotation = π/2.5

        ax_abundance_pred = Axis(
            fig_bar[2, 2],
            # xlabel="Predator Species",
            ylabel="Abundance",
            title="Predators: Abundance",
            yscale=log ? log10 : identity
        )
        pred_indices = 1:length(pred_species_sorted)
        barplot!(ax_abundance_pred, pred_indices, pred_abundance_sorted)
        ax_abundance_pred.xticks = (pred_indices, pred_species_sorted)
        ax_abundance_pred.xticklabelrotation = π/2.5
    end

    display(fig_bar)
end

H_i0
# Find the cell with highest richness excluding NaNs
filtered_richness = filter(!isnan, DA_richness[:])
max_richness = maximum(filtered_richness)
max_richness_idx = findfirst(isequal(max_richness), DA_richness)
max_richness_cell = CartesianIndex(max_richness_idx.I...)
max_richness_idx_in_idx = findfirst(isequal(max_richness_idx), idx)
println("Cell with highest richness excluding NaNs: ", max_richness_cell)

perro = localHatH
perro = sort(perro, rev = true)
barplot(perro)