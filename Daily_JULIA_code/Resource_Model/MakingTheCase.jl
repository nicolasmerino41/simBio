#### THIS SCRIPT TRIES TO GATHER ALL THE INFO, METRICS, PLOTS, I HAVE OBTAINED TIL NOW ####
#### IT DOES NOT CREATE ANYTHING BUT OUTPUTS, HENCE, EVERY FUNCTION NEEDS TO HAVE BEEN ####
#### PREVIOUSLY CALLED BEFORE RUNNING THIS SCRIPT ####

include("CVIvsNRI.jl")
include("CSM.jl")
include("CommunityNPPsaturation.jl")
include("GlobalMetricsRelationship.jl")

df_to_work_with = all_results_list
# --------------------------
# 1) SPECIES_COUNT_DF
# --------------------------
species_count_df = create_species_count_df(
    df_to_work_with,
    only_most_influential = true
)

display(species_count_df)

# --------------------------
# 2) PLOTTING SPECIES METRICS
# --------------------------
# This requires new_all_results_list from Applying_results.jl
# And plot_species_metrics() comes from Influential_species.jl
selected_metric = :clustering # Change to :indegree, :outdegree, :total_degree, :closeness, or :clustering
plot_species_metrics(
    species_count_df,
    new_all_results_list,
    selected_metric::Symbol;
    only_most_influential = true, # If you set create_species_count_df() to only_most_influential = true
    count_or_stand_count = false, # If true, use count; if false, use stand_count
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    by_name_or_by_TL = true,  # if true, color by herb/pred; if false, color continuously by TL from TrophInd
    palette = custom_palette
)
for i in 1:20 
    println(sort(species_count_df, :stand_count, rev=true).species_name[i])
end
# --------------------------
# 3) SPECIES EFFECT ON ECOSYSTEM FUNCTIONING
# --------------------------
# The SEEF function is defined in Functions/SEEF_function.jl
see_not_even = SEEF(all_results_list)
see_even     = SEEF(all_results_list_even_pi)

see_to_plot  = see_not_even
for i in 1:20
    println(sort(see_to_plot, :average_effect, rev=true).species_name[i])
end
# The plot_species_effects function comes from SpeciesEffectOnEcosystemFunctioning.jl
plot_species_effects(
    see_to_plot;
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    log = false,
    standardise_by_H0 = false,
    resolution = (1100, 600),
    by_name_or_by_TL = true,  # if true, assign colors by name; if false, use continuous TL from TrophInd
    palette = custom_palette
)

# --------------------------
# 4) PLOTING AVERAGE EFFECT OF EACH SPECIES VS THEIR METRICS
# --------------------------
# compute_average_species_metrics comes from Functions/Computing_metrics.jl
cell_range = 1:5950
# avg_species_metrics = compute_average_species_metrics(cell_range) # This takes a long time

# The plot_average_effect_vs_metrics function comes from SpeciesEffectOnEcosystemFunctioning.jl
plot_average_effect_vs_metrics(
    see_not_even;
    avg_species_metrics = avg_species_metrics,
    herbivore_names = herbivore_names,
    predator_names = predator_names,
    by_name_or_by_TL = true,  # if true, color by species type; if false, color continuously by TL from TrophInd
    palette = custom_palette
)

# --------------------------
# 4.1) PLOTING SPECIES'S EFFECT VS THEIR METRICS (CELL-SPECIFIC)
# --------------------------
# The plot_species_effect_vs_cell_metrics function comes from SpeciesEffectOnEcosystemFunctioning.jl
for i in 1:20
plot_species_effect_vs_cell_metrics(
    i
)
end
# map_plot(npp_DA; palette = custom_palette)
# --------------------------
# 5) ORIGINAL SENSITIVITY FROM GlobalMetricsRelationship.jl
# --------------------------
cell_sensitivity_df = measure_cell_sensitivity(
    all_results_list;
    capped = true, cap_val = 5.0
)

sensitivity_grid = map_cell_sensitivity(
    cell_sensitivity_df;
    disp = true # Do you want to see the plot or just return the grid?
)

# --------------------------
# 6) PLOTTING THE CORRELATION BETWEEN GLOBAL METRICS AND SENSITIVITY
# --------------------------
plot_global_metrics_vs_sensitivity(
    cell_sensitivity_df;
    save = false
)

# ----------------------------
# 7) Compute CVI and NRI for each cell (From CVIvsNRI.jl)
# ----------------------------
new_cell_stability_df_even_pi = new_measure_cell_stability(all_results_list_even_pi)
new_cell_stability_df_not_even_pi = new_measure_cell_stability(all_results_list)
# ----------------------------
# 8) Mapping functions to visualize cell-level metrics on the grid
# ----------------------------
grid_CVI_even_pi = map_cell_metric(
    new_cell_stability_df_even_pi, :CVI;
    title = "Cell Vulnerability Index (CVI) with even pi",
    standardize_by_NPP = false,
    # capped = true, cap_val = 1.18
)
grid_CVI_not_even_pi = map_cell_metric(
    new_cell_stability_df_not_even_pi, :CVI; 
    title = "Cell Vulnerability Index (CVI) with not even pi",
    capped = true, cap_val = 1.02
)
maximum(filter(!isnan, grid_CVI_not_even_pi))
grid_NRI = map_cell_metric(
    new_cell_stability_df_even_pi, :NRI; 
    title = "Network Robustness Index (NRI)"
) # THIS ONE IS QUITE USELESS BECAUSE YOU CAN GET ALL NRI FROM POINT 9

# ----------------------------
# 9) PLOT THE FULL NRI GRID
# ----------------------------
# This function is defined in Functions/Computing_metrics.jl
compute_and_map_NRI(;
    plot = true, 
    title = "NRI",
    standardise_by_NPP = false,
    resolution = (1000, 700)
)

# ----------------------------
# 10) Scatter plot correlating CVI and NRI across cells
# ----------------------------
scatter_CVI_NRI(
    new_cell_stability_df_not_even_pi;
    info = false
)

# ----------------------------
# 11) COMPUTE CSM AND PLOT IT (From CSM.jl)
# ----------------------------
csm_df = compute_CSM(all_results_list)
csm_grid = map_CSM(
    csm_df; 
    plot = true, palette = custom_palette,
    resolution = (600, 400), title = "CSM"
)

# ----------------------------
# 12) COMMUNITY NPP SATURATION (From CommunityNPPsaturation.jl)
# ----------------------------
npp_saturation_df, npp_saturation_grid = CommunityNPPsaturation(
    Big_P_results_maximised; 
    scatter = true,
    colour_scatter_by_richness_or_H0 = false,
    map = true,
    palette = custom_palette,
    resolution_scatter = (600,400),
    resolution_map = (1000,400),
    scatter_title = "NPP vs. Total Biomass",
    map_title = "Residuals (Observed - Predicted Biomass)",
    NPP_aside = true,
    richness_aside = false,
    evaluate_richness = true
)

# ----------------------------
# 13) THE EFFECT OF RICHNESS
# ----------------------------
DA_richness_birmmals = deserialize("Objects/DA_richness_birmmals.jls")
npp_saturation_df, npp_saturation_grid = CommunityNPPsaturation(
    Big_P_results_maximised;
    scatter = true,
    colour_scatter_by_richness_or_H0 = false, # true = richness, false = H0
    map = true,
    palette = custom_palette,
    resolution_scatter = (600,400),
    resolution_map = (1000,400),
    scatter_title = "NPP vs. Total Biomass",
    map_title = "Residuals (Observed - Predicted Biomass)",
    NPP_aside = false,
    richness_aside = true,
    evaluate_richness = true
)

# ----------------------------
# 14) PLOT CLUSTER MAP
# ----------------------------
# To run function 14 and 15 you need to first run CharacterisingCommunities.jl and obtain the cluster grid
number_clusters = 3
include("CharacterisingCommunities.jl")
plot_cluster_map(
    cluster_grid; 
    resolution=(600,400), 
    title="Community Clusters",
    colormap=custom_palette
)

# ----------------------------
# 15) PLOT CLUSTER SCATTER
# ----------------------------
plot_cluster_scatter(
    df_cells; 
    xvar = :over_under, yvar = :CSM, # You can choose :NPP, :CSM, :CVI, :NRI, :raw_sensitivity, :richness and, :over_under
    cluster_col = :cluster, 
    resolution=(600,600),
    colormap=:Set1, markersize=8
)

# ----------------------------
# 16) PCA EXPLAINS OVER-UNDER
# ----------------------------
pca_explain_over_under(
    df_cells,
    vars = [:raw_sensitivity, :CVI, :NRI, :CSM, :NPP, :richness, :over_under]
)

# ----------------------------
# 17) OVER-UNDER VS LUIS' DATA
# ----------------------------
begin
    fig = Figure(resolution = (1000, 600))
    ax = Axis(fig[1, 1])
    heatmap!(ax, npp_saturation_grid; interpolate = false, colormap = custom_palette)
    ax2 = Axis(fig[1, 2])
    MK.heatmap!(ax2, Matrix(DA_birmmals_with_pi_corrected); colormap = custom_palette)
    ax.yreversed = true
    ax2.yreversed = true
    display(fig)
end
