#### THIS SCRIPT TRIES TO GATHER ALL THE INFO, METRICS, PLOTS, I HAVE OBTAINED TIL NOW ####
#### IT DOES NOT CREATE ANYTHING BUT OUTPUTS, HENCE, EVERY FUNCTION NEEDS TO HAVE BEEN ####
#### PREVIOUSLY CALLED BEFORE RUNNING THIS SCRIPT ####

include("CVIvsNRI.jl")
include("CSM.jl")
include("CommunityNPPsaturation.jl")

df_to_work_with = all_results_list_even_pi
# --------------------------
# 1) SPECIES_COUNT_DF
# --------------------------
species_count_df = create_species_count_df(df_to_work_with)
display(species_count_df)

# --------------------------
# 2) PLOTTING SPECIES METRICS
# --------------------------
# This requires new_all_results_list from Applying_results.jl
# And plot_species_metrics() comes from Influential_species.jl
selected_metric = :indegree # Change to :indegree, :outdegree, :total_degree, :closeness, or :clustering
plot_species_metrics(species_count_df, new_all_results_list, selected_metric)

# --------------------------
# 3) SPECIES EFFECT ON ECOSYSTEM FUNCTIONING
# --------------------------
# The SEEF function is defined in Functions/SEEF_function.jl
see_not_even = SEEF(all_results_list)
see_even     = SEEF(all_results_list_even_pi)

see_to_plot  = see_even
# The plot_species_effects function comes from SpeciesEffectOnEcosystemFunctioning.jl
plot_species_effects(
    see_to_plot;
    log = false, 
    avg_eff_or_avg_eff_stand = true # true for average_effect, false for average_effect_standardized
)

# --------------------------
# 4) PLOTING AVERAGE EFFECT OF EACH SPECIES VS THEIR METRICS
# --------------------------
# compute_average_species_metrics comes from Functions/Computing_metrics.jl
cell_range = 1:5950
avg_species_metrics = compute_average_species_metrics(cell_range)

# The plot_average_effect_vs_metrics function comes from SpeciesEffectOnEcosystemFunctioning.jl
plot_average_effect_vs_metrics(
    see_even;
    avg_species_metrics = avg_species_metrics 
)

# --------------------------
# 5) ORIGINAL SENSITIVITY FROM GlobalMetricsRelationship.jl
# --------------------------
cell_sensitivity_df = measure_cell_sensitivity(
    df_to_work_with;
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
    capped = false, cap_val = 1.18
)
grid_CVI_not_even_pi = map_cell_metric(
    new_cell_stability_df_not_even_pi, :CVI; 
    title = "Cell Vulnerability Index (CVI) with not even pi",
    capped = true, cap_val = 1.18
)
grid_NRI = map_cell_metric(
    new_cell_stability_df_even_pi, :NRI; 
    title = "Network Robustness Index (NRI)"
) # THIS ONE IS QUITE USELESS BECAUSE YOU CAN GET ALL NRI FROM THE FOLLOWING

# ----------------------------
# 9) PLOT THE FULL NRI GRID
# ----------------------------
# This function is defined in Functions/Computing_metrics.jl
compute_and_map_NRI(;
    plot = true, 
    title = "NRI",
    standardise_by_NPP = false
)

# ----------------------------
# 10) Scatter plot correlating CVI and NRI across cells
# ----------------------------
scatter_CVI_NRI(
    new_cell_stability_df_even_pi;
    info = false
)

# ----------------------------
# 11) COMPUTE CSM AND PLOT IT (From CSM.jl)
# ----------------------------
csm_df = compute_CSM(df_to_work_with)
csm_grid = map_CSM(
    csm_df; 
    plot = true, palette = custom_palette,
    resolution = (600, 600), title = "CSM"
)

# ----------------------------
# 12) COMMUNITY NPP SATURATION (From CommunityNPPsaturation.jl)
# ----------------------------
npp_saturation_df = CommunityNPPsaturation(
    Big_P_even_pi_maximised; 
    scatter = true, 
    map = true, NPP_aside = true,
    resolution_scatter = (600,600),
    resolution_map = (1000,400)
)
