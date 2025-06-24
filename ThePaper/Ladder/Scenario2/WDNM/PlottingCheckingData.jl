include("FunctionsPlottingCheckingData.jl")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking60000.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography60000.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography60000PL.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography400000PL.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography100000MOD.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography300000PL.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_recalculating_demography10000ER.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/NF_checking50000ER.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_yes_recalculating_50000ER.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_not_recalculating_50000ER.jls")
R = deserialize("ThePaper/Ladder/Outputs/checking/checking_changing_groups_100000ER.jls")

################ CLEANING ################
desired = [
  :conn, :IS, :scen, :delta, :epsi, :m_val, :g_val, :ite, :pex, :p_min_deg, :mod_gamma,
  :resilience_full, :reactivity_full, :collectivity_full, :tau_full, :mean_tau_full, :sigma_over_min_d_full, :SL_full, :mean_SL_full, :inverse_tau_full, :mean_inverse_tau_full, :analytical_rmed_full, :ssp_analytical_rmed_full,
  :resilience_S1, :reactivity_S1, :collectivity_S1, :tau_S1, :mean_tau_S1, :sigma_over_min_d_S1, :SL_S1, :mean_SL_S1, :inverse_tau_S1, :mean_inverse_tau_S1, :analytical_rmed_S1, :ssp_analytical_rmed_S1,
  :resilience_S2, :reactivity_S2, :collectivity_S2, :tau_S2, :mean_tau_S2, :sigma_over_min_d_S2, :SL_S2, :mean_SL_S2, :inverse_tau_S2, :mean_inverse_tau_S2, :analytical_rmed_S2, :ssp_analytical_rmed_S2,
  :resilience_S3, :reactivity_S3, :collectivity_S3, :tau_S3, :mean_tau_S3, :sigma_over_min_d_S3, :SL_S3, :mean_SL_S3, :inverse_tau_S3, :mean_inverse_tau_S3, :analytical_rmed_S3, :ssp_analytical_rmed_S3,
  :resilience_S4, :reactivity_S4, :collectivity_S4, :tau_S4, :mean_tau_S4, :sigma_over_min_d_S4, :SL_S4, :mean_SL_S4, :inverse_tau_S4, :mean_inverse_tau_S4, :analytical_rmed_S4, :ssp_analytical_rmed_S4
]

G = R[!, desired]
G = R

step_keys = ["_full","_S1","_S2","_S3","_S4","_S5"]
res_cols = Symbol.("resilience" .* step_keys)
G = filter(row -> all(row[c] < 0 for c in res_cols), G)
println("subset size: ", nrow(G))
min_d_cols = Symbol.("sigma_over_min_d" .* step_keys)
G = filter(row -> all(row[c] < 50.0 for c in min_d_cols), G)

G = filter(row -> all(x -> !(x isa AbstractFloat) || (!isnan(x) && !isinf(x)), row), G)

################### FOR SCALAR COMPARISONS ###################
# To show R² to 1:1 line
plot_scalar_correlations(
    G;
    scenarios = [:ER],
    fit_to_1_1_line=false
)

# To show Pearson correlation r
plot_scalar_correlations(
    G;
    scenarios = [:ER],
    fit_to_1_1_line=false
)

################### for vector correlations ###################
begin
    want_fit = true
    colorBy = :conn
    scen = [:ER]
    plot_vector_correlations(
        G;
        variable=:tau, color_by=colorBy,
        fit_to_1_1_line=want_fit,
        scenarios=scen
    )
    
    plot_vector_correlations(
        G;
        variable=:inverse_tau, color_by=colorBy,
        fit_to_1_1_line=want_fit,
        scenarios=scen
    )
    
    plot_vector_correlations(
        G;
        variable=:SL, color_by=colorBy,
        fit_to_1_1_line=want_fit,
        scenarios=scen
    )
    
    plot_vector_correlations(
        G;
        variable=:ssp_analytical_rmed, color_by=colorBy,
        fit_to_1_1_line=want_fit,
        scenarios=scen
    )
end