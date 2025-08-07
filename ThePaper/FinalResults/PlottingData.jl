G = R
# G = filter(row -> row.delta == -5.0, G)
step_keys = ["_full","_S1","_S2","_S3","_S5"]
res_cols = Symbol.("resilience" .* step_keys)
G = filter(row -> all(row[c] < 0 for c in res_cols), G)
println("subset size: ", nrow(G))
min_d_cols = Symbol.("sigma_over_min_d" .* step_keys)
G = filter(row -> all(row[c] < 10.0 for c in min_d_cols), G)
 
G = filter(row -> all(x -> !(x isa AbstractFloat) || (!isnan(x) && !isinf(x)), row), G)
sl_cols = Symbol.("SL" .* step_keys)
G = filter(row -> all(all(x -> x < 1000 && x > 0.0, row[c]) for c in sl_cols), G)
################### FOR SCALAR COMPARISONS ###################
# To show R² to 1:1 line  
for x in [:ER]
    plot_scalar_correlations_glv(
        G;
        scenarios = [:ER, :PL, :MOD],
        fit_to_1_1_line=true,
        save_plot = true,
        resolution = (1100, 950)
    )
end

T = filter(row -> row.scen == :PL, G)
E = filter(row -> row.scen == :ER, G)

G_low_Beq_cv = filter(row -> row.Beq_cv < 0.6, G)
G_high_Beq_cv = filter(row -> row.Beq_cv > 0.6, G)

begin
    want_fit = true
    colorBy = :conn
    scen = [:ER]
    plot_vector_correlations_glv(
        G;
        variable=:SL, color_by=colorBy,
        fit_to_1_1_line=want_fit,
        scenarios=scen,
        save_plot = true,
        resolution = (1000, 450),
        pixels_per_unit = 6.0
    )
end