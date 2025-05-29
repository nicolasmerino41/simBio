step_keys = [
        "S1","S2","S3","S4","S5","S6", "S7", "S8"]
        # "S7"
    # ]
step_names = [
    "Full Model", "Global A (Global ϵ)", " Global AE",
    "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
    "Global A (Global ϵ) Mean B", "Global AE Mean B"
]

res_cols = Symbol.("resilience_" .* step_keys)
df = filter(row -> all(row[c] < 0 for c in res_cols), df)
println("subset size: ", nrow(df))
step_numbers = [1,2,3,4,5,6,7,8] # Numeric steps
n = nrow(df)

long_jdiff = DataFrame(
    Step = repeat(step_numbers, inner=n),
    J_diff = vcat([df[!, Symbol("J_diff_S$(s)")] for s in step_numbers]...)
)
begin
    fig = Figure(size=(800, 650))
    ax = Axis(fig[1,1], ylabel="Jacobian diff (Frobenius norm) logScale", xlabel="Step", title="Jacobian diff per step", yscale=log10,
        titlealign = :center, ylabelsize = 18, xlabelsize = 18, xticklabelrotation = π/4,
        xticks = (step_numbers, step_names), )
    boxplot!(ax, long_jdiff.Step, long_jdiff.J_diff)
    display(fig)    
end

begin
    fig = Figure(size=(800, 650))
    ax = Axis(fig[1,1], ylabel="Jacobian diff (Frobenius norm)", xlabel="Step", title="Jacobian diff per step",
        titlealign = :center, ylabelsize = 18, xlabelsize = 18, xticklabelrotation = π/4,
        xticks = (step_numbers, step_names), )
    boxplot!(ax, long_jdiff.Step, long_jdiff.J_diff)
    display(fig)    
end

