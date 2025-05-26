
using DataFrames

res_cols = Symbol.("resilience_" .* step_keys)
dff = filter(row -> all(row[c] < 0 for c in res_cols), dff)
println("subset size: ", nrow(dff))
step_numbers = [1,3,4,5,6,7,8] # Numeric steps
n = nrow(dff)

long_jdiff = DataFrame(
    Step = repeat(step_numbers, inner=n),
    J_diff = vcat([dff[!, Symbol("J_diff_S$(s)")] for s in step_numbers]...)
)

fig = Figure(size=(800, 400))
ax = Axis(fig[1,1], ylabel="Jacobian diff (Frobenius norm)", xlabel="Step", title="Jacobian diff per step")
boxplot!(ax, long_jdiff.Step, long_jdiff.J_diff)
display(fig)
