using CairoMakie

# Replace this with your actual DataFrame variable
# df = DataFrame(...) 
dff = copy(df)  # Assuming df is your DataFrame
res_cols = Symbol.("resilience_" .* step_keys)
dff = filter(row -> all(row[c] < 0 for c in res_cols), dff)
res_cols = Symbol.("after_pulse_" .* step_keys)
dff = filter(row -> all(row[c] .> 0.9 for c in res_cols), dff)
println("subset size: ", nrow(dff))
# Prepare arrays for plotting
rt_resources = Float64[]
rt_consumers = Float64[]
x_resources = Float64[]
x_consumers = Float64[]

for i in 1:10
    rt_resources = Float64[]
    rt_consumers = Float64[]
    x_resources = Float64[]
    x_consumers = Float64[]
    for row in eachrow(dff)[i:i]  # Assuming dff is a DataFrame
        rt_vec = row[:rt_pulse_full_vector]
        k_xi   = row[:K_Xi_full]
        R_eq   = row[:R_eq]
        C_eq   = row[:C_eq]
        # m_val  = row[:m_val]
        # g_val  = row[:g_val]
        m_val  = row[:m_cons]
        g_val  = row[:r_res]

        # # Resources: positions 1:30
        # for i in 1:30
        #     Ki = k_xi[i]
        #     rBi = R_eq[i]
        #     xval = g_val[i] * rBi / Ki
        #     push!(x_resources, xval)
        #     push!(rt_resources, rt_vec[i])
        # end

        # Consumers: positions 31:50
        for i in 1:20
            xi = k_xi[30+i]
            Bi = C_eq[i]
            xval = m_val[i] * Bi / xi
            push!(x_consumers, xval)
            push!(rt_consumers, rt_vec[30+i])
        end
    end

    begin
        fig = Figure(; size = (800, 500))
        ax = Axis(fig[1, 1];
            xlabel = "Self-Regulation Loss (rB/K for resources, mB/xi for consumers)",
            ylabel = "Return Time",
            title = "Self-Regulation Loss vs Return Time",
            titlealign = :center,
            ylabelsize = 18,
            xlabelsize = 18,
        )
        scatter!(ax, x_resources, rt_resources, color = :blue, markersize = 7, label = "Resources (rB/K)", alpha=0.6)
        scatter!(ax, x_consumers, rt_consumers, color = :red, markersize = 7, label = "Consumers (mB/xi)", alpha=0.6)

        axislegend(ax; position = :rt, framevisible=false)
        display(fig)
    end
end
