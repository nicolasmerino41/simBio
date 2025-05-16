using CairoMakie, Graphs
for i in 1:10
    row = A[i, :]  # select the row/community you want

    A_matrix = row.p_final[8]
    adjacency = abs.(A_matrix) .> 0
    g = SimpleGraph(adjacency)
    degs = degree(g)
    rt_vec = row.rt_press_full_vector

    R = 30  # number of resources
    C = 20  # number of consumers (assumed, adjust if needed)

    # Indices for resources and consumers
    resource_idx = 1:R
    consumer_idx = R+1:R+C

    begin
        fig = Figure(; size=(700, 400))
        ax = Axis(fig[1, 1];
            xlabel = "Degree (number of links)",
            ylabel = "Return Time (rt_pulse_full_vector)",
            title = "Degree vs Return Time (Resources vs Consumers)"
        )
        scatter!(ax, degs[resource_idx], rt_vec[resource_idx], color=:dodgerblue, label="Resources", markersize=11, alpha=0.8)
        scatter!(ax, degs[consumer_idx], rt_vec[consumer_idx], color=:crimson, label="Consumers", markersize=11, alpha=0.8)
        axislegend(ax, position=:rb)
        display(fig)
    end
end