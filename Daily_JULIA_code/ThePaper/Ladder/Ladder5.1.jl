df = new_results10
# define the steps in order
step_keys = ["Full"; ["S$(i)" for i in 1:15]...]  # adjust the 15 if you have more/less steps

# the metrics you want to plot
metrics = [
    :return_time, 
    # :overshoot, 
    # :ire, 
    # :compound_error,
    :after_persistence
]

# for each metric, build a figure
for metric in metrics
    n = length(step_keys[1:end])
    fig = Figure(; size = (1000, 1000), fontsize=12)
    for i in 1:div(n, 4)
        
        for j in 1:4
            s = step_keys[j + (i-1)*4]
            # Create a colormap (e.g., viridis) based on number of consumers
            cmap = cgrad(:viridis, length(unique(df.connectance)))
            color_indices = [findfirst(==(val), sort(unique(df.connectance))) for val in df.connectance]
            color_vals = cmap[color_indices]
            ax = Axis(
                fig[i, j]; title = s,
                xlabel = "Degree CV",
                ylabel = string(metric),
                # title  = "Step $s"
            )
            # scatter degree_cv vs metric_s
            scatter!(
                ax, 
                df.degree_cv, 
                df[!, Symbol("$(metric)_$s")],
                markersize = 6,
                color = df.connectance,
                colormap = :viridis,
                colorrange = (minimum(df.connectance), maximum(df.connectance))
            )
            # add a little 1:1 reference line?
            # you could also compute & annotate a correlation here
        end
    end
    # fig.layoutgap = (10,10)
    display(fig)
end
