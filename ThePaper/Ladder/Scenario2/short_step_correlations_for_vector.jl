function short_step_correlations_vectors(
    df::DataFrame,
    varx::Symbol,
    vary::Symbol;
    compare_to_full = true,
    color_by::Symbol = :conn,
    remove_unstable::Bool = false
)
    # step identifiers & titles
    step_keys  = ["S1","S2","S3","S4","S5","S6","S7","S8"]
    step_names = [
      "Full Model", "Global A (Global ϵ)", "Global AE",
      "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
      "Global A (Global ϵ) Mean B", "Global AE Mean B"
    ]

    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        @info "Filtered unstable: now $(nrow(df)) rows"
    end

    # build column symbols
    full_xcol   = Symbol(string(varx)*"_full")
    full_ycol   = Symbol(string(vary)*"_full")
    step_xcols  = Symbol.(string(varx) .* "_" .* step_keys)
    step_ycols  = Symbol.(string(vary) .* "_" .* step_keys)

    # color scale
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    ns   = length(step_keys)
    cols = 3

    fig = Figure(; size = (960, 640))
    Label(fig[0, 1:cols], uppercase("$(varx) vs $(vary) correlations"), fontsize=18)

    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1

        ax = Axis(fig[r, c];
            title  = step_names[idx],
            xlabel = string(compare_to_full ? full_xcol : step_xcols[idx]),
            ylabel = string(step_ycols[idx]),
        )

        xs = Float64[]
        ys = Float64[]
        cs = Float64[]

        for i in 1:nrow(df)
            xvec = compare_to_full ? df[i, full_xcol] : df[i, step_xcols[idx]]
            yvec = df[i, step_ycols[idx]]
            # @assert length(xvec) == length(yvec) "length mismatch in row $i, step $idx"
            for j in eachindex(xvec)
                push!(xs, xvec[j])
                push!(ys, yvec[j])
                push!(cs, color_vals[i])
            end
        end

        scatter!(ax, xs, ys;
            colormap   = :inferno,
            color      = cs,
            colorrange = (cmin, cmax),
            markersize = 6,
            alpha      = 0.7
        )

        mn = min(minimum(xs), minimum(ys))
        mx = max(maximum(xs), maximum(ys))
        lines!(ax, [mn, mx], [mn, mx], color = :black, linestyle = :dash)

        r_val = cor(xs, ys)
        text!(ax, "r=$(round(r_val, digits=3))",
            position = (mx, mn),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end

    display(fig)
end

#TODO You must load "exploring_min_extinction_100000.jls" once done from drago and run the following plots

short_step_correlations_vectors(
    df, :tau, :rt_press_vector;
    compare_to_full=false,
    color_by=:conn,
    remove_unstable=true
)

short_step_correlations_vectors(
    df, :tau, :rt_pulse_vector;
    compare_to_full=false,
    color_by=:conn,
    remove_unstable=false
)

short_step_correlations_vectors(
    df, :rt_pulse_vector, :rt_pulse_vector;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=true
)



