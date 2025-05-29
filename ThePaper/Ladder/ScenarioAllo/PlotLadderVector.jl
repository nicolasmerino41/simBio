function allo_step_correlations_vectors(
    df::DataFrame,
    varx::Symbol,
    vary::Symbol;
    compare_to_full::Bool    = true,
    color_by::Symbol         = :conn,
    remove_unstable::Bool    = false,
    remove_zeros::Bool       = false,
    equal_axes::Bool         = true
)
    step_keys  = ["S1","S2"]
    step_names = ["Full Model", "Global A (Global ϵ)"]

    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        @info "Filtered unstable: now $(nrow(df)) rows"
    end

    # build symbols for full and step‐columns
    full_xcol   = Symbol(string(varx) * "_full")
    full_ycol   = Symbol(string(vary) * "_full")    # <-- singular
    step_xcols  = Symbol.(string(varx) .* "_" .* step_keys)
    step_ycols  = Symbol.(string(vary) .* "_" .* step_keys)

    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    fig = Figure(size=(960, 480))
    Label(fig[0, 1:2], uppercase("$(varx) vs $(vary) correlations"), fontsize=18)

    for idx in 1:2
        ax = Axis(fig[1, idx];
            title  = step_names[idx],
            xlabel = string(compare_to_full ? full_xcol   : step_xcols[idx]),
            ylabel = string(step_ycols[idx]),
        )

        xs = Float64[]
        ys = Float64[]
        cs = Float64[]

        for i in 1:nrow(df)
            xvec = compare_to_full ? df[i, full_xcol]   : df[i, step_xcols[idx]]
            yvec = df[i, step_ycols[idx]][1:length(xvec)]
            @assert length(xvec) == length(yvec) "length mismatch at row $i, step $idx, 
                xvec: $(length(xvec)), yvec: $(length(yvec))"

            for j in eachindex(xvec)
                xj, yj = xvec[j], yvec[j]
                if remove_zeros && (xj == 0.0 || yj == 0.0)
                    continue
                end
                push!(xs, xj)
                push!(ys, yj)
                push!(cs, color_vals[i])
            end
        end

        scatter!(ax, xs, ys;
            colormap   = :viridis,
            color      = cs,
            colorrange = (cmin, cmax),
            markersize = 6,
            alpha      = 0.7
        )

        xmin, xmax = minimum(xs), maximum(xs)
        ymin, ymax = minimum(ys), maximum(ys)

        if equal_axes
            mn = min(xmin, ymin)
            mx = max(xmax, ymax)
            xlims!(ax, mn, mx)
            ylims!(ax, mn, mx)
            lines!(ax, [mn, mx], [mn, mx], color=:black, linestyle=:dash)
        else
            xlims!(ax, xmin, xmax)
            ylims!(ax, ymin, ymax)
            lines!(ax, [xmin, xmax], [xmin, xmax], color=:black, linestyle=:dash)
        end

        r_val = cor(xs, ys)
        text!(ax, "r=$(round(r_val, digits=3))",
            position = (xmax, ymin),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end

    display(fig)
    return fig
end

#TODO once drago finishes ShortAlloLadder_withMinExtinction_768.jls re-run this again
# df = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder_withMinExtinction_48.jls")
df = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder_withMinExtinction_96.jls")
df = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder_withMinExtinction_96_20species.jls")
df.model_type_n = map(df.model_type) do mt
    if mt == :cascade
        1
    elseif mt == :niche
        2
    else
        error("Unknown model_type: $mt")
    end
end
# df = A

allo_step_correlations_vectors(
    df, :tau, :rt_pulse_vector;
    color_by=:conn, 
    compare_to_full=true,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=false
)

allo_step_correlations_vectors(
    df, :min_delta_K, :min_delta_K;
    color_by=:conn, 
    remove_unstable=false,
    compare_to_full=true,
    remove_zeros=false,
    equal_axes=true
)




U = df[:, [:min_delta_K_full, :min_delta_K_S1, :min_delta_K_S2]]
