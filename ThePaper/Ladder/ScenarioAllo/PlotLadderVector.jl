function short_step_correlations_vectors(
    df::DataFrame,
    var::Symbol;
    color_by::Symbol = :conn,
    remove_unstable::Bool = false
)
    step_keys   = ["S1", "S2"]#, "S3"]
    step_names = ["Full Model", "Global A (Global ϵ)"]#, "Global AE"]

    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        @info "Filtered unstable: now $(nrow(df)) rows"
    end

    full_col   = Symbol(string(var)*"_full")
    panel_cols = Symbol.(string(var) .* "_" .* step_keys)

    full_vals  = df[!, full_col]    # Vector of Vector{Float64}
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    ns   = length(panel_cols)
    cols = 3
    rows = ceil(Int, ns/cols)

    fig = Figure(; size=(1000, 600))
    Label(fig[0, 1:cols], uppercase(string(var, " correlations")), fontsize=18)

    for idx in 1:ns
        r = div(idx-1, cols) + 1
        c = mod(idx-1, cols) + 1
        ax = Axis(fig[r, c];
            title  = step_names[idx],
            xlabel = string(full_col),
            ylabel = string(panel_cols[idx]),
        )

        xs = Float64[]
        ys = Float64[]
        cs = Float64[]

        for (i, fvvec) in enumerate(full_vals)
            vec = df[i, panel_cols[idx]]
            @assert length(fvvec) == length(vec) "length mismatch in row $i"
            for j in eachindex(vec)
                push!(xs, fvvec[j])
                push!(ys,    vec[j])
                push!(cs, color_vals[i])
            end
        end

        scatter!(ax, xs, ys;
            colormap   = :inferno,
            color      = cs,
            colorrange = (cmin, cmax),
            markersize = 8,
            alpha      = 0.6
        )

        mn = min(extrema(xs)..., extrema(ys)...)
        mx = max(extrema(xs)..., extrema(ys)...)
        lines!(ax, [mn, mx], [mn, mx]; color=:black, linestyle=:dash)

        r_val = cor(xs, ys)
        text!(ax, "r=$(round(r_val, digits=3))",
            position=(mx, mn), align=(:right,:bottom), fontsize=10)
    end

    display(fig)
end

#TODO once drago finishes ShortAlloLadder_withMinExtinction_768.jls re-run this again
df = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder_withMinExtinction_48.jls")
df = deserialize("ThePaper/Ladder/ScenarioAllo/Outputs/ShortAlloLadder_withMinExtinction_768.jls")
df = A
short_step_correlations_vectors(
    df, :min_delta_K;
    color_by=:conn, remove_unstable=false
)

short_step_correlations_vectors(
    df, :min_delta_xi;
    color_by=:conn, remove_unstable=false
)



