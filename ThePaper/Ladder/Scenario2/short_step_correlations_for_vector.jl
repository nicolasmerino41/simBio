function short_step_correlations_vectors(
    df::DataFrame,
    varx::Symbol,
    vary::Symbol;
    compare_to_full::Bool = true,
    color_by::Symbol = :conn,
    remove_unstable::Bool = false,
    remove_zeros::Bool = false,
    equal_axes::Bool = true
)
    # step identifiers & titles
    step_keys  = ["S1"]#,"S2","S3","S4","S5","S6","S7","S8", "S9", "S10"]#, "S11"]
    step_names = [
        "Full Model"]    
    # "Full Model", "Global A (Global ϵ)", " Global AE",
    #     "Randomize m_cons ↻", "Randomize ξ̂ ↻", "Randomize K_res ↻",
    #     "Global A (Global ϵ) Mean B", "Global AE Mean B",
    #     "Rewire network", "Rewire network Randomly", "Rewire network with diff C"
    # ]

    if remove_unstable
        res_cols = Symbol.("resilience_" .* step_keys)
        df = filter(row -> all(row[c] < 0 for c in res_cols), df)
        @info "Filtered unstable: now $(nrow(df)) rows"
    end

    # build column symbols
    full_xcol  = Symbol(string(varx)*"_full")
    step_xcols = Symbol.(string(varx) .* "_" .* step_keys)
    step_ycols = Symbol.(string(vary) .* "_" .* step_keys)

    # color scale
    color_vals = df[!, color_by]
    cmin, cmax = extrema(color_vals)

    ns   = length(step_keys)
    cols = 4

    fig = Figure(size=(960, 465))
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

        # compute limits
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
            # still draw unity for reference
            lines!(ax, [xmin, xmax], [xmin, xmax], color=:black, linestyle=:dash)
        end

        # Pearson r
        r_val = cor(xs, ys)
        # r_val = corspearman(xs, ys)
        text!(ax, "r=$(round(r_val, digits=3))",
            position = (xmax, ymin),
            align    = (:right, :bottom),
            fontsize = 10
        )
    end
    Colorbar(fig[4, 4], colormap=:viridis, colorrange=(cmin, cmax))

    display(fig)
end


#TODO You must load "exploring_min_extinction_100000.jls" once done from drago and run the following plots
df = deserialize("ThePaper/Ladder/Outputs/exploring_min_extinction_5000_withsssp_rmed.jls")
df = deserialize("ThePaper/Ladder/Outputs/exploring_min_extinction_96_withsssp_rmed_tinyRmed.jls")

short_step_correlations_vectors(
    df, :ssp_rmed, :rt_pulse_vector;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=false
)
##############################################################################
##############################################################################
############################ TAU #############################################
##############################################################################
short_step_correlations_vectors(
    df, :ssp_rmed, :ssp_rmed;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=true
)

short_step_correlations_vectors(
    df, :rt_press_vector, :rt_press_vector;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=true
)

short_step_correlations_vectors(
    df, :rt_pulse_vector, :rt_pulse_vector;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=false
)

scen_value = Float64[]
for i in 1:nrow(df)
    if df[i, :scen] == :ER
        push!(scen_value, 1)
    elseif df[i, :scen] == :PL
        push!(scen_value, 2)
    else
        push!(scen_value, 3)
    end
end
df.scen_value = scen_value
short_step_correlations_vectors(
    df, :tau, :tau;
    compare_to_full=true,
    color_by=:scen_value,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=true
)
##############################################################################
##############################################################################
##############################################################################
##############################################################################
short_step_correlations_vectors(
    df, :rt_pulse_vector, :rt_pulse_vector;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=true,
    remove_zeros=true,
    equal_axes=false
)

short_step_correlations_vectors(
    df, :min_delta_K, :min_delta_K;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=false,
    remove_zeros=true,
    equal_axes=true
)

short_step_correlations_vectors(
    df, :ssp_rmed, :tau;
    compare_to_full=true,
    color_by=:conn,
    remove_unstable=true,
    remove_zeros=true,
    equal_axes=false
)
