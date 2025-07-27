function reorder_columns(df::DataFrame)
    # Metadata first
    meta_cols = [:conn, :scen, :IS, :delta, :marg, :ite]

    # Metric base names, in desired order
    metric_bases = [
        :resilience, :reactivity, 
        :after_press, :after_pulse,
        :rt_press, :rt_pulse,
        :S,
        :collectivity,
        :SL, :mean_SL, :sigma_over_min_d, :rmed
    ]

    # All metric columns: _full, then _S1–S5
    metric_cols = String[]
    for base in metric_bases
        push!(metric_cols, string(base, "_full"))
        for step in 1:5
            push!(metric_cols, string(base, "_S", step))
        end
    end

    # Final columns
    final_cols = ["p_final", "B_eq", "Beq_cv"]

    # Combined column order
    ordered_names = Symbol.(vcat(meta_cols, metric_cols, final_cols))

    # Subset and return reordered DataFrame
    return df[:, ordered_names]
end
