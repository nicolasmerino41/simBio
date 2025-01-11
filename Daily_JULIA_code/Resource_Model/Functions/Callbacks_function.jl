function build_callbacks(S::Int, R::Int, EXTINCTION_THRESHOLD::Float64, T_ext::Float64, i_star::Int)
    # A) Continuous threshold-based extinctions
    callbacks = []

    # 1) Always ensure positivity
    push!(callbacks, PositiveDomain())

    # 2) Herbivores: set to zero if below EXTINCTION_THRESHOLD
    for x in 1:S
        function threshold_condition(u, t, integrator)
            # event if u[x] < EXTINCTION_THRESHOLD
            return u[x] - EXTINCTION_THRESHOLD
        end
        function threshold_affect!(integrator)
            integrator.u[x] = 0.0
        end
        push!(callbacks, ContinuousCallback(threshold_condition, threshold_affect!))
    end

    # 3) Predators in [S+1..S+R]: set to zero if below EXTINCTION_THRESHOLD
    offset = S
    for α in 1:R
        ind = offset + α
        function threshold_condition(u, t, integrator)
            return u[ind] - EXTINCTION_THRESHOLD
        end
        function threshold_affect!(integrator)
            integrator.u[ind] = 0.0
        end
        push!(callbacks, ContinuousCallback(threshold_condition, threshold_affect!))
    end

    # Build callback set => no forced extinction
    cb_no_trigger = CallbackSet(callbacks...)

    # ============= B) Discrete Forced Extinction of i_star at time T_ext =============
    # 1) Condition function returns boolean once time surpasses T_ext
    function cause_extinction_condition(u, t, integrator)
        # Trigger once t crosses T_ext (use a small tolerance or a direct comparison)
        return (t > T_ext)
    end

    # 2) The effect sets species i_star to zero
    function cause_extinction_affect!(integrator)
        println(">>> Forcibly caused extinction of species index $i_star at time ", integrator.t)
        integrator.u[i_star] = 0.0
    end

    # Build the discrete callback
    forced_extinction_callback = DiscreteCallback(
        cause_extinction_condition,
        cause_extinction_affect!
    )

    # Combine threshold-based callbacks with forced-extinction
    cb_trigger = CallbackSet(callbacks..., forced_extinction_callback)

    return cb_no_trigger, cb_trigger
end