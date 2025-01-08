# ============= A) Continuous Threshold Extinctions =============
callbacks = []

# 1) Always ensure positivity
push!(callbacks, PositiveDomain())

# 2) For each herbivore species i in [1..S], set to zero if below EXTINCTION_THRESHOLD
for x in 1:S
    function threshold_condition(u, t, integrator)
        # Return a float expression whose sign crossing triggers the event
        # We want event if u[x] < EXTINCTION_THRESHOLD => (u[x] - EXTINCTION_THRESHOLD) crosses zero
        return u[x] - EXTINCTION_THRESHOLD
    end
    function threshold_affect!(integrator)
        integrator.u[x] = 0.0
    end
    push!(callbacks, ContinuousCallback(threshold_condition, threshold_affect!))
end

# 3) For each predator species α in [S+1..S+R]
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

# 1) Condition function returns a float: t - T_ext
function cause_extinction_condition(u, t, integrator)
    # The discrete callback triggers the event once (t - T_ext) crosses 0 from negative to positive
    return isapprox(t - T_ext, 0.0, atol=1.0) ? true : false
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

# Combine the threshold callbacks with the forced-extinction callback
# if you want both
cb_trigger = CallbackSet(callbacks..., forced_extinction_callback)