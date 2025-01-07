# Callbacks
callbacks = []
push!(callbacks, PositiveDomain())

for x in 1:S
    condition(u, t, integrator) = u[x] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[x] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

offset = S
for α in 1:R
    ind = offset + α
    condition(u, t, integrator) = u[ind] - EXTINCTION_THRESHOLD
    affect!(integrator) = (integrator.u[ind] = 0.0)
    push!(callbacks, ContinuousCallback(condition, affect!))
end

cb = CallbackSet(callbacks...)