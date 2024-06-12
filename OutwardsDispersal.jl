outdisp = OutwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
indisp = InwardsDispersal(
    formulation=ExponentialKernel(λ=0.0125),
    distancemethod=AreaToArea(30)
)
init = fill(0.0, (10,10))
init[5,5] = 1000000.0

mask = fill(true, (10,10))
for i in 1:10 # Adjust the range as needed depending on mask dimensions
    mask[1, i] = false
end

out = ArrayOutput(init; tspan=1:1000, mask=mask)
i = sim!(out, Ruleset(indisp; boundary=Reflect()))
o = sim!(out, Ruleset(outdisp; boundary=Reflect()))

sum(i[100])
sum(o[100])