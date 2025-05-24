using EcologicalNetworksDynamics, Plots
fw = Foodweb([1 => 2, 2 => 3]) # 1 eats 2, and 2 eats 3.

m = default_model(fw)

B0 = [0.1, 0.1, 0.1] # The 3 species start with a biomass of 0.1.
t = 100 # The simulation will run for 100 time units.
out = simulate(m, B0, t)

Plots.plot(out)

fw1 = Foodweb(:niche; S = 5, C = 0.2)

m = default_model(fw1)
m.FoodWeb

B0 = [1.0, 1.0, 1.0, 1.0, 1.0] # The 5 species start with a biomass of 1.0.
t = 100 # The simulation will run for 1000 time units.
out1 = simulate(m, B0, t)

Plots.plot(out1)

fw1