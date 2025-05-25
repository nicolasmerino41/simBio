using EcologicalNetworksDynamics, DifferentialEquations, ForwardDiff

# 1) Build a niche food‐web and default model
fw = Foodweb(:niche; S = 50, C = 0.3)
m  = default_model(fw)

# 2) Extract parameters from m
S       = m.species.number             # total S
r       = m.r                          # intrinsic growth rates r_i
K       = m.K                          # carrying capacities K_i
x       = m.x                          # metabolic rates x_i
d       = m.d                          # background mortality d_i
y       = m.y                          # max ingestion y_i
e       = m.efficiency                 # ε[i,j]
A       = m.trophic.matrix             # A[i,j]
ω       = m.w                          # ω[i,j]
B0_ref  = m.half_saturation_density    # B0_i
i0      = m.intraspecific_interference # i0_i
h       = m.hill_exponent              # h

# Pack into one tuple
p = (S, r, K, x, d, y, e, ω, B0_ref, i0, h)

# 3) Define Yodzis–Innes allometric ODE
function allo_ode!(du, u, p, t)
    S, r, K, x, d, y, e, ω, B0_ref, i0, h = p

    # zero derivatives
    @inbounds for i in 1:S
        du[i] = 0.0
    end

    # precompute B^h
    Bh = u .^ h

    @inbounds for i in 1:S
        # (1a) growth + metabolism + mortality
        Gi = r[i] > 0 ? (1 - u[i]/K[i]) : 0.0
        du[i] = r[i]*u[i]*Gi - x[i]*u[i] - d[i]*u[i]

        gain = 0.0
        loss = 0.0

        # trophic terms
        @inbounds for j in 1:S
            # i eats j?
            if ω[i,j] > 0
                denom = B0_ref[i]^h + i0[i]*Bh[i]
                for k in 1:S
                    if ω[i,k] > 0
                        denom += ω[i,k]*Bh[k]
                    end
                end
                Fij = ω[i,j]*Bh[j] / denom
                gain += x[i]*y[i]*u[i]*Fij
            end

            # j eats i?
            if ω[j,i] > 0
                denom = B0_ref[j]^h + i0[j]*Bh[j]
                for k in 1:S
                    if ω[j,k] > 0
                        denom += ω[j,k]*Bh[k]
                    end
                end
                Fji = ω[j,i]*Bh[i] / denom
                loss += x[j]*y[j]*u[j]*Fji / e[j,i]
            end
        end

        du[i] += gain - loss
    end
end

# 4) Solve our own ODE
cb = build_callbacks(50, EXTINCTION_THRESHOLD)
u0   = rand(S)  # initial biomasses
prob = ODEProblem(allo_ode!, u0, (0.0, 500.0), p)
sol1 = solve(prob, Tsit5(); callback = cb, abstol=1e-6, reltol=1e-6)

plot_simulation(sol1)

# 5) Solve via EcologicalNetworksDynamics.simulate
sol2 = simulate(m, u0, 500.0)
plot_simulation(sol2)

# 6) Compare final biomasses
B1 = sol1.u[end]
B2 = sol2.u[end]

@info "Maximum absolute difference: ", maximum(abs.(B1 .- B2))

survival = sum(B1 .> 0.0) / S
@info "Survival rate: ", survival

survival2 = sum(B2 .> 0.0) / S
@info "Survival rate (EcologicalNetworksDynamics): ", survival2