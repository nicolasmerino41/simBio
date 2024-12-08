using LinearAlgebra
using DifferentialEquations
using CairoMakie
using DiffEqCallbacks

# begin
    ########################################################
    # Parameters and Structures
    ########################################################
    mutable struct Herbivore
        id::Int
        m::Float64      # mortality rate m_i
        H0::Float64     # characteristic density H_i^0
        g::Float64      # growth rate g_i (will be computed)
        H_init::Float64 # initial abundance
    end

    mutable struct Predator
        id::Int
        m::Float64       # mortality rate m_α
        P_init::Float64  # initial abundance
    end

    ########################################################
    # Example Parameterization
    ########################################################

    # Number of species
    S = 10      # number of herbivore species
    R = 10      # number of predator species

    # Mean parameters for herbivores (example values)
    M_mean = 0.1     # mean herbivore mortality
    H0_mean = 100.0   # mean herbivore characteristic density
    NPP = 1000.0      # total Net Primary Production

    # Predator parameters (example values)
    m_alpha_mean = 0.05  # predator mortality mean
    epsilon = 0.1        # assimilation efficiency

    exponent_abundance = 0.5  # exponent for SAD (0=even, >0=more skewed)

    # Interaction parameters
    # mu_ij for herbivores (competition matrix)
    mu = 0.1
    mu_matrix = fill(mu, S, S)
    for i in 1:S
        mu_matrix[i, i] = 0.0 # no self interaction in mu_ij, as defined in doc
    end

    # V matrix from the doc = (I + μ)^{-1} or a constructed competition matrix
    # Here we just create a symmetric positive-definite competition matrix V
    V = Matrix{Float64}(I, S, S)
    for i in 1:S
        for j in 1:S
            if i != j
                V[i, j] = 1/(1+mu) # Simple choice; you can refine as needed
            end
        end
    end

    # a_{i\alpha} predation rate coefficients
    a_matrix = rand(S, R) .* 0.001 # Random small predation coefficients

    # A_{αβ} predator-predator interaction matrix 
    # We'll choose a stable structure, for instance, a negative diagonal dominance
    A = -0.01 * Matrix(I, R, R)  # predators have weak negative self-regulation
    for α in 1:R
        for β in 1:R
            if α != β
                A[α, β] = 0.0 # no predator-predator facilitation
            end
        end
    end

    # Predator mortality rates and initial conditions
    m_alpha = [m_alpha_mean for _ in 1:R]
    P_init_values = fill(10.0, R)  # initial predator densities

    # Construct herbivore and predator lists
    herbivores_list = Herbivore[]
    for i in 1:S
        m_i = rand(Normal(M_mean, M_mean/10))
        H0_i = H0_mean
        H_init_i = H0_mean
        push!(herbivores_list, Herbivore(i, m_i, H0_i, 0.0, H_init_i))
    end

    predator_list = Predator[]
    for α in 1:R
        push!(predator_list, Predator(α, m_alpha[α], P_init_values[α]))
    end

    ########################################################
    # Computing p_i, g_i, and scaling x from theory
    #
    # For simplicity, we'll assume p_i and g_i come from the
    # given NPP and from a chosen SAD (H_hat). In the doc, 
    # p_i is derived from H_hat and μ. We'll just pick a simple
    # SAD and then compute p_i as shown.
    ########################################################
    # Generate a non-even SAD (H_hat_i)
    function generate_Hi_hat(S::Int, exponent::Float64)
        ranks = collect(1:S)
        abundances = ranks .^ (-exponent)
        abundances /= sum(abundances)
        return abundances
    end
    H_hat = generate_Hi_hat(S, exponent_abundance)

    # 3. Define h_i
    # According to the doc, h_i > 0 and sum_i h_i = 1.
    # Let's choose h_i = 1/S for simplicity:
    h_i = fill(1.0/S, S)

    # Compute pi_hat
    # pi_hat = (Hi_hat + sum_j mu_ij Hj_hat) / h_i
    pi_hat = similar(H_hat)
    for i in 1:S
        interaction_sum = 0.0
        for j in 1:S
            interaction_sum += mu_matrix[i, j]*H_hat[j]
        end
        pi_hat[i] = (H_hat[i] + interaction_sum) / h_i[i]
    end

    # Normalize to get p_i
    p_i = pi_hat ./ sum(pi_hat)

    # Compute inner product <p | V p>
    Vp = V * p
    inner_product = dot(p, Vp)

    # Compute x
    barH = H0_mean   # mean characteristic density (approx)
    barM = M_mean     # mean mortality rate (approx)
    x = sqrt(NPP/(barH*barM*inner_product))

    # Compute g_i and assign to herbivores
    for i in 1:S
        herbivores_list[i].g = x * herbivores_list[i].m * p[i]
    end

    ########################################################
    # Compute C_{ij} and G_i from the doc
    ########################################################

    A_inv = inv(A)

    C = zeros(S, S)
    G = zeros(S)

    for i in 1:S
        for j in 1:S
            val = 0.0
            for α in 1:R
                for β in 1:R
                    val += epsilon * a_matrix[i, α] * A_inv[α, β] * a_matrix[j, β]
                end
            end
            C[i, j] = val
        end
    end

    for i in 1:S
        val = 0.0
        for β in 1:R
            for α in 1:R
                val += a_matrix[i, α]*A_inv[α, β]*m_alpha[β]
            end
        end
        G[i] = val
    end

    ########################################################
    # Define μ_ij + (C_ij H_i^0)/m_i for herbivores
    ########################################################

    M_modified = copy(mu_matrix)
    for i in 1:S
        for j in 1:S
            M_modified[i, j] = mu_matrix[i, j] + (C[i, j]*herbivores_list[i].H0)/herbivores_list[i].m
        end
    end

    ########################################################
    # ODE System
    ########################################################

    function ecosystem_dynamics!(du, u, p, t)
        herbivores_list, predator_list, mu_matrix, M_modified, a_matrix, A, C, G, epsilon = p

        S = length(herbivores_list)
        R = length(predator_list)

        H = @view u[1:S]
        P = @view u[S+1:S+R]

        duH = zeros(S)
        duP = zeros(R)

        for i in 1:S
            if H[i] > 0.0
                sp = herbivores_list[i]
                m_i = sp.m
                g_i = sp.g
                H_i0 = sp.H0

                interaction_sum = sum(M_modified[i, j]*H[j] for j in 1:S)
                bracket = ((g_i + G[i])/m_i - 1) - (H[i] + interaction_sum)/H_i0
                duH[i] = H[i]*m_i*bracket
            else
                duH[i] = 0.0
            end
        end

        for α in 1:R
            if P[α] > 0.0
                pred = predator_list[α]
                m_alpha = pred.m

                predation_sum = sum(a_matrix[j, α]*H[j] for j in 1:S)
                predator_interactions = sum(A[α, β]*P[β] for β in 1:R)

                duP[α] = P[α]*(epsilon*predation_sum - m_alpha + predator_interactions)
            else
                duP[α] = 0.0
            end
        end

        du[1:S] = duH
        du[S+1:S+R] = duP
    end

    ########################################################
    # Initial Conditions and Problem Setup
    ########################################################

    H_init_values = [sp.H_init for sp in herbivores_list]
    P_init_values = [pred.P_init for pred in predator_list]

    u0 = vcat(H_init_values, P_init_values)
    tspan = (0.0, 1000.0)

    params = (herbivores_list, predator_list, mu_matrix, M_modified, a_matrix, A, C, G, epsilon)

    EXTINCTION_THRESHOLD = 1.0
    callbacks = []
    push!(callbacks, PositiveDomain())

    for i in 1:S
        condition(u, t, integrator) = u[i] - EXTINCTION_THRESHOLD
        affect!(integrator) = integrator.u[i] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    offset = S
    for α in 1:R
        idx = offset + α
        condition(u, t, integrator) = u[idx] - EXTINCTION_THRESHOLD
        affect!(integrator) = integrator.u[idx] = 0.0
        push!(callbacks, ContinuousCallback(condition, affect!))
    end

    cb = CallbackSet(callbacks...)

    prob = ODEProblem(ecosystem_dynamics!, u0, tspan, params)
    sol = solve(prob, Tsit5(); callback=cb, reltol=1e-6, abstol=1e-6)

    ########################################################
    # Results and Plot
    ########################################################

    times = sol.t
    herbivore_data = sol[1:S, :]
    predator_data = sol[S+1:S+R, :]

    num_survived_herb = count(herbivore_data[:, end] .> EXTINCTION_THRESHOLD)
    num_survived_pred = count(predator_data[:, end] .> EXTINCTION_THRESHOLD)

    println("$num_survived_herb/$S herbivore species survived.")
    println("$num_survived_pred/$R predator species survived.")

    fig = Figure(resolution=(800,600))
    ax = Axis(fig[1,1], xlabel="Time", ylabel="Density", title="Herbivores and Predators Dynamics")

    for i in 1:S
        lines!(ax, times, herbivore_data[i, :], label="H$i")
    end
    for α in 1:R
        lines!(ax, times, predator_data[α, :], label="P$α", linestyle=:dash)
    end

    axislegend(ax, position=:rt)
    display(fig)
# end

