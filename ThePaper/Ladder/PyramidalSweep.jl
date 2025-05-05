using DifferentialEquations, Random, Statistics, DataFrames, CairoMakie
import Base.Threads: @threads

# assume: make_A, calibrate_params, compute_jacobian, trophic_ode!,
#          build_callbacks, EXTINCTION_THRESHOLD, analytical_vs_sim
function pyramidal_sweep(;
    S::Int = 30,
    conn::Float64 = 0.1,
    IS::Float64   = 1.0,
    d_val::Float64= 1.0,
    m_val::Float64= 0.1,
    eps_val::Float64 = 0.1,
    C_ratio::Float64 = 0.2,
    scenario::Symbol = :ER,
    abundance_type::Symbol = :Normal,
    skew_vals::Vector{Float64} = [0.01,0.1,1.0,10.0],
    max_calib::Int = 10,
    abundance_mean::Float64 = 1.0,
    iterations::Int = 10
)
    R = Int(round(S*(1-C_ratio)))
    C = S - R
    cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    results = Vector{NamedTuple}()
    
    for skew in skew_vals, i in 1:iterations
        # predeclare so they exist after the while
        p = nothing
        R_eq_final = nothing
        C_eq_final = nothing

        # prepare placeholders
        m_cons = fill(m_val, C)
        d_res  = fill(d_val, R)

        # prepare placeholders
        xi_cons, r_res = fill(NaN, C), fill(NaN, R)
        tries = 0
        # try up to max_calib times
        while (any(isnan, xi_cons) || any(isnan, r_res)) && tries < max_calib
            tries += 1
            # (re)build A & ε
            A = make_A(zeros(S,S), R, conn, scenario) .* IS
            ε = clamp.(rand(Normal(eps_val, eps_val), S, S), 0.0, 1.0)

            # draw equilibrium with this skew
            if abundance_type == :Log
                R_eq = abs.(rand(LogNormal(
                    log(abundance_mean) - (abundance_mean^2)/2,
                    abundance_mean
                ), R))
                C_eq = abs.(rand(LogNormal(
                    log(abundance_mean*skew) - ((abundance_mean*skew)^2)/2,
                    abundance_mean*skew
                ), C))
            else
                R_eq = abs.(rand(Normal(abundance_mean, abundance_mean*0.1), R))
                C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew*0.1), C))
            end
            p_s = (R, C, m_cons, d_res, ε, A)
            # calibrate
            xi_cons, r_res = calibrate_params(
                R_eq, C_eq, p_s;
                xi_threshold=0.7,
                constraints=true
            )
            # println("hey")
            if !any(isnan, xi_cons)
                # capture the working values
                p           = (R, C, fill(m_val,C), xi_cons, r_res, fill(d_val,R), ε, A)
                R_eq_final  = R_eq
                C_eq_final  = C_eq
                break
            end
        end

        # if calibration never succeeded, skip
        if p === nothing
            @warn "calibration failed for skew=$skew"
            continue
        end

        # now do the perturbation‐vs‐analytic test
        fixed = vcat(R_eq_final, C_eq_final)
        sens_corr, extinct = analytical_vs_sim(p, fixed; cb=cb, δξ=1.0)
        if extinct
            @warn "extinction on perturbation at skew=$skew"
            continue
        end

        push!(results, (
            skew          = skew,
            pyramid_ratio = mean(R_eq_final)/mean(C_eq_final),
            sens_corr     = sens_corr
        ))
    end

    return DataFrame(results)
end

# run the sweep
A = pyramidal_sweep(;
    S=30, conn=0.1, IS=0.1, d_val=1.0, m_val=0.1, eps_val=0.1,
    C_ratio=0.1, scenario=:PL, abundance_type=:Normal,
    skew_vals=[0.001,0.01,0.05,0.1,0.5,1.0,10.0], max_calib=100,
    abundance_mean=100.0, iterations=10
)

# plot sens_corr vs pyramid_ratio
fig = Figure(resolution=(600,400))
ax  = Axis(fig[1,1],
           xlabel="Mean R / Mean C",
           ylabel="Sensitivity corr.",
           title="Analytic vs Sim Correlation across Pyramid Shapes")
scatter!(ax, df.pyramid_ratio, df.sens_corr; markersize=8)
lines!(ax, df.pyramid_ratio, df.sens_corr; linestyle=:dash)
display(fig)
