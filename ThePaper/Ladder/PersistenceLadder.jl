using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
import Base.Threads: @threads
function persistence_sweep(;
    S::Int=30,
    C_vals=[3,6,9,12],
    conn_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    IS_vals=[0.1,1.0],
    skew_vals=[0.01,0.1,0.5],
    abundance_types=[:Normal,:Log],
    scenarios=[:ER,:PL,:MOD],
    # leave pex and modg fixed at 1.0
    eps_scales=[0.1],
    d_vals=[0.1, 1.0],
    m_vals=[0.1, 0.3],
    abundance_mean=10.0,
    number_of_combinations=1000,
    tspan=(0.,500.0),
    t_perturb=250.0,
    delta_vals=[1.0, 3.0],
    max_calib=10,
)
    locki = ReentrantLock()
    results = NamedTuple[]

    # full cartesian product
    combos = collect(Iterators.product(
        C_vals, conn_vals, IS_vals, skew_vals, abundance_types, scenarios, delta_vals, d_vals, m_vals
    ))
    global cb = build_callbacks(S, EXTINCTION_THRESHOLD)

    @threads for (C, conn, IS, skew, abund_type, scenario, delta, d, m) in combos[sample(1:length(combos), min(length(combos), number_of_combinations), replace=false)]
        R = S - C

        # 1) try calibrating
        xi_cons, r_res = fill(NaN,C), fill(NaN,R)
        p_final = nothing
        R_eq, C_eq = nothing, nothing
        tries = 0

        while any(isnan, xi_cons) && tries < max_calib
            tries += 1

            # build A and ε
            A = make_A(zeros(S,S), R, conn, scenario;
                        pareto_exponent=1.0, mod_gamma=1.0) .* IS
            ε_scale = rand(eps_scales)
            ε = clamp.(rand(Normal(ε_scale,ε_scale), S, S), 0, 1)

            # draw equilibrium
            if abund_type == :Log
            R_eq = abs.(rand(LogNormal(log(abundance_mean)-abundance_mean^2/2,
                                        abundance_mean), R))
            C_eq = abs.(rand(LogNormal(log(abundance_mean*skew)-
                                        (abundance_mean*skew)^2/2,
                                        abundance_mean*skew), C))
            else
            R_eq = abs.(rand(Normal(abundance_mean, abundance_mean*skew), R))
            C_eq = abs.(rand(Normal(abundance_mean*skew, abundance_mean*skew), C))
            end

            # calibrate
            xi_cons, r_res = calibrate_params(
            R_eq, C_eq,
            (R,C, fill(m,C), fill(d,R), ε, A);
            xi_threshold=0.7, constraints=true
            )

            if !any(isnan, xi_cons)
            p_final = (R, C, fill(m,C), xi_cons, r_res, fill(d,R), ε, A)
            break
            end
        end

        isnothing(p_final) && continue

        # 2) local stability + survival
        fixed = vcat(R_eq, C_eq)
        D,M = compute_jacobian(fixed,p_final)
        maximum(real, eigvals(D*M)) >= 0 && continue
        ok, B_eq = survives!(fixed,p_final; tspan=tspan, cb=cb)
        # !ok && continue

        # 3) full‐model persistence
        _,_,_, pers_full, _, _ = simulate_press_perturbation(
            B_eq, p_final, tspan, t_perturb, delta;
            solver=Tsit5(), cb=cb)

        # 4) ladder‐step persistence
        pers_steps = Dict(i => NaN for i in 1:16)
        for step in 1:16
            A_s, ε_s = transform_for_ladder_step(step, p_final[8], p_final[7])
            p_s = (R,C, p_final[3], p_final[4], p_final[5], p_final[6], ε_s, A_s)
            _,_,_, pstep, _, _ = simulate_press_perturbation(
            B_eq, p_s, tspan, t_perturb, delta;
            solver=Tsit5(), cb=cb, full_or_simple=false)
            pers_steps[step] = pstep
        end

        record = (
            C=C, conn=conn, IS=IS, skew=skew,
            abundance=abund_type, scenario=scenario,
            delta=delta, pers_full=pers_full,
            d=d, m=m,
            [ Symbol("pers_step_$i") => pers_steps[i] for i in 1:16 ]...,
            p_final=p_final,R_eq=R_eq,C_eq=C_eq
        )

        lock(locki) do
            push!(results, record)
        end
    end

    DataFrame(results)
end

# run the sweep
dfp = persistence_sweep(;
    S=30,
    C_vals=[3,6,9,12],
    conn_vals=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    IS_vals=[0.1,1.0, 2.0],
    skew_vals=[0.01,0.1,0.5, 1.0],
    abundance_types=[:Normal,:Log],
    scenarios=[:ER,:PL,:MOD],
    # leave pex and modg fixed at 1.0
    eps_scales=[0.1],
    d_vals=[0.1, 1.0],
    m_vals=[0.1, 0.3],
    abundance_mean=10.0,
    number_of_combinations=10000,
    tspan=(0.,500.0),
    t_perturb=250.0,
    delta_vals=[1.0, 3.0, 5.0],
    max_calib=10
)
# braoder_dfp = CSV.write("ThePaper/Ladder/Outputs/broader_persistence_sweep.csv", dfp)
limited_dfp = CSV.File("ThePaper/Ladder/Outputs/persistence_sweep.csv") |> DataFrame
broader_dfp = CSV.File("ThePaper/Ladder/Outputs/broader_persistence_sweep.csv") |> DataFrame
# plotting helper
function get_grid_position(step)
    # 4×4 grid: place step 1 in (1,1), step 2 in (1,2), … step 4 in (1,4),
    # step 5 in (2,1), etc.
    row = 1 + div(step-1, 4)
    col = 1 + mod(step-1, 4)
    return row, col
end

begin
    fig = Figure(; size=(1100,600))

    for step in 1:16
    r,c = get_grid_position(step)
    ax = Axis(fig[r,c];
        title="Step $step vs Full",
        xlabel="pers_full", ylabel="pers_step_$step")

    xs = dfp.pers_full
    ys = dfp[!, Symbol("pers_step_$step")]

    scatter!(ax, xs, ys; markersize=5, alpha=0.6)

    # 1:1 line
    lines!(ax, 0:1, 0:1; color=:black, linestyle=:dash)

    # annotate r in bottom right
    bad_idx = isnan.(xs) .| isnan.(ys)
    correlation = cor(xs[.!bad_idx], ys[.!bad_idx])
    text!(ax, "r=$(round(correlation, digits=2))",;
            position = (0.95,0.05), 
            align = (:right,:bottom), 
            fontsize=8)
    end

    display(fig)
end

function plot_persistence_grids(
    df::DataFrame;
    color_by::Symbol,
    facet_by::Symbol,
    steps::UnitRange=1:16,
    ncols::Int=4
)   
    df = filter(row -> row.pers_full > 0, df)
    facets = unique(df[!, facet_by])
    cats   = unique(df[!, color_by])
    palette = distinguishable_colors(length(cats))
    color_map = Dict(cats[i] => palette[i] for i in eachindex(cats))
    
    for f in facets
        sub = df[df[!, facet_by] .== f, :]
        nsteps = length(steps)
        nrows  = ceil(Int, nsteps/ncols)
        fig = Figure(; size=(900, 600))
        Label(fig[0, 2], "$facet_by = $f", fontsize = 24, tellwidth = false)

        for (idx, step) in enumerate(steps)
            row = 1 + div(idx-1, ncols)
            col = 1 + mod(idx-1, ncols)
            ax = Axis(fig[row, col];
                      title  = "step $step",
                      xlabel = "pers_full",
                      ylabel = "pers_step_$step")

            xs = sub.pers_full
            ys = sub[!, Symbol("pers_step_$step")]

            # now build a color array by indexing into our Dict
            cols = [ color_map[val] for val in sub[!, color_by] ]
            limits!(ax, (0.5, 1.05), (0.5, 1.05))
            scatter!(ax, xs, ys; color=cols, markersize=6, alpha=0.7)

            lines!(ax, 0:1, 0:1; color=:black, linestyle=:dash)

            bad_idx = isnan.(xs) .| isnan.(ys)
            correlation = cor(xs[.!bad_idx], ys[.!bad_idx])
            text!(ax, "r=$(round(correlation, digits=2))",;
                    position = (0.9,0.6), 
                    align = (:right,:bottom), 
                    fontsize=8)
        end
        Colorbar(fig[1, 5], limits = (minimum(cats), maximum(cats)), colormap = distinguishable_colors(length(unique(cats))))
        # # build a legend once, using small proxy scatters
        # legend_plots = [ scatter(fig, [], []; color=palette[i], label=string(cats[i])) 
        #                  for i in eachindex(cats) ]
        # Legend(fig[nrows, ncols], legend_plots, string.(cats);
        #        title = string(color_by))

        display(fig)
    end
end

plot_persistence_grids(
  dfp;
  color_by = :IS,
  facet_by = :C,
  steps    = 1:16,
  ncols    = 4
)