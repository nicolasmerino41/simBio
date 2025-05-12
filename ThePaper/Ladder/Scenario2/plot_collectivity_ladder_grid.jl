using CairoMakie, Random

function plot_collectivity_ladder_grid(
    ncases::Int=16;
    S::Int=50, C::Int=20,
    conn_vals=[0.05,0.1,0.2,0.4],
    IS_vals=[0.1,1.0,2.0, 3.0],
    scenarios=[:ER,:PL,:MOD],
    eps_scales=[0.1,0.5,1.0]
)
    R = S - C

    # build all possible combos and sample ncases of them
    combos = collect(Iterators.product(conn_vals, IS_vals, scenarios, eps_scales))
    picked = sample(combos, ncases; replace=false)

    fig = Figure(; size = (1000,650))
    for (i, (conn, IS, scen, eps_scale)) in enumerate(picked)
        # 1) full A & ε
        A_full = make_A(zeros(S,S), R, conn, scen; IS = IS)
        ε_full = clamp.(randn(S,S) .* eps_scale, 0, 1)

        # 2) full φ
        φ_full = compute_collectivity(A_full, ε_full)

        # 3) φ through the ladder steps 1–16 (we’ll index 0→16 so 0=full)
        φs = Float64[]
        for step in 0:16
            if step == 0
                push!(φs, φ_full)
            else
                Aₛ, εₛ = transform_for_ladder_step(step, A_full, ε_full)
                push!(φs, compute_collectivity(Aₛ, εₛ))
            end
        end

        # 4) subplot in a 4×4 grid
        row = div(i-1, 4) + 1
        col = mod(i-1, 4) + 1
        ax = Axis(fig[row, col];
            title = "case $i: conn=$(conn), IS=$(IS), $scen, ε-scale=$(eps_scale)",
            xlabel = "Ladder step",
            ylabel = "φ (collectivity)",
            xlabelsize = 8,
            ylabelsize = 8,
            titlesize = 10,
            xticks = 1:16,
            xticklabelsize = 8
        )
        limits!(ax, (1, 16), (0.0, 10.0))

        lines!(ax, 0:16, φs; color = :blue)
        scatter!(ax, 0:16, φs; color = :blue, markersize=6)
        # dashed line showing φ_full
        hlines!(ax, [φ_full]; linestyle = :dash, color = :black)
    end

    fig
end

# Call it!
plot_collectivity_ladder_grid()
