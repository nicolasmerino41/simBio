function plot_simulation(sol; size_plot = (1200,600))
    fig = Figure(; size=size_plot)
    ax = Axis(fig[1, 1])
    for i in 1:size(sol[:,:], 1)
        lines!(ax, sol.t, sol[i, :])
    end
    display(fig)
end

# ─────────────────────────────────────────────────────────────────────────────
# 4) Modified plot_simulation + plot_community
# ─────────────────────────────────────────────────────────────────────────────
function plot_simulation(sol; Rb, Ri, Rt, size_plot=(1200,600))
    S = length(sol.u[1])
    levels = vcat(fill(:basal, Rb), fill(:inter, Ri), fill(:top, Rt))
    colors = Dict(:basal=>"forestgreen", :inter=>"royalblue", :top=>"tomato")

    fig = Figure(; size=size_plot)
    ax  = Axis(fig[1,1]; title="Dynamics (Rb=$Rb, Ri=$Ri, Rt=$Rt)",
                       xlabel="time", ylabel="Biomass")

    for i in 1:S
        lvl = levels[i]
        ts  = sol.t
        Bs  = getindex.(sol.u, i)
        lines!(ax, ts, Bs, color=colors[lvl], label=string(lvl))
    end

    # axislegend(ax; position=:rt)
    display(fig)
    return fig
end

function plot_simulation2(sol, sol2; Rb, Ri, Rt, size_plot=(1200,600))
    S = length(sol.u[1])
    levels = vcat(fill(:basal, Rb), fill(:inter, Ri), fill(:top, Rt))
    colors = Dict(:basal=>"forestgreen", :inter=>"royalblue", :top=>"tomato")

    fig = Figure(; size=size_plot)
    ax  = Axis(fig[1,1]; title="Dynamics (Rb=$Rb, Ri=$Ri, Rt=$Rt)",
                       xlabel="time", ylabel="Biomass")

    for i in 1:S
        lvl = levels[i]
        ts  = sol.t
        Bs  = getindex.(sol.u, i)
        lines!(ax, ts, Bs, color=colors[lvl], label=string(lvl))
    end

    ax2 = Axis(fig[1,2]; title="Dynamics (Rb=$Rb, Ri=$Ri, Rt=$Rt)",
                       xlabel="time", ylabel="Biomass")
    for i in 1:S
        lvl = levels[i]
        ts  = sol2.t
        Bs  = getindex.(sol2.u, i)
        lines!(ax2, ts, Bs, color=colors[lvl], label=string(lvl))
    end

    # axislegend(ax; position=:rt)
    display(fig)
    return fig
end

function plot_community2(df::DataFrame, combo_id::Int; T=200.0)
    row = df[df.combo_id .== combo_id, :]
    @assert nrow(row)==1 "combo_id $combo_id not found"
    row = first(eachrow(row))

    # re-run the ODE from saved params
    p    = (r=row.r, K=row.K, A=row.A_mat)
    cb = build_callbacks(row.Rb+row.Ri+row.Rt, EXTINCTION_THRESHOLD)
    prob = ODEProblem(glv!, row.Bstar, (0.0,T), p)
    sol  = solve(prob; callback = cb, abstol=1e-8, reltol=1e-8)
    
    p    = (r=row.r, K=row.K, A=row.A_mat)
    prob = ODEProblem(glv!, row.Bstar .* 0.9, (0.0,T), p)
    sol2 = solve(prob; callback = cb, abstol=1e-8, reltol=1e-8)
    
    return plot_simulation2(sol, sol2; Rb=row.Rb, Ri=row.Ri, Rt=row.Rt)
end

function plot_community(df::DataFrame, combo_id::Int; T=200.0)
    row = df[df.combo_id .== combo_id, :]
    @assert nrow(row)==1 "combo_id $combo_id not found"
    row = first(eachrow(row))

    # re-run the ODE from saved params
    p    = (r=row.r, K=row.K, A=row.A_mat)
    cb = build_callbacks(row.Rb+row.Ri+row.Rt, EXTINCTION_THRESHOLD)
    prob = ODEProblem(glv!, row.Bstar, (0.0,T), p)
    sol  = solve(prob; callback = cb, abstol=1e-8, reltol=1e-8)
    
    return plot_simulation(sol; Rb=row.Rb, Ri=row.Ri, Rt=row.Rt)
end

# ─────────────────────────────────────────────────────────────────────────────
# Example: plot the first 5 feasible communities
# ─────────────────────────────────────────────────────────────────────────────
for cid in df.combo_id[1:10]
    plot_community(df, cid)
end

for cid in df.combo_id[1:10]
    plot_community2(df, cid)
end


df = deserialize("ThePaper/Ladder/GLV/glvOutputs/communities.jls")[1:250, :]
