using CairoMakie

# Replace this with your actual DataFrame variable
# df = DataFrame(...) 
dff = copy(df)  # Assuming df is your DataFrame
res_cols = Symbol.("resilience_" .* step_keys)
dff = filter(row -> all(row[c] < 0 for c in res_cols), dff)
res_cols = Symbol.("after_pulse_" .* step_keys)
dff = filter(row -> all(row[c] .> 0.9 for c in res_cols), dff)
println("subset size: ", nrow(dff))

for i in 10:20
    rt_resources = Float64[]
    rt_consumers = Float64[]
    x_resources = Float64[]
    x_consumers = Float64[]
    for row in eachrow(dff)[i:i]  # Assuming dff is a DataFrame
        # rt_vec = row[:ssp_rmed_full]
        rt_vect = row[:rt_pulse_vector_full]
        @info "rt_vec = $rt_vect"
        tau_full = row[:tau_full]
        @info "tau_full = $tau_full"
        k_xi   = row[:K_Xi_full]
        R_eq   = row[:R_eq]
        C_eq   = row[:C_eq]
        # m_val  = row[:m_val]
        # g_val  = row[:g_val]
        m_val  = row[:m_cons]
        g_val  = row[:r_res]

        # # Resources: positions 1:30
        # for i in 1:30
        #     Ki = k_xi[i]
        #     rBi = R_eq[i]
        #     xval = g_val[i] * rBi / Ki
        #     push!(x_resources, xval)
        #     push!(rt_resources, rt_vec[i])
        # end

        # Consumers: positions 31:50
        for i in 1:20
            xi = k_xi[30+i]
            Bi = C_eq[i]
            xval = m_val[i] * Bi / xi
            push!(x_consumers, xval)
            push!(rt_consumers, rt_vect[30+i])
        end
    end

    begin
        fig = Figure(; size = (800, 500))
        ax = Axis(fig[1, 1];
            xlabel = "Self-Regulation Loss (rB/K for resources, mB/xi for consumers)",
            ylabel = "Return Time",
            title = "Self-Regulation Loss vs Return Time",
            titlealign = :center,
            ylabelsize = 18,
            xlabelsize = 18,
        )
        scatter!(ax, x_resources, rt_resources, color = :blue, markersize = 7, label = "Resources (rB/K)", alpha=0.6)
        scatter!(ax, x_consumers, rt_consumers, color = :red, markersize = 7, label = "Consumers (mB/xi)", alpha=0.6)

        axislegend(ax; position = :rt, framevisible=false)
        display(fig)
    end
end

using CairoMakie, Statistics, DataFrames, Logging

"""
    plot_selfreg_vs(df::DataFrame, rowidx::Int;
                   yvar::Symbol,
                   groups::Vector{Symbol} = [:resources, :consumers],
                   R::Int, C::Int)

For a single row of `df` (indexed by `rowidx`), plots self‐regulation loss (SL) on the x‐axis
against the per‐species quantity in `yvar` (a Vector{<:Real}) on the y‐axis.  

- `yvar` should name a column of `df` containing a length-(R+C) vector (e.g. `:rt_pulse_full_vector` or `:tau_full`).  
- `groups` chooses which trophic levels to plot; any subset of `[:resources, :consumers]`.  
- `R` and `C` are the numbers of resources and consumers.  
"""
function plot_SLloss_vs(df::DataFrame, rowidx::Int;
                          yvar::Symbol,
                          groups::Vector{Symbol} = [:resources, :consumers],
                          R::Int, C::Int)

    row = df[rowidx, :]
    SL       = row[:SL_vect]               # length R+C
    yvec     = row[yvar]                   # length R+C
    R_eq     = row[:R_eq]
    C_eq     = row[:C_eq]
    K_Xi     = row[:K_Xi_full]             # length R+C
    r_res    = row[:r_res]                 # length R
    m_cons   = row[:m_cons]                # length C

    xres = Float64[]
    yres = Float64[]

    xcon = Float64[]
    ycon = Float64[]

    # resources are indices 1:R
    if :resources in groups
        for i in 1:R
            Ki  = K_Xi[i]
            Bi  = R_eq[i]
            ξi  = nothing # unused for resources
            xval = r_res[i] * Bi / Ki
            push!(xres, SL[i])
            push!(yres,   yvec[i])
        end
    end

    # consumers are indices R+1:R+C
    if :consumers in groups
        for k in 1:C
            i   = R + k
            ξi  = K_Xi[i]
            Bi  = C_eq[k]
            xval = m_cons[k] * Bi / ξi
            push!(xcon, SL[i])
            push!(ycon,   yvec[i])
        end
    end

    fig = Figure(; size=(800,500))
    ax  = Axis(fig[1,1],
               xlabel = "Self‐regulation loss",
               ylabel = string(yvar),
               title  = "$(yvar) vs SL for row $rowidx")

    if :resources in groups
        scatter!(ax, xres, yres;
                 color=:blue, label="Resources", markersize=6, alpha=0.7)
    end
    if :consumers in groups
        scatter!(ax, xcon, ycon;
                 color=:red,  label="Consumers", markersize=6, alpha=0.7)
    end

    axislegend(ax, position=:rt, framevisible=false)
    display(fig)
    return fig
end

# Example usage:
for i in 1:15
    plot_SLloss_vs(
        dff, 2; yvar=:ssp_rmed_full,
        groups=[:resources, :consumers],
        R=0, C=20
    )
end