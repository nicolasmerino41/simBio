######### Europe Species Df PART #########
# Read the first row to get column names
temp_df = CSV.File("JULIA Europe/europe_species_df.csv")
temp_df_df = temp_df |> DataFrame
for row in axes(temp_df_df, 1), col in 2:1150
    if temp_df_df[row, col] == "NA"
        temp_df_df[row, col] = "0"
    end
end

parsed_temp_df = parse.(Float64, temp_df_df[!, 2:end])
species_names = names(parsed_temp_df)
europe_species_df = hcat(temp_df_df[!, 1], parsed_temp_df)
europe_species_df = rename!(europe_species_df, :x1 => :Cell)

########## RASTER PART (maskk) ##########
maskk = Raster("JULIA Europe/mask.tif")
parent_mask = parent(maskk)
missingval(maskk)
maskk = replace_missing(maskk, NaN)
good_indices = findall(!isnan, maskk)
info_df = DataFrame(row = Int64[], col = Int64[], value = Float64[])
for i in CartesianIndices(maskk)
    if !isnan(maskk[i])
        push!(info_df, (i[1], i[2], maskk[i]))
    end
end
maskID = CSV.File("JULIA Europe/mask_ID.csv") |> DataFrame
maskID = maskID[!, [2,4]]

###### Abundance RASTER OF MyStructs1149 ######
raster_data_DA = [MyStructs1149(SVector{1149, Float64}(fill(0.0, 1149))) for _ in 1:671*589]
raster_DA = Raster(reshape(raster_data_DA, 671, 589), dims=dims(maskk))
i = good_indices[2] 
for i in good_indices
    row = i[1]
    col = i[2]
    value = info_df[findfirst((info_df.row .== row) .& (info_df.col .== col)), :value]
    pn = maskID[findfirst(maskID.Value .== value), :PageName]
    data_row = findfirst(europe_species_df.Cell .== pn)
    data = Vector(europe_species_df[data_row, 2:end])
    data_svector = SVector{1149, Float64}(data)
    raster_DA[row, col] = MyStructs1149(data_svector)
end

function Makie.convert_arguments(t::Type{<:Makie.Image}, A::AbstractRaster{<:MyStructs1149, 2})
    # Count presence based on the threshold
    richness = map(mystruct -> count(i -> mystruct.a[i] > 0.5, 1:length(mystruct.a)), A)
    return Makie.convert_arguments(t, richness)
end
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractRaster{<:MyStructs1149, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end

fig = Figure(resolution = (800, 600))
ax = Axis(fig[1,1])
hidexdecorations!(ax; grid=false)
hideydecorations!(ax; grid=false)
image!(ax, raster_DA; interpolate=false, colormap=custom_palette)
image!(ax, europe_sum; interpolate=false, colormap=custom_palette)
display(fig)

############### Europe_sum ###############
europe_sum_DA = [false for _ in 1:671*589]
europe_sum = Raster(reshape(europe_sum_DA, 671, 589), dims=dims(maskk))
for i in good_indices
    europe_sum[i] = true
end
########################
rule = Cell{:state, :state}() do data, state, I
    # if any(isinf, state.a) || any(isnan, state.a)
    #     @error "state has NA values"
    #     println(I)
    # end
    return MyStructs1149(SVector{1149, Float64}(state.a .+ 1))
end
pepe = (;
    state = raster_DA,
    europe_sum = europe_sum
)

array_output = ResultOutput(
    pepe, tspan = 1:100;
    mask = europe_sum
)

makie_output = MakieOutput(pepe, tspan = 1:10;
    fps = 10, ruleset = Ruleset(rule),
    mask = europe_sum) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:biomass, :simulated_richness, :europe_sum]
    titles = ["Biomass", "Simulated Richness", "europe_sum"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[i, j]; title=titles[(i-1)*2 + j]) for i in 1:1, j in 1:2]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state]; interpolate=false, colormap=custom_palette)
        elseif key == :simulated_richness
            Makie.image!(ax, frame[:state]; colormap=custom_palette)
        elseif key == :europe_sum
            Makie.heatmap!(ax, frame[:europe_sum]; interpolate=false, colormap=custom_palette)
        end
        hidexdecorations!(ax; grid=false)
        hideydecorations!(ax; grid=false)
        ax.title = title  # Set the title for each axis
        ax.titlegap[] = 5  # Adjust the title gap to make it smaller
        ax.titlesize[] = 12 # Set the title font size
        ax.titlecolor[] = RGBA(0, 0, 0, 1)  # Set the title color to black
        ax.yreversed[] = false
    end
end