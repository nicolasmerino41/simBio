hf = Rasters.Raster("Rasters/wildareas-v3-2009-human-footprint.tif", lazy = true)

hf_reprojected = Rasters.resample(hf; to=utmraster)

hf_spain = map(x -> Float32(x), hf_reprojected)
hf_spain = map(x -> x == 128.0 ? 0.0 : x, hf_spain)

inverted_hf = 1 ./ deepcopy(hf_spain)
inverted_hf = map(x -> x == Inf ? 0.0 : x, inverted_hf)

map_plot(inverted_hf; palette = custom_palette, flip = false);

######## FOR DA'S ##########
inverted_hf_DA = DimArray(Matrix(inverted_hf), (Dim{:a}(1:125), Dim{:b}(1:76)))
# Ensure values below 0.1 (but greater than 0.0) in inverted_hf_DA are set to 0.1
inverted_hf_DA .= ifelse.((inverted_hf_DA .> 0.0) .& (inverted_hf_DA .< 0.1), 0.1, inverted_hf_DA)

function adjust_inverted_hf_DA(lambda_DA, inverted_hf_DA)
    inverted_hf_DA1 = deepcopy(inverted_hf_DA)
    rows, cols = size(lambda_DA.multiplicative) # Any lambda is ok, it's just for size matching
    adjusted_hf_DA = deepcopy(inverted_hf_DA1)  # Start with a copy of the original matrix

    # First pass: Calculate averages
    for row in 1:rows
        for col in 1:cols
            if !isnan(lambda_DA.multiplicative[row, col])
                neighbors = get_neighbors(inverted_hf_DA1, row, col)
                non_zero_neighbors = filter(x -> x != 0.0, neighbors)

                if !isempty(non_zero_neighbors)
                    adjusted_hf_DA[row, col] = mean(non_zero_neighbors)
                else
                    adjusted_hf_DA[row, col] = 0.1  # or some other default value if all neighbors are zero
                end
            end
        end
    end

    # Second pass: Set NA cells in lambda_DA to 0.0 in adjusted_hf_DA
    for row in 1:rows
        for col in 1:cols
            if isnan(lambda_DA.multiplicative[row, col])
                adjusted_hf_DA[row, col] = 0.0
            end
        end
    end

    return adjusted_hf_DA
end

# Adjust the inverted_hf_DA raster based on lambda_DA
adjusted_inverted_hf_DA = adjust_inverted_hf_DA(lambda_DA, inverted_hf_DA)

k_DA_hf_multiplicative = k_DA.DA_multiplicative .* adjusted_inverted_hf_DA
k_DA_hf_additive = k_DA.DA_additive .* adjusted_inverted_hf_DA
k_DA_hf_geometric = k_DA.DA_geometric .* adjusted_inverted_hf_DA
k_DA_hf_minimum = k_DA.DA_min .* adjusted_inverted_hf_DA
k_DA_hf_harmonic = k_DA.DA_harmonic .* adjusted_inverted_hf_DA

########## FOR RASTER ###########
inverted_hf .= ifelse.((inverted_hf .> 0.0) .& (inverted_hf .< 0.1), 0.1, inverted_hf)

function adjust_inverted_hf(lambda_raster, inverted_hf)
    inverted_hf1 = deepcopy(inverted_hf)
    rows, cols = size(lambda_raster.multiplicative) # Any lambda is ok, it's just for size matching
    adjusted_hf = deepcopy(inverted_hf1)  # Start with a copy of the original matrix

    # First pass: Calculate averages
    for row in 1:rows
        for col in 1:cols
            if !isnan(lambda_raster.multiplicative[row, col])
                neighbors = get_neighbors(inverted_hf1, row, col)
                non_zero_neighbors = filter(x -> x != 0.0, neighbors)

                if !isempty(non_zero_neighbors)
                    adjusted_hf[row, col] = mean(non_zero_neighbors)
                else
                    adjusted_hf[row, col] = 0.1  # or some other default value if all neighbors are zero
                end
            end
        end
    end

    # Second pass: Set NA cells in lambda_raster to 0.0 in adjusted_hf
    for row in 1:rows
        for col in 1:cols
            if isnan(lambda_raster.multiplicative[row, col])
                adjusted_hf[row, col] = 0.0
            end
        end
    end

    return adjusted_hf
end

# Adjust the inverted_hf raster based on lambda_raster
adjusted_inverted_hf = adjust_inverted_hf(lambda_raster, inverted_hf)

k_raster_hf_multiplicative = k_raster.raster_multiplicative .* adjusted_inverted_hf
k_raster_hf_additive = k_raster.raster_additive .* adjusted_inverted_hf
k_raster_hf_geometric = k_raster.raster_geometric .* adjusted_inverted_hf
k_raster_hf_minimum = k_raster.raster_min .* adjusted_inverted_hf
k_raster_hf_harmonic = k_raster.raster_harmonic .* adjusted_inverted_hf

