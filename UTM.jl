using Rasters, Shapefile, ArchGDAL, Plots
using CSV, DataFrames
meta_path = "C:\\Users\\MM-1\\OneDrive\\PhD\\Metaweb Modelling"

spain_vector = Shapefile.Table(joinpath(meta_path, "ATLAS_data", "UTMgrid_IberianPeninsula", "UTM_GRID.shp")) |> DataFrame
plot(spain_vector.geometry, color = :white, fill = true, alpha = 0.5)

Amph = CSV.read(joinpath(meta_path, "Atlas_data", "DB_Amphibians_IP.txt"), delim='\t', DataFrame) |> DataFrame

spain_vector_with_amph = leftjoin(spain_vector, Amph, on = :UTMCODE)
spain_vector_with_amph.rowsum = sum.(eachrow(spain_vector_with_amph[:, 3:end])) 

plot(spain_vector_with_amph.geometry, color=:thermal, clim=(minimum(spain_vector_with_amph.rowsum), maximum(spain_vector_with_amph.rowsum)), fill=true, alpha=0.5)

using PlotShapefiles
PlotShapefiles.choropleth(spain_vector_with_amph.geometry, spain_vector_with_amph[:rowsum], colormap("blues"))