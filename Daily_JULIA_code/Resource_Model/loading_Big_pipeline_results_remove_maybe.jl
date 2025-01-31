Big_P_results = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/29-1/Big_pipeline_results_drago_cheato.csv", DataFrame)
Big_P_results = Big_P_results[.!ismissing.(Big_P_results.cell_id), :]
unique_cells = string.(unique(Big_P_results.cell_id))

begin
    # Initialize an empty DataFrame with the same structure as Big_P_results
    Big_P_results_maximised = similar(Big_P_results, 0)

    for cell in unique(Big_P_results.cell_id)
        Big_P_results_cell = Big_P_results[Big_P_results.cell_id .== cell, :]
        ind = argmax(Big_P_results_cell.survival_rate)  # Get index of max survival_rate
        push!(Big_P_results_maximised, Big_P_results_cell[ind, :])  # Append the selected row
    end

    sr = DataFrame(survival_rate = Big_P_results_maximised.survival_rate)
    initial_part = Big_P_results_maximised[:, 1:3]
    final_part = hcat(Big_P_results_maximised[:, 4:14], Big_P_results_maximised[:, 16:end])
    Big_P_results_maximised = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = Big_P_results.survival_rate)
    initial_part = Big_P_results[:, 1:3]
    final_part = hcat(Big_P_results[:, 4:14], Big_P_results[:, 16:end])
    Big_P_results = hcat(initial_part, sr, final_part)
end
A_fixed = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/29-1/Big_pipeline_results_drago_fixed.csv", DataFrame)

begin    
    # Extract the `sp_removed_name` column
    sp_removed_name = Big_P_results_maximised.sp_removed_name
    # Convert to a DataFrame (ensuring it's properly structured for aggregation)
    sp_removed_name_df = DataFrame(x = sp_removed_name)
    # Count occurrences of each unique species name
    sp_removed_name_count = DF.combine(groupby(sp_removed_name_df, :x), nrow => :count)
    # Substitute the missing name by "none"
    sp_removed_name_count[ismissing.(sp_removed_name_count.x), :x] .= "none"
    # Sort in descending order of count
    sort!(sp_removed_name_count, :count, rev=true)

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title="Species removal frequency", xlabel="Species", ylabel="Count")
    MK.barplot!(ax, 1:length(sp_removed_name_count.x), sp_removed_name_count.count)
    ax.xticks = (1:length(sp_removed_name_count.x), sp_removed_name_count.x)
    ax.xticklabelrotation = π/2.5
    ax.xticklabelsize = 8
    display(fig)
end

extract_info_of_a_species("Capreolus capreolus")