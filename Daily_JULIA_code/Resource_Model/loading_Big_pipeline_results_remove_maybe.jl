Big_P_results = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/29-1/Big_pipeline_results_drago_cheato.csv", DataFrame)
Big_P_results = Big_P_results[.!ismissing.(Big_P_results.cell_id), :]
Big_P_results = Big_P_results[Big_P_results.cell_id .!= "tia cetti", :]
unique_cells = string.(unique(Big_P_results.cell_id))

#### We needed to loaded this way cause there was column problems but it works fine
lines = readlines("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/30-1/Big_pipeline_results_with_even_pi_express.csv")
data = [join(split(line, ",")[1:24], ",") for line in lines if length(split(line, ",")) == 24]
# Now parse from an IOBuffer.
Big_P_even_pi = CSV.read(IOBuffer(join(data, "\n")), DataFrame, header=true)
Big_P_even_pi = Big_P_even_pi[.!ismissing.(Big_P_even_pi.cell_id), :][:, 1:24]
unique_cells_pi = string.(unique(Big_P_even_pi.cell_id))

#### We needed to loaded this way cause there was column problems but it works fine
lines = readlines("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/10-2/Big_pipeline_results_(special)_not_even_pi.csv")
data = [join(split(line, ",")[1:24], ",") for line in lines if length(split(line, ",")) == 24]
# Now parse from an IOBuffer.
SPECIAL_Big_P_not_even_pi = CSV.read(IOBuffer(join(data, "\n")), DataFrame, header=true)
SPECIAL_Big_P_not_even_pi = SPECIAL_Big_P_not_even_pi[.!ismissing.(SPECIAL_Big_P_not_even_pi.cell_id), :][:, 1:24]
unique_cells_pi = string.(unique(SPECIAL_Big_P_not_even_pi.cell_id))

begin
    # Initialize an empty DataFrame with the same structure as Big_P_results
    Big_P_results_maximised = similar(Big_P_results, 0)
    Big_P_even_pi_maximised = similar(Big_P_even_pi, 0)
    SPECIAL_Big_P_not_even_pi_maximised = similar(SPECIAL_Big_P_not_even_pi, 0)

    for cell in unique(Big_P_results.cell_id)
        Big_P_results_cell = Big_P_results[Big_P_results.cell_id .== cell, :]
        ind = argmax(Big_P_results_cell.survival_rate)  # Get index of max survival_rate
        push!(Big_P_results_maximised, Big_P_results_cell[ind, :])  # Append the selected row
    end

    for cell in unique(Big_P_even_pi.cell_id)
        Big_P_even_pi_cell = Big_P_even_pi[Big_P_even_pi.cell_id .== cell, :]
        ind = argmax(Big_P_even_pi_cell.survival_rate)  # Get index of max survival_rate
        push!(Big_P_even_pi_maximised, Big_P_even_pi_cell[ind, :])  # Append the selected row
    end

    for cell in unique(SPECIAL_Big_P_not_even_pi.cell_id)
        SPECIAL_Big_P_not_even_pi_cell = SPECIAL_Big_P_not_even_pi[SPECIAL_Big_P_not_even_pi.cell_id .== cell, :]
        ind = argmax(SPECIAL_Big_P_not_even_pi_cell.survival_rate)  # Get index of max survival_rate
        push!(SPECIAL_Big_P_not_even_pi_maximised, SPECIAL_Big_P_not_even_pi_cell[ind, :])  # Append the selected row
    end

    sr = DataFrame(survival_rate = Big_P_results_maximised.survival_rate)
    initial_part = Big_P_results_maximised[:, 1:3]
    final_part = hcat(Big_P_results_maximised[:, 4:14], Big_P_results_maximised[:, 16:end])
    Big_P_results_maximised = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = Big_P_even_pi_maximised.survival_rate)
    initial_part = Big_P_even_pi_maximised[:, 1:3]
    final_part = hcat(Big_P_even_pi_maximised[:, 4:14], Big_P_even_pi_maximised[:, 16:end])
    Big_P_even_pi_maximised = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = SPECIAL_Big_P_not_even_pi_maximised.survival_rate)
    initial_part = SPECIAL_Big_P_not_even_pi_maximised[:, 1:3]
    final_part = hcat(SPECIAL_Big_P_not_even_pi_maximised[:, 4:14], SPECIAL_Big_P_not_even_pi_maximised[:, 16:end])
    SPECIAL_Big_P_not_even_pi_maximised = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = Big_P_results.survival_rate)
    initial_part = Big_P_results[:, 1:3]
    final_part = hcat(Big_P_results[:, 4:14], Big_P_results[:, 16:end])
    Big_P_results = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = Big_P_even_pi.survival_rate)
    initial_part = Big_P_even_pi[:, 1:3]
    final_part = hcat(Big_P_even_pi[:, 4:14], Big_P_even_pi[:, 16:end])
    Big_P_even_pi = hcat(initial_part, sr, final_part)

    sr = DataFrame(survival_rate = SPECIAL_Big_P_not_even_pi.survival_rate)
    initial_part = SPECIAL_Big_P_not_even_pi[:, 1:3]
    final_part = hcat(SPECIAL_Big_P_not_even_pi[:, 4:14], SPECIAL_Big_P_not_even_pi[:, 16:end])
    SPECIAL_Big_P_not_even_pi = hcat(initial_part, sr, final_part)
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