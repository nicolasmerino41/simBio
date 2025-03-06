A_hol = deserialize("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF_hollingII_with_callback.jls")
A_csv = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF_hollingII_with_callback.csv", DataFrame)
A_csv.Column12

# Load the CSV file
A_csv = CSV.read("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF_hollingII_with_callback.csv", DataFrame)

# Function to clean and convert each row
function parse_abundance(row)
    cleaned = replace(row, "None" => "NaN")  # Replace "None" with NaN
    return parse.(Float64, split(cleaned, ";"))  # Convert to array of Float64
end

# Apply the function to the column
A_csv.Column12 = [parse_abundance(row) for row in A_csv.Column12]
A_csv.Column13 = [parse_abundance(row) for row in A_csv.Column13]

# Verify transformation
println(A_csv.Column12[47])  # Should now be an array of numbers

begin
    
    mu = 0.002199
    mu_pred = 0.00186861
    epsilon = 01274528885054175
    m_alpha = 0.020223714919694744
    h_run, sol = herbivore_run(
        1, 
        mu, mu_pred, epsilon, true, m_alpha;
        include_predators=true, time_end=2000.0, plot=true,
        do_you_want_params = false,
        do_you_want_sol = true,
        NPP=Float64(npp_DA_relative_to_1000[idx[1][1], idx[1][2]]), artificial_pi=false,
        alpha=0.25,
        ignore_inf_error = true,
        hollingII = true, h = 0.1,
        H_init = A_hol[2, :].H_eq[1],
        P_init = A_hol[2, :].P_eq[1]
    )
end
num_zeros = count(x -> x == 0.0, sol[:, end])
println("There are $num_zeros zeros in the final state.")

