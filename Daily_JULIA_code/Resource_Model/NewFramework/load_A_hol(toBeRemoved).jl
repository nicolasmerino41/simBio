A_hol = deserialize("Daily_JULIA_code/Resource_Model/Best_params_&_other_outputs/21-02/ga_results_NF_hollingII.jls")

h_run, sol = herbivore_run(
    1, 
    A_hol[2, :].mu, A_hol[2, :].mu_predation, A_hol[2, :].epsilon_val, true, A_hol[2, :].m_alpha;
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

num_zeros = count(x -> x == 0.0, sol[:, end])
println("There are $num_zeros zeros in the final state.")

