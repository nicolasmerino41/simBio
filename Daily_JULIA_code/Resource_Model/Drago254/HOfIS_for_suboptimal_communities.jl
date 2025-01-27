################################################################################
# Instead of storing the best configurations (max survivors) in a global DataFrame,
# we write each cell's best result (i.e. the line of best configuration) directly
# to a CSV file immediately after processing that cell.
################################################################################
include("prior.jl")
include("prior2.jl")
include("DA_birmmals_with_pi.jl")
include("generate_competition_matrix.jl")
include("species_dict.jl")

include("ecosystem_dynamics!.jl")
include("FI_functions.jl")
include("extract_H0_DA.jl")
include("attempt_setup_community.jl")
include("Callbacks_function.jl")
include("npp_DA_relative_to_1000.jl")
# Import SpinLock for thread-safety
using Base.Threads: SpinLock
using Combinatorics

# --------------------------------------------------
# GLOBALS & LOCK
# --------------------------------------------------
const file_lock = SpinLock()
const RESULTS_CSV = "best_scenarios.csv"
const MAX_ITERS           = 50000
const SURVIVAL_THRESHOLD  = 0.0

# Helper to write a row to CSV under a lock
function write_row_to_csv(rownt::NamedTuple; filepath::String=RESULTS_CSV)
    lock(file_lock)
    try
        df = DataFrame([rownt])  # single-row DataFrame
        if !isfile(filepath)
            # write with header
            CSV.write(filepath, df)
        else
            # append without header
            CSV.write(filepath, df; append=true, writeheader=false)
        end
    finally
        unlock(file_lock)
    end
end

# --------------------------------------------------
# RUN SIMULATION
# --------------------------------------------------
function run_simulation(u0, params; max_time=500.0, extinction_threshold=1e-6)
    S2, R2 = params[1], params[2]

    prob = ODEProblem(ecosystem_dynamics!, u0, (0.0, max_time), params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        solve(prob, Tsit5(); reltol=1e-7, abstol=1e-8)
    end

    if sol.t[end] < max_time || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        return (sr=0.0, survived_herb=0, survived_pred=0,
                H_bio=0.0, P_bio=0.0, total_bio=0.0,
                H_end=zeros(S2), P_end=zeros(R2))
    end

    H_end = sol[1:S2, end]
    P_end = sol[S2+1:S2+R2, end]

    # zero out
    H_end[H_end .< extinction_threshold] .= 0.0
    P_end[P_end .< extinction_threshold] .= 0.0

    survived_herb = count(H_end .> 0.0)
    survived_pred = count(P_end .> 0.0)
    total_surv = survived_herb + survived_pred
    total_sp = S2 + R2
    sr = total_sp>0 ? total_surv / total_sp : 0.0

    H_bio = sum(H_end)
    P_bio = sum(P_end)

    return (
        sr         = sr,
        survived_herb = survived_herb,
        survived_pred = survived_pred,
        H_bio      = H_bio,
        P_bio      = P_bio,
        total_bio  = H_bio + P_bio,
        H_end      = H_end,
        P_end      = P_end
    )
end

# --------------------------------------------------
# FIND BEST PARAM ALL SPECIES
# --------------------------------------------------
function find_best_param_all_species(cell, param_combos, sp_names, localNPP, localH0_vector)

    best_survival = 0.0
    best_solutions = NamedTuple[]  # store ties for best sr

    for combo in param_combos[1:min(end, MAX_ITERS)]
        mu_val, mu_pred, eps_val, sym_comp = combo

        setup = attempt_setup_community(
            idx[cell][1], idx[cell][2],
            mu_val, mu_pred, eps_val, sym_comp;
            localNPP=localNPP, localH0_vector=localH0_vector, species_names=sp_names
        )
        if setup === nothing
            continue
        end

        (S2, R2, H_i0, m_i, g_i, G, M_mod, a_matrix, A, eps_vec, m_alpha) =
            (setup.S, setup.R, setup.H_i0, setup.m_i, setup.g_i, setup.G,
             setup.M_modified, setup.a_matrix, setup.A, setup.epsilon_vector, setup.m_alpha)

        if S2+R2==0 || R2>length(H_i0)
            continue
        end

        H_init = H_i0
        P_init = H_init[1:R2] ./ 10.0
        u0 = vcat(H_init, P_init)

        params = (S2, R2, H_i0, m_i, g_i, G, M_mod, a_matrix, A, eps_vec, m_alpha)
        sim = run_simulation(u0, params)
        sr = sim.sr

        giHi = sum(setup.g_i .* sim.H_end)  # final herb biomass * final g_i
        ratio = giHi / localNPP
        ratio_ok = true
        # ratio_ok = (ratio > 0.5) && (ratio < 10.0)

        if sr > best_survival && ratio_ok
            best_survival = sr
            best_solutions = [
                (
                    cell_id = cell,
                    mu = mu_val,
                    mu_predation = mu_pred,
                    epsilon_val = eps_val,
                    sym_comp = sym_comp,
                    S2 = S2,
                    R2 = R2,
                    survival_rate = sr,
                    survived_herbivores = sim.survived_herb,
                    survived_predators  = sim.survived_pred,
                    H_biomass = sim.H_bio,
                    P_biomass = sim.P_bio,
                    total_biomass = sim.total_bio,
                    g_iH_i = giHi,
                    g_iH_i_over_NPP = ratio
                )
            ]
        elseif isapprox(sr, best_survival; atol=1e-10) && ratio_ok
            push!(best_solutions, (
                cell_id = cell,
                mu = mu_val,
                mu_predation = mu_pred,
                epsilon_val = eps_val,
                sym_comp = sym_comp,
                S2 = S2,
                R2 = R2,
                survival_rate = sr,
                survived_herbivores = sim.survived_herb,
                survived_predators  = sim.survived_pred,
                H_biomass = sim.H_bio,
                P_biomass = sim.P_bio,
                total_biomass = sim.total_bio,
                g_iH_i = giHi,
                g_iH_i_over_NPP = ratio
            ))
        end

        # If sr=1.0 + ratio_ok, break early
        if isapprox(sr, 1.0; atol=1e-10) && ratio_ok
            break
        end
    end
    return (best_survival, best_solutions)
end

# --------------------------------------------------
# SIMULATE REMOVAL
# --------------------------------------------------
function simulate_removal(best_sol, removal_list::Vector{String}, sp_nm_all::Vector{String}, localNPP, localH0_vector)
    cell = best_sol.cell_id

    # Re-run the same param set
    setup = attempt_setup_community(
        idx[cell][1], idx[cell][2],
        best_sol.mu, best_sol.mu_predation, best_sol.epsilon_val, best_sol.sym_comp;
        localNPP=localNPP, localH0_vector=localH0_vector, species_names=sp_nm_all
    )
    if setup === nothing
        return nothing
    end

    (S2, R2, H_i0, m_i, g_i, G, M_mod, a_matrix, A, eps_vec, m_alpha,
     herb_list, pred_list) =
     (setup.S, setup.R, setup.H_i0, setup.m_i, setup.g_i, setup.G,
      setup.M_modified, setup.a_matrix, setup.A, setup.epsilon_vector, setup.m_alpha,
      setup.herbivore_list, setup.predator_list)

    H_init = copy(H_i0)
    for (i_h, h_name) in enumerate(herb_list)
        if h_name in removal_list
            H_init[i_h] = 0.0
        end
    end
    P_init = H_init[1:R2] ./ 10.0
    for (i_p, p_name) in enumerate(pred_list)
        if p_name in removal_list
            P_init[i_p] = 0.0
        end
    end

    u0 = vcat(H_init, P_init)
    params = (S2, R2, H_i0, m_i, g_i, G, M_mod, a_matrix, A, eps_vec, m_alpha)

    sim = run_simulation(u0, params)
    sr = sim.sr
    giHi = sum(g_i .* sim.H_end)
    ratio = giHi / localNPP
    ratio_ok = (ratio > 0.5) && (ratio < 10.0)

    if ratio_ok
        return (
            cell_id = cell,
            mu = best_sol.mu,
            mu_predation = best_sol.mu_predation,
            epsilon_val = best_sol.epsilon_val,
            sym_comp = best_sol.sym_comp,
            S2 = S2,
            R2 = R2,
            survival_rate = sr,
            survived_herbivores = sim.survived_herb,
            survived_predators  = sim.survived_pred,
            H_biomass = sim.H_bio,
            P_biomass = sim.P_bio,
            total_biomass = sim.total_bio,
            g_iH_i = giHi,
            g_iH_i_over_NPP = ratio,
            sp_removed = removal_list
        )
    else
        return nothing
    end
end

# --------------------------------------------------
# MAIN PARALLEL DRIVER
# --------------------------------------------------
function run_pipeline(cells_to_process, param_combinations)
    Threads.@threads for cell in cells_to_process
        @info "Processing cell $cell on thread $(Threads.threadid())..."
        local_i, local_j = idx[cell][1], idx[cell][2]
        sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi[local_i, local_j])
        localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
        localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)

        # Possibly remove no-prey predators
        predator_check = check_predator_has_prey(sp_nm)
        if !predator_check[1]
            filter!(s->!(s in predator_check[3]), sp_nm)
        end

        # --- (A) Find best param for all species
        best_sr_all, best_solutions_all =
            find_best_param_all_species(cell, param_combinations, sp_nm, localNPP, localH0_vector)

        if best_sr_all < 1e-12
            println(best_solutions_all)
            @warn "No feasible param found for cell $cell with sr = $best_sr_all."
            continue
        end

        found_1_all = any(sol->isapprox(sol.survival_rate,1.0; atol=1e-10), best_solutions_all)

        # Write them out: scenario_type="all_species", sp_removed="none"
        for sol in best_solutions_all
            rownt = (
                cell_id             = sol.cell_id,
                mu                  = sol.mu,
                mu_predation        = sol.mu_predation,
                epsilon_val         = sol.epsilon_val,
                sym_comp            = sol.sym_comp,
                S2                  = sol.S2,
                R2                  = sol.R2,
                survival_rate       = sol.survival_rate,
                survived_herbivores = sol.survived_herbivores,
                survived_predators  = sol.survived_predators,
                H_biomass           = sol.H_biomass,
                P_biomass           = sol.P_biomass,
                total_biomass       = sol.total_biomass,
                g_iH_i              = sol.g_iH_i,
                g_iH_i_over_NPP     = sol.g_iH_i_over_NPP,
                sp_removed          = "none",
                scenario_type       = "all_species"
            )
            write_row_to_csv(rownt)
        end

        # If any all-species scenario yields sr=1.0, skip removals
        if found_1_all
            @info "Cell $cell => full survival with all species & ratio_ok. Skipping removals."
            continue
        end

        # --- (B) Single removal
        param_check = attempt_setup_community(
            local_i, local_j,
            best_solutions_all[1].mu,
            best_solutions_all[1].mu_predation,
            best_solutions_all[1].epsilon_val,
            best_solutions_all[1].sym_comp;
            localNPP=localNPP, species_names=sp_nm,
            localH0_vector=localH0_vector
        )
        if param_check === nothing
            # could happen if something about re-setup fails
            continue
        end

        full_sp_list_herb = param_check.herbivore_list
        full_sp_list_pred = param_check.predator_list
        sp_ordered = vcat(full_sp_list_herb, full_sp_list_pred)

        best_sr_single = 0.0
        best_solutions_single = NamedTuple[]

        for s_combo in combinations(sp_ordered, 1)
            scenario = simulate_removal(best_solutions_all[1], collect(s_combo), sp_ordered, localNPP, localH0_vector)
            if scenario === nothing
                continue
            end
            sr = scenario.survival_rate

            if sr > best_sr_single
                best_sr_single = sr
                best_solutions_single = [scenario]
            elseif isapprox(sr, best_sr_single; atol=1e-10)
                push!(best_solutions_single, scenario)
            end
        end

        if best_sr_single > best_sr_all
            found_1_single = any(sol->isapprox(sol.survival_rate,1.0; atol=1e-10), best_solutions_single)
            for sol in best_solutions_single
                rownt = (
                    cell_id             = sol.cell_id,
                    mu                  = sol.mu,
                    mu_predation        = sol.mu_predation,
                    epsilon_val         = sol.epsilon_val,
                    sym_comp            = sol.sym_comp,
                    S2                  = sol.S2,
                    R2                  = sol.R2,
                    survival_rate       = sol.survival_rate,
                    survived_herbivores = sol.survived_herbivores,
                    survived_predators  = sol.survived_predators,
                    H_biomass           = sol.H_biomass,
                    P_biomass           = sol.P_biomass,
                    total_biomass       = sol.total_biomass,
                    g_iH_i              = sol.g_iH_i,
                    g_iH_i_over_NPP     = sol.g_iH_i_over_NPP,
                    sp_removed          = join(sol.sp_removed, ","),
                    scenario_type       = "single"
                )
                write_row_to_csv(rownt)
            end

            if found_1_single
                @info "Cell $cell => full survival found by single removal => skip pairs/triplets."
                continue
            end
        end

        # --- (C) Pair removal
        best_sr_pair = max(best_sr_all, best_sr_single)
        best_solutions_pair = NamedTuple[]

        for s_combo in combinations(sp_ordered, 2)
            scenario = simulate_removal(best_solutions_all[1], collect(s_combo), sp_ordered, localNPP, localH0_vector)
            if scenario === nothing
                continue
            end
            sr = scenario.survival_rate

            if sr > best_sr_pair
                best_sr_pair = sr
                best_solutions_pair = [scenario]
            elseif isapprox(sr, best_sr_pair; atol=1e-10)
                push!(best_solutions_pair, scenario)
            end
        end

        if best_sr_pair > max(best_sr_all,best_sr_single)
            found_1_pair = any(sol->isapprox(sol.survival_rate,1.0; atol=1e-10), best_solutions_pair)
            for sol in best_solutions_pair
                rownt = (
                    cell_id             = sol.cell_id,
                    mu                  = sol.mu,
                    mu_predation        = sol.mu_predation,
                    epsilon_val         = sol.epsilon_val,
                    sym_comp            = sol.sym_comp,
                    S2                  = sol.S2,
                    R2                  = sol.R2,
                    survival_rate       = sol.survival_rate,
                    survived_herbivores = sol.survived_herbivores,
                    survived_predators  = sol.survived_predators,
                    H_biomass           = sol.H_biomass,
                    P_biomass           = sol.P_biomass,
                    total_biomass       = sol.total_biomass,
                    g_iH_i              = sol.g_iH_i,
                    g_iH_i_over_NPP     = sol.g_iH_i_over_NPP,
                    sp_removed          = join(sol.sp_removed, ","),
                    scenario_type       = "pair"
                )
                write_row_to_csv(rownt)
            end

            if found_1_pair
                @info "Cell $cell => full survival found by pair removal => skip triplets."
                continue
            end
        end

        # --- Triplets (if desired) are commented out below ...
    end
end

# -----------------------------------------
# EXAMPLE USAGE
# -----------------------------------------
mu_vals = [0.5]
mu_pred_vals = range(0.0, 0.08, length=100)
eps_vals = range(0.5, 1.0, length=10)
sym_comp_vals = [true]

param_combos = [
    (mu, mp, eps, sc) for mu in mu_vals
                     for mp in mu_pred_vals
                     for eps in eps_vals
                     for sc in sym_comp_vals
]

cells_to_process = 1:8

@info "Starting pipeline..."
run_pipeline(cells_to_process, param_combos)
@info "All done. Results are in $(RESULTS_CSV)."

# -----------------------------------------