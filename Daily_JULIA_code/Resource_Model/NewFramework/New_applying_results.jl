using CSV, DataFrames, Evolutionary, DifferentialEquations, Random, Logging, Serialization
using Base.Threads: SpinLock, @threads

# ===========================
# GLOBAL CONSTANTS & LOCKS
# ===========================
const EXTINCTION_THRESHOLD = 1e-6
const art_pi = false      # adjust as needed
const file_lock = SpinLock()

# ---------------------------
# Thread-safe file writing
# ---------------------------
function write_result_jls(result, filename::String)
    lock(file_lock)
    try
        if !isfile(filename)
            open(filename, "w") do io
                serialize(io, [result])
            end
        else
            open(filename, "r") do io
                data = deserialize(io)
                push!(data, result)
                open(filename, "w") do io2
                    serialize(io2, data)
                end
            end
        end
    finally
        unlock(file_lock)
    end
end

# ---------------------------
# Collateral extinctions helper
# ---------------------------
function collateral_extinctions(active_removed::Vector{Int}, full_mask::Vector{Int}, current_mask::Vector{Int}, names::Vector{String})
    inds = [i for i in 1:length(names) if full_mask[i] == 0 && current_mask[i] == 1 && !(i in active_removed)]
    if isempty(inds)
        return ("none", ["none"])
    else
        return (join(string.(inds), ","), [names[i] for i in inds])
    end
end

# =============================================================================
# Helper: Compute extinct indices from equilibrium arrays
# =============================================================================
# We assume the full species list is ordered as all herbivores followed by all predators.
function extinct_indices(H_eq::Vector{Float64}, P_eq::Vector{Float64})
    extinct_H = [i for i in 1:length(H_eq) if H_eq[i] < EXTINCTION_THRESHOLD]
    extinct_P = [i for i in 1:length(P_eq) if P_eq[i] < EXTINCTION_THRESHOLD]
    return vcat(extinct_H, [length(H_eq) + i for i in extinct_P])
end

# =============================================================================
# Modified removal scenario that accepts a removal mask (a vector of indices to force to zero)
# =============================================================================
function run_removal_scenario_with_mask(cell::Int, sp_list::Vector{String},
    S_full::Int, R_full::Int, H_init_full::Vector{Float64}, P_init_full::Vector{Float64},
    baseline_params, removal_mask::Vector{Int};
    callbacks = false)
    
    sp_removed_names = sp_list[removal_mask]  # species forced to zero
    H_init = copy(H_init_full)
    P_init = copy(P_init_full)
    for i in removal_mask
        if i <= S_full
            H_init[i] = 0.0
        else
            pred_idx = i - S_full
            P_init[pred_idx] = 0.0
        end
    end
    
    u0_removal = vcat(H_init, P_init)
    prob_removal = ODEProblem(new_dynamics!, u0_removal, (0.0, 500.0), baseline_params)
    logger = SimpleLogger(stderr, Logging.Error)
    sol_removal = with_logger(logger) do
        if callbacks
            solve(prob_removal, Tsit5(); callback=cb_no_trigger, abstol=1e-8, reltol=1e-6)
        else
            solve(prob_removal, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
    end
    if sol_removal.t[end] < 500.0 || any(isnan, sol_removal.u[end]) || any(isinf, sol_removal.u[end])
        @warn "Removal simulation failed for cell $cell because <500 or infs detected."
        return nothing
    end

    H_end_rem = copy(sol_removal[1:S_full, end])
    P_end_rem = (R_full > 0) ? copy(sol_removal[S_full+1:S_full+R_full, end]) : Float64[]
    H_end_rem[H_end_rem .< EXTINCTION_THRESHOLD] .= 0.0
    if R_full > 0
        P_end_rem[P_end_rem .< EXTINCTION_THRESHOLD] .= 0.0
    end
    survived_herb_rem = count(H_end_rem .> 0.0)
    survived_pred_rem = (R_full > 0) ? count(P_end_rem .> 0.0) : 0
    sr_rem = (S_full + R_full) > 0 ? (survived_herb_rem + survived_pred_rem) / (S_full + R_full) : 0.0
    total_bio_rem = sum(H_end_rem) + sum(P_end_rem)
    return (
        S_rem = S_full,
        R_rem = R_full,
        H_end = H_end_rem,
        P_end = P_end_rem,
        survival_rate = sr_rem,
        total_biomass = total_bio_rem,
        survived_herb = survived_herb_rem,
        survived_pred = survived_pred_rem,
        sp_removed = sp_removed_names
    )
end

# =============================================================================
# Modified Baseline and Removal Experiment for a Cell (NEW VERSION)
# =============================================================================
function run_experiments_for_cell_new(cell::Int, best_row::DataFrameRow; callbacks = false)
    results_df = DataFrame()
    
    # 1) Get cell indices and full species list (all species; no filtering here)
    local_i, local_j = idx[cell][1], idx[cell][2]
    sp_full = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
    # 2) Run new_attempt_setup_community with ALL species.
    mu_val      = best_row.mu
    mu_pred_val = best_row.mu_predation
    eps_val     = best_row.epsilon_val
    m_alpha_val = best_row.m_alpha
    localNPP    = best_row.localNPP
    alpha_val   = 0.25
    hollingII   = true
    h_val       = 0.1

    setup_results = new_attempt_setup_community(
        local_i, local_j,
        mu_val, mu_pred_val, eps_val, true, m_alpha_val;
        localNPP = localNPP,
        species_names = sp_full,
        artificial_pi = true,
        alpha = alpha_val,
        hollingII = hollingII,
        h = h_val
    )
    if setup_results === nothing
        @warn "Cell $cell: Setup failed using full species list."
        return results_df
    end

    # 3) Extract parameters from setup_results.
    S_new = setup_results.S
    R_new = setup_results.R
    H_i0 = setup_results.H_i0
    m_i = setup_results.m_i
    g_i = setup_results.g_i
    beta = setup_results.beta
    M_mod = setup_results.M_modified
    A_star = setup_results.A_star
    a_matrix = setup_results.a_matrix
    A = setup_results.A
    epsilon_vector = setup_results.epsilon
    m_alpha_vec = setup_results.m_alpha
    A_pred = transpose(a_matrix)
    P0 = m_alpha_vec
    h_used = setup_results.h

    # 4) Compute the equilibrium (baseline) using nlsolve.
    x0 = vcat(H_i0, fill(0.1, R_new))
    params_tuple = (S_new, R_new, H_i0, m_i, g_i, beta, M_mod, A_star, A_pred, P0, A, m_alpha_vec, h_used)
    sol_eq = nlsolve((F, x) -> equilibrium_system!(F, x, params_tuple), x0)
    if !sol_eq.f_converged
        @warn "Cell $cell: Equilibrium solver did not converge."
        return results_df
    end

    H_eq = copy(sol_eq.zero[1:S_new])
    P_eq = (R_new > 0) ? copy(sol_eq.zero[S_new+1:S_new+R_new]) : Float64[]
    H_eq[H_eq .< EXTINCTION_THRESHOLD] .= 0.0
    if R_new > 0
        P_eq[P_eq .< EXTINCTION_THRESHOLD] .= 0.0
    end

    # 5) Run baseline simulation using these equilibrium values as initial conditions.
    H_init = H_eq
    P_init = (R_new > 0) ? P_eq : Float64[]
    u0 = vcat(H_init, P_init)
    prob = ODEProblem(new_dynamics!, u0, (0.0, 500.0), params_tuple)
    logger = SimpleLogger(stderr, Logging.Error)
    sol = with_logger(logger) do
        if callbacks
            solve(prob, Tsit5(); callback=cb_no_trigger, abstol=1e-8, reltol=1e-6)
        else
            solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
    end

    if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
        @warn "Cell $cell: Baseline simulation failed."
        return results_df
    end

    H_end = copy(sol[1:S_new, end])
    P_end = (R_new > 0) ? copy(sol[S_new+1:S_new+R_new, end]) : Float64[]
    H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
    if R_new > 0
        P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
    end

    survived_herb = count(H_end .> 0.0)
    survived_pred = (R_new > 0) ? count(P_end .> 0.0) : 0
    total_surv = survived_herb + survived_pred
    survival_rate = (S_new + R_new) > 0 ? total_surv / (S_new + R_new) : 0.0
    total_bio = sum(H_end) + sum(P_end)

    # 6) Save baseline results.
    push!(results_df, (
        cell = cell,
        sp_removed = "none",
        sp_id_removed = "0",
        survival_rate = round(survival_rate, digits=3),
        total_biomass = round(total_bio, digits=3),
        h_biomass = round(sum(H_end), digits=3),
        p_biomass = round(sum(P_end), digits=3),
        herb_pred_ratio = (sum(H_end)==0.0) ? NaN : round(sum(P_end)/sum(H_end), digits=3),
        herbivore_survival_rate = (S_new > 0) ? round(survived_herb/S_new, digits=3) : 0.0,
        predator_survival_rate  = (R_new > 0) ? round(survived_pred/R_new, digits=3) : 0.0,
        delta_total_biomass = 0.0,
        H_biomass_vector = H_end,
        P_biomass_vector = P_end,
        H_full_minus_H = fill(0.0, S_new),
        P_full_minus_P = fill(0.0, R_new),
        ind_ext_num = "none",
        ind_ext_name = ["none"],
        survived_herbs = survived_herb,
        survived_preds = survived_pred
    ))

    # 7) Determine indices to always remove (those already extinct).
    extinct_mask = extinct_indices(H_eq, P_eq)
    @info "Cell $cell: Already extinct indices: $(extinct_mask)"
    candidate_indices = setdiff(1:(S_new+R_new), extinct_mask)

    # 8) For each candidate additional removal, form the removal mask as the union of extinct indices and that candidate.
    @info "Candidate species for removal in cell $cell: $(candidate_indices)"

    for cand in candidate_indices
        # Ensure only extinct species + the one selected for removal are affected
        removal_mask = vcat(extinct_mask, cand)
        
        removal_out = run_removal_scenario_with_mask(
            cell,
            sp_full,   # full species list (order as in the cell)
            S_new, R_new,
            H_init, P_init,
            params_tuple,
            removal_mask;
            callbacks = callbacks
        )
    
        if removal_out === nothing
            @warn "Cell $cell: Removal scenario for species $(cand) failed. Skipping."
            continue
        end
    
        H_end_rem = removal_out.H_end
        P_end_rem = removal_out.P_end
        sr_rem = removal_out.survival_rate
        total_bio_rem = removal_out.total_biomass
        sp_removed_names = removal_out.sp_removed
    
        full_ext_mask = map(x -> x > 0.0 ? 0 : 1, vcat(H_end, P_end))
        current_ext_mask = map(x -> x > 0.0 ? 0 : 1, vcat(H_end_rem, P_end_rem))
        (ind_ext_num_str, ind_ext_name_str) = collateral_extinctions(removal_mask, full_ext_mask, current_ext_mask, sp_full)
    
        H_diff = H_end .- H_end_rem
        P_diff = P_end .- P_end_rem
        delta_total_biomass = total_bio - total_bio_rem
    
        push!(results_df, (
            cell = cell,
            sp_removed = join(string.(sp_removed_names), ";"),
            sp_id_removed = string(cand),  # Now correctly records one species per row
            survival_rate = round(sr_rem, digits=3),
            total_biomass = round(total_bio_rem, digits=3),
            h_biomass = round(sum(H_end_rem), digits=3),
            p_biomass = round(sum(P_end_rem), digits=3),
            herb_pred_ratio = (sum(H_end_rem)==0.0) ? NaN : round(sum(P_end_rem)/sum(H_end_rem), digits=3),
            herbivore_survival_rate = (S_new > 0) ? round(removal_out.survived_herb/S_new, digits=3) : 0.0,
            predator_survival_rate  = (R_new > 0) ? round(removal_out.survived_pred/R_new, digits=3) : 0.0,
            delta_total_biomass = round(delta_total_biomass, digits=3),
            H_biomass_vector = H_end_rem,
            P_biomass_vector = P_end_rem,
            H_full_minus_H = H_diff,
            P_full_minus_P = P_diff,
            ind_ext_num = ind_ext_num_str,
            ind_ext_name = ind_ext_name_str,
            survived_herbs = removal_out.survived_herb,
            survived_preds = removal_out.survived_pred
        ))
    end    

    return results_df
end

# =============================================================================
# DRIVER FUNCTION FOR ALL CELLS
# =============================================================================
function run_keystone_removal(new_config_df::DataFrame; jls_filename="all_results.jls", callbacks = false)
    all_results_list = Vector{DataFrame}()
    grouped = groupby(new_config_df, :cell_id)
    for subdf in grouped
        best_row = subdf[argmax(subdf.survival_rate), :]
        cell = best_row.cell_id
        if callbacks
            @info "Making callbacks for cell $cell"
            local_i, local_j = idx[cell][1], idx[cell][2]
            # Extract species names from the cell data.
            sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
            local_S, local_R = identify_n_of_herbs_and_preds(sp_nm)
            cb_no_trigger, cb_trigger = build_callbacks(local_S, local_R, EXTINCTION_THRESHOLD, T_ext, 1)
        end
        @info "Processing cell $cell using best configuration"
        cell_results = run_experiments_for_cell_new(cell, best_row; callbacks = callbacks)
        if cell_results !== nothing
            push!(all_results_list, cell_results)
            write_result_jls(cell_results, jls_filename)
        end
    end
    return all_results_list
end

# =============================================================================
# DRIVER CALL
# =============================================================================
A_results = run_keystone_removal(A_hol; jls_filename="Daily_Julia_code/Resource_Model/Best_params_&_other_outputs/NEW_keystone_results.jls")

AA_results = A_results[1]#[:, 2][2]