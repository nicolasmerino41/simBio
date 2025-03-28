include("Scripts/prior.jl")
include("Scripts/prior2.jl")
include("Scripts/DA_birmmals_with_pi.jl")
include("Scripts/generate_competition_matrix.jl")
include("Scripts/species_dict.jl")

include("Scripts/ecosystem_dynamics!.jl")
include("Scripts/FI_functions.jl")
include("Scripts/extract_H0_DA.jl")
include("Scripts/attempt_setup_community.jl")
include("Scripts/Callbacks_function.jl")
include("Scripts/npp_DA_relative_to_1000.jl")
include("Scripts/attempt_feasibility.jl")

#### We needed to loaded this way cause there was column problems but it works fine
lines = readlines("Results/Big_pipeline_results_(special)_not_even_pi.csv")
data = [join(split(line, ",")[1:24], ",") for line in lines if length(split(line, ",")) == 24]
# Now parse from an IOBuffer.
Big_P_results = CSV.read(IOBuffer(join(data, "\n")), DataFrame, header=true)
Big_P_results = Big_P_results[.!ismissing.(Big_P_results.cell_id), :][:, 1:24]
unique_cells_pi = string.(unique(Big_P_results.cell_id))

#### We needed to loaded this way cause there was column problems but it works fine
ga_cell_results = CSV.read("C:/Users/MM-1/Downloads/ga_cell_results.csv", DataFrame, header=false)

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

# Import SpinLock for thread-safety
using Base.Threads: SpinLock
using Base.Threads: @threads
# NOW WE'LL IMPLEMENT THE SAME APPROACH WE USED FOR Searching_keystone_species.jl AND SpeciesEffectOnEcosystemFunctioning.jl
const file_lock = ReentrantLock()
    
    """
        write_result_jls(result::Any, filename::String)
    
    Thread-safe function to append `result` to a list stored in `filename` as a `.jls` file.
    
    - Locks `file_lock` before reading/writing to avoid race conditions.
    - Reads existing data (if any) from `filename`, pushes `result`, then overwrites the file.
    """
    function write_result_jls(result, filename::String)
        lock(file_lock)
        try
            if !isfile(filename)
                # File doesn't exist. Create new array with this first result
                open(filename, "w") do io
                    serialize(io, [result])
                end
            else
                # File already exists: read data, push new result, re-serialize
                existing_data = nothing
                open(filename, "r") do io
                    existing_data = deserialize(io)
                end
                
                push!(existing_data, result)  # append the new row or DataFrame
    
                # Overwrite the file with the updated data
                open(filename, "w") do io
                    serialize(io, existing_data)
                end
            end
        finally
            unlock(file_lock)
        end
    end
    
    ###############################################################################
    # 2) COLLATERAL EXTINCTIONS
    ###############################################################################
    function collateral_extinctions(
        active_removed::Int,
        full_mask::Vector{Int},
        current_mask::Vector{Int},
        names::Vector{String}
    )
        # "Newly extinct" => was alive in the baseline (full_mask[i]==0)
        # but now extinct (current_mask[i]==1), excluding actively_removed index
        inds = [i for i in 1:length(names) if full_mask[i]==0 && current_mask[i]==1 && i!=active_removed]
        if isempty(inds)
            return ("none", ["none"])
        else
            return (join(string.(inds), ","), [names[i] for i in inds])
        end
    end
    
    ###############################################################################
    # 3) BASELINE SCENARIO
    ###############################################################################
    """
        run_baseline_scenario(cell, mu_val, mu_pred_val, eps_val, sym_competition, sp_removed_best)
    
    Runs the "full" scenario for a cell, except that if `sp_removed_best` is not `nothing`,
    that species is *entirely removed from the system*. 
    Returns a NamedTuple with final arrays, or `nothing` if the solver fails.
    """
    function run_baseline_scenario(
        cell::Int,
        mu_val::Float64,
        mu_pred_val::Float64,
        eps_val::Float64,
        sym_competition::Bool,
        sp_removed_best::Union{Nothing,AbstractString}
    )
        # 1) Identify i,j
        local_i, local_j = idx[cell][1], idx[cell][2]
    
        # 2) Build species list
        sp_nm = extract_species_names_from_a_cell(DA_birmmals_with_pi_corrected[local_i, local_j])
    
        # Possibly remove predators with no prey
        predator_has_prey = check_predator_has_prey(sp_nm)
        if !predator_has_prey[1]
            filter!(name -> !(name in predator_has_prey[3]), sp_nm)
        end
    
        # Exclude sp_removed_best if specified
        if !isnothing(sp_removed_best)
            filter!(name -> name != sp_removed_best, sp_nm)
        end
    
        # Compute final S, R
        (S_full, R_full) = identify_n_of_herbs_and_preds(sp_nm)
        if (S_full + R_full) == 0
            return nothing
        end
    
        # 3) NPP & H0
        localNPP       = Float64(npp_DA_relative_to_1000[local_i, local_j])
        localH0_vector = Vector{Float64}(H0_DA[local_i, local_j].a)
        
        # 4) Attempt setup
        setup_results = attempt_setup_community(
            local_i, local_j,
            mu_val, mu_pred_val, eps_val, sym_competition;
            localNPP      = localNPP,
            localH0_vector= localH0_vector,
            species_names = sp_nm,
            artificial_pi = true
        )
        if setup_results === nothing
            return nothing
        end
    
        # Destructure from setup
        (S_new, R_new, H_i0, m_i, g_i, G, M_modified,
         a_matrix, A, epsilon_vector, m_alpha) = (
            setup_results.S, setup_results.R, setup_results.H_i0, setup_results.m_i,
            setup_results.g_i, setup_results.G, setup_results.M_modified,
            setup_results.a_matrix, setup_results.A, setup_results.epsilon_vector,
            setup_results.m_alpha
        )
    
        if (S_new + R_new) == 0 || R_new > length(H_i0)
            return nothing
        end
    
        # 5) Baseline initial condition
        H_init = H_i0
        P_init = H_init[1:R_new] ./ 10.0
        u0     = vcat(H_init, P_init)
    
        # 6) Build baseline params
        baseline_params = (S_new, R_new, H_i0, m_i, g_i, G,
                           M_modified, a_matrix, A,
                           epsilon_vector, m_alpha)
    
        # 7) Solve the ODE
        prob   = ODEProblem(ecosystem_dynamics!, u0, (0.0, 500.0), baseline_params)
        logger = SimpleLogger(stderr, Logging.Error)
        sol    = with_logger(logger) do
            solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
    
        if sol.t[end] < 500.0 || any(isnan, sol.u[end]) || any(isinf, sol.u[end])
            return nothing
        end
    
        # 8) Final states
        H_end = copy(sol[1:S_new, end])
        P_end = copy(sol[S_new+1:S_new+R_new, end])
        H_end[H_end .< EXTINCTION_THRESHOLD] .= 0.0
        P_end[P_end .< EXTINCTION_THRESHOLD] .= 0.0
    
        survived_herb = count(H_end .> 0.0)
        survived_pred = count(P_end .> 0.0)
        survival_rate = (S_new + R_new) > 0 ? (survived_herb + survived_pred)/(S_new + R_new) : 0.0
        total_bio     = sum(H_end) + sum(P_end)
    
        # Return many objects directly:
        return (
            sp_list = sp_nm,
            S_full  = S_new,
            R_full  = R_new,
            H_full  = H_end,
            P_full  = P_end,
            localNPP = localNPP,
            localH0_vector = localH0_vector,
            g_i = g_i,
            H_i0 = H_i0,
            P_init = P_init,
            H_end = H_end,
            P_end = P_end,
            survival_rate = survival_rate,
            total_biomass = total_bio,
            H_bio = sum(H_end),
            P_bio = sum(P_end),
            full_surv_herb = survived_herb,
            full_surv_pred = survived_pred,
            baseline_params = baseline_params
    
        )
    end
    
    ###############################################################################
    # 4) ONE‐BY‐ONE REMOVALS (Keeping dimension the same)
    ###############################################################################
    """
        run_removal_scenario_baselineDims(
            cell, baseline_out, i_sp
        )
    
    Given the baseline scenario for a cell (which already excluded `sp_removed_best` if any),
    we do NOT remove that species from the system. Instead, we set its initial abundance to 0
    in the ODE initial condition, so the dimension (S_full+R_full) remains the same.
    
    - baseline_out: NamedTuple from `run_baseline_scenario`, containing:
        :sp_list (Vector{String}) with dimension S_full+R_full
        :S_full, :R_full
        :H_end_full, :P_end_full
        :g_i, etc.
    - i_sp: the 1-based index of the species to remove (1..(S_full+R_full))
    
    Returns another NamedTuple with final states, or `nothing` if it fails.
    """
    function run_removal_scenario_baselineDims(
        cell::Int,
        # Baseline data
        sp_list::Vector{String},
        S_full::Int,
        R_full::Int,
        H_init_full::Vector{Float64},
        P_init_full::Vector{Float64},
        baseline_params,  # The ODE params from the baseline
        i_sp::Int
    )
        # Identify species being removed
        sp_removed = sp_list[i_sp]
    
        # Build new initial state by copying
        H_init = copy(H_init_full)
        P_init = copy(P_init_full)
    
        if i_sp <= S_full
            H_init[i_sp] = 0.0
        else
            pred_idx = i_sp - S_full
            P_init[pred_idx] = 0.0
        end
    
        # Combine
        u0_removal = vcat(H_init, P_init)
    
        # Solve
        prob_removal = ODEProblem(ecosystem_dynamics!, u0_removal, (0.0, 500.0), baseline_params)
        logger = SimpleLogger(stderr, Logging.Error)
        sol_removal = with_logger(logger) do
            solve(prob_removal, Tsit5(); abstol=1e-8, reltol=1e-6)
        end
    
        if sol_removal.t[end]<500.0 || any(isnan, sol_removal.u[end]) || any(isinf, sol_removal.u[end])
            return nothing
        end
    
        # Extract final states
        H_end_rem = copy(sol_removal[1:S_full, end])
        P_end_rem = copy(sol_removal[S_full+1:S_full+R_full, end])
        H_end_rem[H_end_rem .< EXTINCTION_THRESHOLD] .= 0.0
        P_end_rem[P_end_rem .< EXTINCTION_THRESHOLD] .= 0.0
    
        survived_herb_rem = count(H_end_rem .> 0.0)
        survived_pred_rem = count(P_end_rem .> 0.0)
        sr_rem = (S_full+R_full) > 0 ? (survived_herb_rem+survived_pred_rem)/(S_full+R_full) : 0.0
        total_bio_rem = sum(H_end_rem) + sum(P_end_rem)
    
        # Return multiple values as a NamedTuple or just a tuple. Let's do a NamedTuple for clarity here:
        return (
            S_rem         = S_full,
            R_rem         = R_full,
            H_end         = H_end_rem,
            P_end         = P_end_rem,
            survival_rate = sr_rem,
            total_biomass = total_bio_rem,
            survived_herb = survived_herb_rem,
            survived_pred = survived_pred_rem,
            sp_removed    = sp_removed
        )
    end
    
    ###############################################################################
    # 5) MAIN: RUN BASELINE & THEN REMOVE SPECIES ONE BY ONE
    ###############################################################################
    function run_experiments_for_cell(
        cell::Int,
        mu_val::Float64,
        mu_pred_val::Float64,
        eps_val::Float64,
        sym_comp::Bool,
        sp_removed_best::Union{Nothing,String}
    )
        results_df = DataFrame()
    
        # A) Baseline
        # We call run_baseline_scenario, which returns NOTHING or 17 objects:
        baseline_out = run_baseline_scenario(
            cell, mu_val, mu_pred_val, eps_val, sym_comp, sp_removed_best
        )
        if baseline_out === nothing
            @warn "Cell $cell: baseline scenario failed."
            return results_df
        end
        # println("Baseline_out is", typeof(baseline_out))
        # Destructure the 17 objects:
        sp_list         = baseline_out.sp_list
        S_full          = baseline_out.S_full
        R_full          = baseline_out.R_full
        localNPP        = baseline_out.localNPP
        localH0_vector  = baseline_out.localH0_vector
        g_i             = baseline_out.g_i
        H_init_full     = baseline_out.H_i0
        P_init_full     = baseline_out.P_init
        H_end_full      = baseline_out.H_end
        P_end_full      = baseline_out.P_end
        full_sr         = baseline_out.survival_rate
        full_bio        = baseline_out.total_biomass
        H_bio           = baseline_out.H_bio
        P_bio           = baseline_out.P_bio
        full_surv_herb  = baseline_out.full_surv_herb
        full_surv_pred  = baseline_out.full_surv_pred
        baseline_params = baseline_out.baseline_params
    
        # Build a mask for baseline
        full_ext_mask = map(x -> x>0.0 ? 0 : 1, vcat(H_end_full, P_end_full))
    
        # Now push the "none removed" scenario
        push!(results_df, (
            cell = cell,
            sp_removed = "none",
            sp_id_removed = 0,
            survival_rate = round(full_sr, digits=3),
            total_biomass = round(full_bio, digits=3),
            h_biomass = round(H_bio, digits=3),
            p_biomass = round(P_bio, digits=3),
            herb_pred_ratio = (H_bio==0.0) ? NaN : round(P_bio/H_bio, digits=3),
            herbivore_survival_rate = (S_full>0) ? round(full_surv_herb/S_full, digits=3) : 0.0,
            predator_survival_rate  = (R_full>0) ? round(full_surv_pred/R_full, digits=3) : 0.0,
            delta_total_biomass = 0.0,
            H_biomass_vector = H_end_full,
            P_biomass_vector = P_end_full,
            H_full_minus_H = fill(0.0, S_full),
            P_full_minus_P = fill(0.0, R_full),
            ind_ext_num = "none",
            ind_ext_name = ["none"],
            survived_herbs = full_surv_herb,
            survived_preds = full_surv_pred
        ))
    
        # B) Remove each species one by one
        for i_sp in 1:(S_full + R_full)
            # If sp_removed_best was originally excluded, it won't appear in sp_list anyway
            if !isnothing(sp_removed_best) && sp_list[i_sp] == sp_removed_best
                continue
            end
    
            # Run removal scenario
            removal_out = run_removal_scenario_baselineDims(
                cell,
                sp_list,
                S_full,
                R_full,
                H_init_full,
                P_init_full,
                baseline_params,
                i_sp
            )
            if removal_out === nothing
                @warn "Cell $cell: removal scenario for sp=$(sp_list[i_sp]) failed. Skipping."
                continue
            end
    
            # Unpack
            S_rem         = removal_out.S_rem
            R_rem         = removal_out.R_rem
            H_end_rem     = removal_out.H_end
            P_end_rem     = removal_out.P_end
            sr_rem        = removal_out.survival_rate
            total_bio_rem = removal_out.total_biomass
            sp_removed    = removal_out.sp_removed
            surv_herb_rem = removal_out.survived_herb
            surv_pred_rem = removal_out.survived_pred
    
            # Build ext_mask for the removal scenario
            current_ext_mask = map(x -> x>0.0 ? 0 : 1, vcat(H_end_rem, P_end_rem))
            (ind_ext_num_str, ind_ext_name_str) = collateral_extinctions(
                i_sp, full_ext_mask, current_ext_mask, sp_list
            )
    
            # Differences
            H_full_minus_H = H_end_full .- H_end_rem
            P_full_minus_P = P_end_full .- P_end_rem
            delta_total_biomass = full_bio - total_bio_rem
    
            push!(results_df, (
                cell = cell,
                sp_removed = sp_removed,
                sp_id_removed = i_sp,
                survival_rate = round(sr_rem, digits=3),
                total_biomass = round(total_bio_rem, digits=3),
                h_biomass = round(sum(H_end_rem), digits=3),
                p_biomass = round(sum(P_end_rem), digits=3),
                herb_pred_ratio = (sum(H_end_rem)==0.0) ? NaN : round(sum(P_end_rem)/sum(H_end_rem), digits=3),
                herbivore_survival_rate = (S_rem>0) ? round(surv_herb_rem/S_rem, digits=3) : 0.0,
                predator_survival_rate  = (R_rem>0) ? round(surv_pred_rem/R_rem, digits=3) : 0.0,
                delta_total_biomass = round(delta_total_biomass, digits=3),
                H_biomass_vector = H_end_rem,
                P_biomass_vector = P_end_rem,
                H_full_minus_H = H_full_minus_H,
                P_full_minus_P = P_full_minus_P,
                ind_ext_num = ind_ext_num_str,
                ind_ext_name = ind_ext_name_str,
                survived_herbs = surv_herb_rem,
                survived_preds = surv_pred_rem
            ))
        end
    
        return results_df
    end
    
    
    ###############################################################################
    # 6) DRIVER FUNCTION FOR ALL CELLS
    ###############################################################################
    function run_keystone_removal(Big_P_results_maximised::DataFrame; jls_filename="all_results.jls")
        all_results_list = Vector{DataFrame}()
    
        Threads.@threads for row in eachrow(Big_P_results_maximised)
            cell            = row.cell_id
            mu_val          = row.mu
            mu_pred_val     = row.mu_predation
            eps_val         = row.epsilon_val
            sym_comp        = row.symmetrical_competition
            sp_removed_best = ismissing(row.sp_removed_name) ? nothing : String(row.sp_removed_name)
    
            @info "CELL=$cell => Baseline (excluding $sp_removed_best), then one-by-one removals"
            cell_df = run_experiments_for_cell(cell, mu_val, mu_pred_val, eps_val, sym_comp, sp_removed_best)
            
            @info "Finished cell $cell"
    
            # 1) Save the results of this cell to .jls
            write_result_jls(cell_df, jls_filename)
    
            # 2) Also push into in-memory structure if you still want it
            push!(all_results_list, cell_df)
        end
    
        return all_results_list
    end
    
    @time SPECIAL_all_results_list = run_keystone_removal(Big_P_results_maximised[1:10, :]; jls_filename="Results/SPECIAL_all_results_list_not_even_pi.jls")
    
    