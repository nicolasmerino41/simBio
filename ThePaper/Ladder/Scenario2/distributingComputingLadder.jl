# ───────────────────────────────────────────────────────────────────────────────
# (0) HOW MANY NODES / THREADS?
# ───────────────────────────────────────────────────────────────────────────────
using Distributed
using ClusterManagers
using CSV, DataFrames, Logging, Random
using Base.Threads: SpinLock, @threads
num_nodes = 16
num_threads = 2
queue = "long"
# -------------------------------------------------------------------
# 0) Launching workers: Suppose we want 16 total processes
#    (one per node?), each with 48 threads:
#    *** If you only requested 1 process per node in your Slurm job,
#    change SlurmManager(16) to SlurmManager(1) etc. as needed. ***
# -------------------------------------------------------------------
addprocs(SlurmManager(num_nodes), exeflags=["--project", "--threads=$num_threads"])

# ───────────────────────────────────────────────────────────────────────────────
# (1) BROADCAST YOUR ENTIRE PIPELINE TO ALL WORKERS
# ───────────────────────────────────────────────────────────────────────────────
@everywhere begin
  # exactly as you had at top of your script:
  using DifferentialEquations, Random, LinearAlgebra, Statistics, DataFrames, Graphs
  import Base.Threads: @threads

  include("Ladder4.1.jl")            # your existing file
  include("ComputingLadder.jl") # file containing your ComputingLadder(...) function

  # we’ll write out to this CSV:
  const OUTPUT_CSV = "Final_results_distributed.csv"

  # thread‐safe lock for CSV appends
  const file_lock = SpinLock()

  function write_row!(row::NamedTuple)
    lock(file_lock)
    try
      df = DataFrame([row])
      CSV.write(OUTPUT_CSV, df; append = isfile(OUTPUT_CSV))
    finally
      unlock(file_lock)
    end
  end
end

# ───────────────────────────────────────────────────────────────────────────────
# (2) EXTRACT ALL THE PARAMETER COMBINATIONS
#     (exact same ranges you were using inside ComputingLadder)
# ───────────────────────────────────────────────────────────────────────────────
# const PARAM_COMBOS = collect(Iterators.product(
#   0.01:0.04:0.9,              # conn_vals
#   [0.01,0.1,1.0,2.0],         # IS_vals
#   [:ER,:PL,:MOD],             # scenarios
#   [0.1,0.3,0.5,2.0],          # delta_vals
#   [1.0,0.1,0.5],              # eps_scales
#   [0.1,0.2,0.3,0.4,0.5],      # mortality_vals
#   [0.5,1.0,3.0,5.0,7.0]       # growth_vals
# ))
# shuffle!(PARAM_COMBOS)

PARAM_COMBOS = collect(Iterators.product(
  0.01:0.04:0.9,              # conn_vals
  [0.01],         # IS_vals
  [:ER],             # scenarios
  [0.1],          # delta_vals
  [1.0],              # eps_scales
  [0.1],      # mortality_vals
  [0.5]       # growth_vals
))
shuffle!(PARAM_COMBOS)

# split into one “chunk” per worker process
nworkers = nprocs() - 1
const CHUNKS = [ PARAM_COMBOS[i:nworkers:end] for i in 1:nworkers ]

# ───────────────────────────────────────────────────────────────────────────────
# (3) EACH WORKER THREADS OVER ITS CHUNK, CALLS YOUR ComputingLadder, WRITES OUT
# ───────────────────────────────────────────────────────────────────────────────
@everywhere function process_chunk(chunk)
  for (conn, IS, scen, delta, eps, m_val, g_val) in chunk
    @threads for _ in 1:1
      # call exactly your function, with the same keyword argument names:
      row = ComputingLadder(
        50, 20;
        conn_vals    = [conn],
        IS_vals      = [IS],
        scenarios    = [scen],
        delta_vals   = [delta],
        eps_scales   = [eps],
        tspan        = (0.,500.),
        tpert        = 250.0,
        number_of_combinations = 1,
        threshold_steps = 3,
        B_term       = false,
        mortality_vals = [m_val],
        growth_vals    = [g_val]
      )

      # your ComputingLadder returns a DataFrame; we expect exactly one row
      if nrow(row) == 1
        write_row!(row[1, :])   # append that row
      elseif nrow(row) > 1
        @warn "Chunk produced >1 rows—writing them all"
        for r in eachrow(row)
          write_row!(NamedTuple(r))
        end
      else
        # infeasible / no output for this combo
      end
    end
  end
  return nothing
end

# ───────────────────────────────────────────────────────────────────────────────
# (4) FIRE OFF pmap—ONE CHUNK PER WORKER
# ───────────────────────────────────────────────────────────────────────────────
pmap(process_chunk, CHUNKS)

@info "All done!  See $(first(CHUNKS))  chunks processed → $(OUTPUT_CSV)"
