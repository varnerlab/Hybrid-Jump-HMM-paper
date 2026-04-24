# =============================================================================
# 03-Compose-And-Evaluate.jl
#
# Thin wrapper around `run_composer_experiment` (src/Pipeline.jl) that runs
# the single-seed default evaluation: all six composers, real SPY market
# path, floor and threshold from config.toml. Outputs to
# data/results.jld2 and data/results-summary.csv.
#
# Composers 4 and 6 (residual_jumphmm, garch_t) are skipped automatically if
# their fit caches are absent, so the script runs cleanly at any stage.
#
# For multi-seed uncertainty quantification see scripts/03c-Seed-Sweep.jl.
# For sensitivity sweeps (R² threshold, floor f) see scripts/03d-Sensitivity-Sweep.jl.
# For the clipping stress test see scripts/03b-Stress-Eval.jl.
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg = load_config()
run_composer_experiment(cfg)
