# =============================================================================
# 03c-Seed-Sweep.jl
#
# Seed-to-seed uncertainty quantification (peer-review P1, R2): re-run the
# composer evaluation across the seeds listed in cfg["simulation"]["seed_list"]
# and persist each seed's results separately. Aggregation across seeds
# (means, SDs, 90% CIs on pass rates) is done in scripts/04-Tables.jl when it
# detects per-seed files.
#
# Each seed's run calls `run_composer_experiment(cfg; seed = s,
# output_suffix = "-seed-$s")` and writes to data/results-seed-$s.jld2 /
# results-seed-$s.csv. The base-seed run (s = cfg["simulation"]["seed"]) is
# also produced under the same suffix for uniformity; the default
# data/results.jld2 from scripts/03 is not overwritten.
#
# Runtime: roughly n_seeds × (scripts/03 runtime). With 8 seeds and 6
# composers, budget ~4 h; parallelization is straightforward if needed.
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg        = load_config()
seed_list  = get(cfg["simulation"], "seed_list", [Int(cfg["simulation"]["seed"])])
@info "Seed sweep" n_seeds = length(seed_list) seeds = seed_list

for s in seed_list
    @info ">>> Starting seed sweep leg: seed = $s"
    run_composer_experiment(cfg; seed = Int(s), output_suffix = "-seed-$(s)")
end

@info "Seed sweep complete. Aggregation: scripts/04-Tables.jl will pick up per-seed files."
