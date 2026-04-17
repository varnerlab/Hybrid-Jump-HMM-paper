# =============================================================================
# 03-Compose-And-Evaluate.jl
#
# For each non-market ticker, draws n_paths replications from its JumpHMM
# marginal, composes a synthetic series via three composers (naive,
# Gaussian SIM, hybrid) over the same paired innovation, and scores each
# composed series against the real growth-rate history.
#
# The market path used for composition is the REAL SPY history (length T_eff).
# This isolates composer behavior from market-path randomness; a follow-up
# script can re-run with simulated market paths.
#
# Inputs:
#   data/universe.jld2          → tickers, growth_rates (T × N)
#   data/marginals.jld2         → Dict{String, JumpHiddenMarkovModel}
#   data/sim-calibration.jld2   → DataFrame
#
# Outputs:
#   data/results.jld2           → per-(ticker × composer × replication) metrics
#   data/results-summary.csv    → flat scoreboard for inspection
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg            = load_config()
market_ticker  = cfg["universe"]["market_ticker"]
n_paths        = Int(cfg["simulation"]["n_paths"])
seed           = Int(cfg["simulation"]["seed"])
f              = Float64(cfg["hybrid"]["idiosyncratic_floor"])
R²_thresh      = Float64(cfg["hybrid"]["r2_preserve_threshold"])

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading universe / marginals / calibration..."
ud   = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
md   = load(joinpath(_PATH_TO_DATA, "marginals.jld2"))
cd   = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
marginals = md["marginals"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m  = G[:, market_idx]
σ²_m = var(G_m)
T_eff = length(G_m)
@info "Composition setup" market = market_ticker T_eff = T_eff σ_m = sqrt(σ²_m) n_paths = n_paths

# ── 2. Sweep composers across the universe ──────────────────────────────────
Random.seed!(seed)

# pre-allocate result storage; each row = (ticker, composer, rep, ...metrics)
rows = NamedTuple[]

# helper: persist a metric NamedTuple alongside identifiers
function record!(rows, ticker::String, composer::String, rep::Int,
                 β_eff::Float64, flag::String,
                 metrics::NamedTuple)
    push!(rows, merge(
        (ticker = ticker, composer = composer, rep = rep,
         beta_eff = β_eff, flag = flag),
        metrics
    ))
end

n_assets = nrow(calib)
@info "Composing and scoring across $n_assets non-market tickers × $n_paths replications × 3 composers..."

for (i, row) in enumerate(eachrow(calib))
    ticker = row.ticker
    α      = row.alpha
    β      = row.beta
    R²     = row.r2_real
    σ_εr   = row.sigma_eps_real

    asset_idx = findfirst(==(ticker), tickers)
    G_real    = G[:, asset_idx]
    @assert length(G_real) == T_eff

    model = marginals[ticker]

    # one batch call: n_paths simulated paths of length T_eff
    sim_result = simulate(model, T_eff; n_paths = n_paths, seed = seed + i)

    if i % 10 == 0 || i == 1
        @info "  ticker $i / $n_assets: $ticker (β=$(round(β, digits=3)), R²=$(round(R², digits=3)))"
    end

    for r in 1:n_paths
        ε̃     = Float64.(sim_result.paths[r].observations)
        σ²_g  = var(ε̃)

        # --- naive ---
        g_n = compose_naive(α, β, G_m, ε̃)
        record!(rows, ticker, "naive", r, β, string(NAIVE),
                score_asset(g_n, G_real, G_m))

        # --- Gaussian SIM ---
        g_g = compose_gaussian_sim(α, β, σ_εr, G_m, Random.default_rng())
        record!(rows, ticker, "gaussian", r, β, string(GAUSSIAN_SIM),
                score_asset(g_g, G_real, G_m))

        # --- hybrid ---
        g_h, β_eff, flag = compose_hybrid(α, β, R², G_m, ε̃, σ²_m, σ²_g;
                                          f = f, R²_threshold = R²_thresh)
        record!(rows, ticker, "hybrid", r, β_eff, string(flag),
                score_asset(g_h, G_real, G_m))
    end
end

results = DataFrame(rows)
@info "Result rows" total = nrow(results) per_composer = nrow(results) ÷ 3

# ── 3. Persist ──────────────────────────────────────────────────────────────
results_path = joinpath(_PATH_TO_DATA, "results.jld2")
@info "Caching results to $results_path"
jldsave(results_path; results = results, config = cfg)

csv_path = joinpath(_PATH_TO_DATA, "results-summary.csv")
CSV.write(csv_path, results)
@info "Results CSV dumped to $csv_path"
