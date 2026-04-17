# =============================================================================
# 01-Fit-Marginals.jl
#
# Loads the JDIQ training universe (424 US equities/ETFs) via
# VLQuantitativeFinancePackage, filters to tickers with the required minimum
# observation count, and fits a per-ticker JumpHMM marginal.
#
# Outputs:
#   data/universe.jld2     → (tickers, prices, growth_rates)
#   data/marginals.jld2    → Dict{String, JumpHiddenMarkovModel}
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg = load_config()
min_obs = Int(cfg["universe"]["min_obs_required"])
rf      = Float64(cfg["hmm"]["risk_free_rate"])
dt      = Float64(cfg["hmm"]["dt"])

# ── 1. Load and filter universe ─────────────────────────────────────────────
@info "Loading universe (min_obs = $min_obs days)..."
tickers, prices = load_universe(min_obs)

# ── 2. Compute growth-rate matrix ───────────────────────────────────────────
@info "Computing annualized excess growth rates..."
G = growth_rate_matrix(prices; rf = rf, dt = dt)
@info "Growth-rate matrix" size = size(G)

universe_path = joinpath(_PATH_TO_DATA, "universe.jld2")
@info "Caching universe to $universe_path"
jldsave(universe_path; tickers = tickers, prices = prices, growth_rates = G)

# ── 3. Fit per-ticker marginals ─────────────────────────────────────────────
@info "Fitting per-ticker JumpHMM marginals..."
marginals = fit_per_ticker_marginals(prices, tickers; cfg = cfg)
@info "Fit $(length(marginals)) per-ticker marginals"
