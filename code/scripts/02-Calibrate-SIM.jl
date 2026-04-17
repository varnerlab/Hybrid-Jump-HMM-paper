# =============================================================================
# 02-Calibrate-SIM.jl
#
# Runs OLS per ticker against the SPY market path, producing the calibrated
# (alpha, beta, R^2_real, sigma_eps_real) used as the target by the composers.
#
# Inputs:
#   data/universe.jld2  → tickers, growth_rates
#
# Outputs:
#   data/sim-calibration.jld2  → DataFrame(ticker, alpha, beta, r2_real, sigma_eps_real)
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg = load_config()
market_ticker = cfg["universe"]["market_ticker"]

universe_path = joinpath(_PATH_TO_DATA, "universe.jld2")
isfile(universe_path) || error("Run 01-Fit-Marginals.jl first; missing $universe_path")

@info "Loading universe artifact..."
ud = load(universe_path)
tickers = ud["tickers"]
G       = ud["growth_rates"]

@info "Calibrating SIM against $market_ticker for $(length(tickers) - 1) non-market tickers..."
calib = calibrate_sim(G, tickers, market_ticker)

# Quick sanity summary
@info "Calibration summary" n_assets = nrow(calib) min_R² = minimum(calib.r2_real) max_R² = maximum(calib.r2_real)
@info "β quartiles" q25 = quantile(calib.beta, 0.25) median = median(calib.beta) q75 = quantile(calib.beta, 0.75) max = maximum(calib.beta)

calib_path = joinpath(_PATH_TO_DATA, "sim-calibration.jld2")
@info "Caching calibration to $calib_path"
jldsave(calib_path; calibration = calib)

# Also dump a CSV for inspection
csv_path = joinpath(_PATH_TO_DATA, "sim-calibration.csv")
CSV.write(csv_path, calib)
@info "Calibration CSV dumped to $csv_path"
