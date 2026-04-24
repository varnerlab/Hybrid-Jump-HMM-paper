# =============================================================================
# 01b-Fit-Residual-Marginals.jl
#
# Alternative-baseline fit: one JumpHMM marginal per non-market ticker, fit
# on the OLS residual series e_i(t) = G_i(t) - α_i - β_i · G_m(t) instead of
# the full growth-rate series G_i(t). Downstream composition with the fitted
# residual marginal proceeds with no variance correction (s = 1), because
# the draw is already idiosyncratic by construction.
#
# Used by the `compose_residual_jumphmm` baseline added to Composers.jl in
# response to peer-review point P1 (R2/R3: the "fit-to-residuals" alternative
# was not benchmarked in the original submission).
#
# Inputs:
#   data/universe.jld2          → tickers, growth_rates (T × N)
#   data/sim-calibration.jld2   → DataFrame with (α, β, R², σ_eps_real)
#
# Outputs:
#   data/marginals-residuals.jld2 → Dict{String, JumpHiddenMarkovModel}
#
# Notes:
#   JumpHMM's fit() expects a price series and converts to annualized excess
#   growth rates internally via excess_growth_rates(P; rf, dt) =
#   (1/dt) log(P_t/P_{t-1}) - rf. To fit the residual series while keeping
#   that interface unchanged, we invert the transform: set P_0 = 1.0 and
#   P_t = P_{t-1} · exp(e_t · dt) for t ≥ 1. Then excess_growth_rates on this
#   pseudo-price series reproduces the residuals exactly (with rf = 0).
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg           = load_config()
market_ticker = cfg["universe"]["market_ticker"]
rf            = Float64(cfg["hmm"]["risk_free_rate"])
N_states      = Int(cfg["hmm"]["N"])
ν             = Float64(cfg["hmm"]["nu"])
dt            = Float64(cfg["hmm"]["dt"])
@assert rf == 0.0 "residual pseudo-price inversion assumes rf = 0 at fit time"

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading universe and calibration for residual fitting..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
cd = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m        = G[:, market_idx]
T          = length(G_m)
@info "Residual-fit setup" n_tickers = nrow(calib) T = T dt = dt

# ── 2. Fit JumpHMM per ticker on OLS residuals ──────────────────────────────
cache = joinpath(_PATH_TO_DATA, "marginals-residuals.jld2")
if isfile(cache)
    @info "Residual marginals already cached at $cache — skipping fit."
    exit(0)
end

marginals_resid = Dict{String,JumpHiddenMarkovModel}()

for (i, row) in enumerate(eachrow(calib))
    ticker = row.ticker
    α, β   = row.alpha, row.beta
    j      = findfirst(==(ticker), tickers)
    j === nothing && error("ticker $ticker missing from universe")

    # residual series e_t = G_i(t) - α - β G_m(t), length T
    e = G[:, j] .- α .- β .* G_m

    # pseudo-prices P with P[1] = 1, P[t+1] = P[t] · exp(e[t] · dt) so that
    # excess_growth_rates(P; rf=0, dt) reproduces `e`.
    pseudo_prices = Vector{Float64}(undef, T + 1)
    pseudo_prices[1] = 1.0
    for t in 1:T
        pseudo_prices[t + 1] = pseudo_prices[t] * exp(e[t] * dt)
    end

    @info "Fitting residual marginal $i / $(nrow(calib)): $ticker" β = round(β, digits=3)
    marginals_resid[ticker] = fit(JumpHiddenMarkovModel, pseudo_prices;
                                   rf = rf, N = N_states, ν = ν, dt = dt)
end

@info "Caching residual marginals to $cache"
jldsave(cache; marginals = marginals_resid)
@info "Done — $(length(marginals_resid)) residual marginals fit and cached."
