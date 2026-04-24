# =============================================================================
# 06-VaR-Backtest.jl
#
# Downstream-use validation (peer-review P2, R3): per-ticker VaR backtest at
# α ∈ {0.95, 0.99} comparing composed paths against the real growth-rate
# history. Translates the paper's marginal-fidelity metric (KS pass rates)
# into a consequential VaR-accuracy number: how often does the real path
# breach the synthetic-ensemble VaR threshold, and is the exceedance rate
# statistically consistent with the nominal coverage level (Kupiec test)?
#
# The script pulls already-composed synthetic paths. It does NOT re-compose;
# if scripts/03 has not been run for the expanded composer set, run it first.
#
# Scoring loop: for each (ticker × composer × replication), compute one
# synthetic path's VaR threshold, count breaches on the real path, compute
# Kupiec p-value. Aggregate across replications with ticker- and
# composer-level median/mean/SE.
#
# Inputs:
#   data/universe.jld2          → real growth rates
#   data/marginals.jld2         → full-return marginals (for full composition)
#   data/marginals-residuals.jld2 → (optional) residual marginals
#   data/garch-t-sims.jld2      → (optional) GARCH-t sims
#   data/sim-calibration.jld2   → (α, β, R², σ_eps_real)
#
# Outputs:
#   data/var-backtest.jld2      → per-(ticker × composer × rep × α) records
#   data/var-backtest.csv
#   data/var-backtest-summary.csv  → aggregated per (composer, α)
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg            = load_config()
market_ticker  = cfg["universe"]["market_ticker"]
n_paths        = Int(cfg["simulation"]["n_paths"])
seed           = Int(cfg["simulation"]["seed"])
f              = Float64(cfg["hybrid"]["idiosyncratic_floor"])
R²_thresh      = Float64(cfg["hybrid"]["r2_preserve_threshold"])
block_length   = Float64(get(get(cfg, "bootstrap", Dict()), "mean_block_length", 50))
Δt             = Float64(cfg["hmm"]["dt"])
α_levels       = Float64[0.95, 0.99]

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading universe / marginals / calibration for VaR backtest..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
md = load(joinpath(_PATH_TO_DATA, "marginals.jld2"))
cd = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
marginals = md["marginals"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m        = G[:, market_idx]
σ²_m       = var(G_m)
T          = length(G_m)

residual_cache = joinpath(_PATH_TO_DATA, "marginals-residuals.jld2")
marginals_resid = isfile(residual_cache) ? load(residual_cache)["marginals"] : nothing

garch_cache = joinpath(_PATH_TO_DATA, "garch-t-sims.jld2")
garch_sims = isfile(garch_cache) ? load(garch_cache)["sims"] : nothing

@info "VaR backtest setup" T = T α_levels = α_levels residual = (marginals_resid !== nothing) garch = (garch_sims !== nothing)

# ── 2. Per-composer synthetic path generators ───────────────────────────────
# Each returns a one-day-return series of length T for a given (ticker, rep).

function synth_naive(ticker, r, sim_full)
    α, β = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta)
    ε̃ = Float64.(sim_full.paths[r].observations)
    return one_day_returns(compose_naive(α, β, G_m, ε̃), Δt)
end

function synth_gaussian(ticker, r, rng)
    α, β, σ_εr = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta),
                  calib_lookup(ticker, :sigma_eps_real)
    return one_day_returns(compose_gaussian_sim(α, β, σ_εr, G_m, rng), Δt)
end

function synth_hybrid(ticker, r, sim_full)
    α, β, R² = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta),
                calib_lookup(ticker, :r2_real)
    ε̃ = Float64.(sim_full.paths[r].observations)
    g, _, _ = compose_hybrid(α, β, R², G_m, ε̃, σ²_m, var(ε̃);
                              f = f, R²_threshold = R²_thresh)
    return one_day_returns(g, Δt)
end

function synth_residual(ticker, r, sim_resid)
    α, β = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta)
    ε̃ = Float64.(sim_resid.paths[r].observations)
    return one_day_returns(compose_residual_jumphmm(α, β, G_m, ε̃), Δt)
end

function synth_bootstrap(ticker, r, real_residuals, rng)
    α, β = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta)
    g = compose_block_bootstrap(α, β, G_m, real_residuals, block_length, rng)
    return one_day_returns(g, Δt)
end

function synth_garch(ticker, r)
    α, β = calib_lookup(ticker, :alpha), calib_lookup(ticker, :beta)
    ε̃ = Float64.(garch_sims[ticker][r])
    return one_day_returns(compose_garch_t(α, β, G_m, ε̃), Δt)
end

# lightweight per-ticker calibration lookup
let c = Dict{Tuple{String,Symbol},Float64}()
    for row in eachrow(calib)
        c[(row.ticker, :alpha)] = row.alpha
        c[(row.ticker, :beta)]  = row.beta
        c[(row.ticker, :r2_real)] = row.r2_real
        c[(row.ticker, :sigma_eps_real)] = row.sigma_eps_real
    end
    global function calib_lookup(t::AbstractString, k::Symbol)
        return c[(String(t), k)]
    end
end

# ── 3. Main loop ────────────────────────────────────────────────────────────
Random.seed!(seed)
rows = NamedTuple[]

function record!(rows, ticker, composer, rep, α_level, res)
    push!(rows, (ticker = ticker, composer = composer, rep = rep,
                 alpha_level = α_level, var = res.var, n_breach = res.n_breach,
                 rate = res.rate, expected_rate = res.expected_rate,
                 kupiec_p = res.kupiec_p))
end

for (i, row) in enumerate(eachrow(calib))
    ticker = row.ticker
    α, β   = row.alpha, row.beta
    j      = findfirst(==(ticker), tickers)
    G_real = G[:, j]
    r_real = one_day_returns(G_real, Δt)
    real_residuals = G_real .- α .- β .* G_m

    sim_full  = simulate(marginals[ticker], T; n_paths = n_paths, seed = seed + i)
    sim_resid = marginals_resid === nothing ? nothing :
        simulate(marginals_resid[ticker], T;
                 n_paths = n_paths, seed = seed + i + 1_000_000)

    if i % 25 == 0 || i == 1
        @info "  VaR ticker $i / $(nrow(calib)): $ticker"
    end

    for r in 1:n_paths
        for α_level in α_levels
            record!(rows, ticker, "naive", r, α_level,
                    var_backtest(synth_naive(ticker, r, sim_full), r_real, α_level))

            rng_g = MersenneTwister(seed + i * 2_000_000 + r * 2_000 + 1)
            record!(rows, ticker, "gaussian", r, α_level,
                    var_backtest(synth_gaussian(ticker, r, rng_g), r_real, α_level))

            record!(rows, ticker, "hybrid", r, α_level,
                    var_backtest(synth_hybrid(ticker, r, sim_full), r_real, α_level))

            if sim_resid !== nothing
                record!(rows, ticker, "residual_jumphmm", r, α_level,
                        var_backtest(synth_residual(ticker, r, sim_resid), r_real, α_level))
            end

            rng_b = MersenneTwister(seed + i * 3_000_000 + r * 3_000 + 1)
            record!(rows, ticker, "block_bootstrap", r, α_level,
                    var_backtest(synth_bootstrap(ticker, r, real_residuals, rng_b),
                                 r_real, α_level))

            if garch_sims !== nothing && haskey(garch_sims, ticker)
                record!(rows, ticker, "garch_t", r, α_level,
                        var_backtest(synth_garch(ticker, r), r_real, α_level))
            end
        end
    end
end

results = DataFrame(rows)
@info "VaR result rows" total = nrow(results)

# ── 4. Aggregated summary ───────────────────────────────────────────────────
summary = combine(groupby(results, [:composer, :alpha_level]),
                  :rate     => mean   => :mean_rate,
                  :rate     => std    => :sd_rate,
                  :kupiec_p => mean   => :mean_kupiec_p,
                  :kupiec_p => (x -> mean(x .> 0.05)) => :kupiec_pass_rate)

@info "VaR summary" summary

# ── 5. Persist ──────────────────────────────────────────────────────────────
jldsave(joinpath(_PATH_TO_DATA, "var-backtest.jld2"); results = results, config = cfg)
CSV.write(joinpath(_PATH_TO_DATA, "var-backtest.csv"), results)
CSV.write(joinpath(_PATH_TO_DATA, "var-backtest-summary.csv"), summary)
@info "VaR backtest persisted."
