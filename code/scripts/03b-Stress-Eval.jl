# =============================================================================
# 03b-Stress-Eval.jl
#
# Clipping-branch stress test. Re-runs the composer evaluation under market
# paths with inflated variance (factor × Var[SPY]) so the hybrid composer's
# HYBRID_CLIPPED branch is exercised on high-β tickers. Peer-review P1:
# R1 and R2 noted that clipping never triggered on the real SPY path and
# therefore the entire clipping subsection in the method was unvalidated.
#
# For each factor in cfg["stress"]["factors"], we scale the SPY growth-rate
# series about its mean (preserving temporal structure, only inflating
# variance), then run the same paired-innovation protocol as scripts/03 but
# restricted to the three composers that depend on the market variance
# (naive, gaussian_sim, hybrid). The other composers either resample real
# residuals (block_bootstrap) or sample from a residual-fit generator that
# has no knowledge of the market path, so they are unaffected by the stress.
#
# Inputs:
#   data/universe.jld2          → tickers, growth_rates (T × N)
#   data/marginals.jld2         → Dict{String, JumpHiddenMarkovModel}
#   data/sim-calibration.jld2   → DataFrame
#
# Outputs:
#   data/results-stress.jld2    → per-(ticker × composer × replication × factor) metrics
#   data/results-stress.csv
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg            = load_config()
market_ticker  = cfg["universe"]["market_ticker"]
n_paths        = Int(cfg["simulation"]["n_paths"])
seed           = Int(cfg["simulation"]["seed"])
f              = Float64(cfg["hybrid"]["idiosyncratic_floor"])
R²_thresh      = Float64(cfg["hybrid"]["r2_preserve_threshold"])
factors        = stress_factors(cfg)

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading universe / marginals / calibration..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
md = load(joinpath(_PATH_TO_DATA, "marginals.jld2"))
cd = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
marginals = md["marginals"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m_real   = G[:, market_idx]
T_eff      = length(G_m_real)
σ_m_real   = std(G_m_real)
@info "Stress setup" market = market_ticker T_eff = T_eff σ_m_real = σ_m_real factors = factors

# ── 2. Sweep composers and stress factors ──────────────────────────────────
Random.seed!(seed)

rows = NamedTuple[]

function record!(rows, ticker, composer, rep, factor, β, β_eff, flag, metrics)
    push!(rows, merge(
        (ticker = ticker, composer = composer, rep = rep,
         factor = factor, beta = β, beta_eff = β_eff, flag = flag),
        metrics
    ))
end

for factor in factors
    G_m  = scale_market(G_m_real, factor)
    σ²_m = var(G_m)
    @info "Stress run" factor = factor σ_m = sqrt(σ²_m)

    for (i, row) in enumerate(eachrow(calib))
        ticker = row.ticker
        α, β   = row.alpha, row.beta
        R²     = row.r2_real
        σ_εr   = row.sigma_eps_real

        asset_idx = findfirst(==(ticker), tickers)
        G_real    = G[:, asset_idx]

        model = marginals[ticker]
        sim_result = simulate(model, T_eff;
                              n_paths = n_paths,
                              seed = seed + i + Int(round(factor * 10_000)))

        if i % 50 == 0 || i == 1
            @info "  factor $(factor) ticker $i / $(nrow(calib)): $ticker (β=$(round(β, digits=3)))"
        end

        for r in 1:n_paths
            ε̃    = Float64.(sim_result.paths[r].observations)
            σ²_g = var(ε̃)

            g_n = compose_naive(α, β, G_m, ε̃)
            record!(rows, ticker, "naive", r, factor, β, β, string(NAIVE),
                    score_asset(g_n, G_real, G_m_real))

            g_g = compose_gaussian_sim(α, β, σ_εr, G_m, Random.default_rng())
            record!(rows, ticker, "gaussian", r, factor, β, β, string(GAUSSIAN_SIM),
                    score_asset(g_g, G_real, G_m_real))

            g_h, β_eff, flag = compose_hybrid(α, β, R², G_m, ε̃, σ²_m, σ²_g;
                                               f = f, R²_threshold = R²_thresh)
            record!(rows, ticker, "hybrid", r, factor, β, β_eff, string(flag),
                    score_asset(g_h, G_real, G_m_real))
        end
    end
end

results = DataFrame(rows)
@info "Stress result rows" total = nrow(results)

# Quick clipping summary for log visibility
clipped = filter(r -> r.composer == "hybrid" && r.flag == string(HYBRID_CLIPPED), results)
@info "Clipping incidence" n_clipped_rows = nrow(clipped) per_factor = combine(groupby(clipped, :factor), nrow => :n)

# ── 3. Persist ──────────────────────────────────────────────────────────────
results_path = joinpath(_PATH_TO_DATA, "results-stress.jld2")
jldsave(results_path; results = results, config = cfg)
csv_path = joinpath(_PATH_TO_DATA, "results-stress.csv")
CSV.write(csv_path, results)
@info "Stress results cached" jld = results_path csv = csv_path
