# =============================================================================
# 08-Synthetic-Tracker-Eval.jl
#
# R²-preserve branch correctness check on synthetic tracker assets.
#
# Motivation: in the 423-ticker S&P 500 universe only QQQ (R²=0.861) and SPYG
# (R²=0.933) exceed the R²_preserve = 0.80 threshold; the highest individual
# stock (BLK) sits at R²=0.645. The empirical gap between single stocks and
# index trackers is structural, so the r2-preserve branch has n=2 empirical
# support. To validate the branch under controlled conditions, we construct
# synthetic tracker assets with a calibrated R² and β, run them through
# compose_hybrid, and check that β̂ and R² recover to their target values
# across the full intended R² range.
#
# Construction: for a target (β_true, R²_true),
#   σ_ε_target² = β_true² · σ_m² · (1 - R²_true) / R²_true
#   g_i(t)      = α_true + β_true · g_m(t) + σ_ε_target · ξ(t),  ξ ~ N(0, 1)
# This matches Eq. (9) of the paper exactly; the JumpHMM fit on this series
# should recover the same relationship.
#
# Outputs:
#   data/synth-tracker.csv
#   data/synth-tracker.jld2
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg           = load_config()
market_ticker = cfg["universe"]["market_ticker"]
seed          = Int(cfg["simulation"]["seed"])
f             = Float64(cfg["hybrid"]["idiosyncratic_floor"])
R²_thresh     = Float64(cfg["hybrid"]["r2_preserve_threshold"])
n_paths       = Int(cfg["simulation"]["n_paths"])

# ── 1. Load the real SPY market path as the reference g_m ───────────────────
@info "Loading market path for synthetic tracker construction..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
tickers = ud["tickers"]
G       = ud["growth_rates"]
market_idx = findfirst(==(market_ticker), tickers)
G_m  = G[:, market_idx]
σ²_m = var(G_m)
μ_m  = mean(G_m)
T    = length(G_m)

# ── 2. Target R² grid and matching β values ─────────────────────────────────
R²_grid = [0.80, 0.85, 0.90, 0.95, 0.99]
β_grid  = [0.80, 1.00, 1.20]
α_true  = 0.0

rng = MersenneTwister(seed)
rows = NamedTuple[]

function record!(rows, β_true, R²_true, rep, β̂, R²̂, ks_p, kurt, flag_out, β_eff_out, var_g)
    push!(rows, (beta_true = β_true, r2_true = R²_true, rep = rep,
                 beta_hat = β̂, r2_hat = R²̂, ks_p = ks_p, kurt = kurt,
                 flag = flag_out, beta_eff = β_eff_out, var_g = var_g))
end

for β_true in β_grid, R²_true in R²_grid
    σ²_ε_target = β_true^2 * σ²_m * (1.0 - R²_true) / R²_true
    σ_ε_target  = sqrt(σ²_ε_target)
    @info "Synthetic tracker" β_true=β_true R²_true=R²_true σ_ε_target=round(σ_ε_target, digits=3)

    for rep in 1:n_paths
        # build a synthetic tracker path with the exact target R²
        ξ = randn(rng, T)
        g_real = α_true .+ β_true .* G_m .+ σ_ε_target .* ξ

        # synthetic generator draw "ε̃" = fresh Gaussian innovation matched in
        # variance to the generator-target convention (σ²_gen = σ²_ε_target since
        # the "generator marginal" here IS the residual distribution)
        σ²_gen = σ²_ε_target
        ε̃     = σ_ε_target .* randn(rng, T)

        # run the hybrid composer; branch should be R2_PRESERVE because
        # R²_true ≥ R²_thresh = 0.80
        g_h, β_eff, flag = compose_hybrid(α_true, β_true, R²_true, G_m, ε̃,
                                          σ²_m, σ²_gen;
                                          f = f, R²_threshold = R²_thresh)
        α̂, β̂, R²̂ = sim_recovery(g_h, G_m)
        ks_p = ks_pvalue(g_h, g_real)

        record!(rows, β_true, R²_true, rep, β̂, R²̂, ks_p, excess_kurtosis(g_h),
                string(flag), β_eff, var(g_h))
    end
end

df = DataFrame(rows)
@info "Synthetic tracker result rows" total = nrow(df)

# Aggregate summary per (β_true, R²_true)
summary = combine(groupby(df, [:beta_true, :r2_true]),
    :beta_hat => (x -> median(x)) => :median_beta_hat,
    :beta_hat => (x -> std(x))    => :sd_beta_hat,
    :r2_hat   => (x -> median(x)) => :median_r2_hat,
    :r2_hat   => (x -> std(x))    => :sd_r2_hat,
    :ks_p     => (p -> 100 * mean(p .> 0.05)) => :ks_pass_pct,
    :flag     => (f -> join(unique(f), ",")) => :flags,
)
@info "Synthetic tracker summary"
show(summary, allcols = true, allrows = true); println()

CSV.write(joinpath(_PATH_TO_DATA, "synth-tracker.csv"), df)
CSV.write(joinpath(_PATH_TO_DATA, "synth-tracker-summary.csv"), summary)
jldsave(joinpath(_PATH_TO_DATA, "synth-tracker.jld2");
        results = df, summary = summary, config = cfg)
@info "Synthetic tracker results persisted."
