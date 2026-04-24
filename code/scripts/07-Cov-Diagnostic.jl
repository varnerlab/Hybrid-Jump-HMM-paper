# =============================================================================
# 07-Cov-Diagnostic.jl
#
# Empirical cross-term diagnostic: sample per-ticker innovations from each
# generator and report Cov(ε̃, g_m) / (σ_gen · σ_m) to verify that the
# independence assumption used in the variance-correction derivation
# (method.tex:57-63) holds in practice. Peer-review P2 (R2).
#
# Produces a histogram per composer across the 423 tickers, written to
# figs/cov-diagnostic.pdf, plus a CSV of the raw (ticker × composer) values.
#
# Inputs:
#   data/universe.jld2              → tickers, growth_rates
#   data/marginals.jld2             → full-return JumpHMM marginals
#   data/marginals-residuals.jld2   → (optional) residual-fit marginals
#   data/sim-calibration.jld2       → (α, β, R², σ_eps_real)
#
# Outputs:
#   data/cov-diagnostic.csv
#   figs/cov-diagnostic.pdf
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg           = load_config()
market_ticker = cfg["universe"]["market_ticker"]
n_paths       = Int(cfg["simulation"]["n_paths"])
seed          = Int(cfg["simulation"]["seed"])

@info "Loading universe / marginals / calibration for Cov diagnostic..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
md = load(joinpath(_PATH_TO_DATA, "marginals.jld2"))
cd = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
marginals = md["marginals"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m        = G[:, market_idx]
σ_m        = std(G_m)
T          = length(G_m)

residual_cache = joinpath(_PATH_TO_DATA, "marginals-residuals.jld2")
marginals_resid = isfile(residual_cache) ?
    load(residual_cache)["marginals"] : nothing

rows = NamedTuple[]

function cov_norm(ε, g_m)
    σ_ε = std(ε)
    σ_m_local = std(g_m)
    return (σ_ε > 0.0 && σ_m_local > 0.0) ?
        cov(ε, g_m) / (σ_ε * σ_m_local) : 0.0
end

for (i, row) in enumerate(eachrow(calib))
    ticker = row.ticker
    α, β   = row.alpha, row.beta

    # Full-return JumpHMM draw
    sim_full = simulate(marginals[ticker], T; n_paths = n_paths, seed = seed + i)
    cs_full = [cov_norm(Float64.(sim_full.paths[r].observations), G_m) for r in 1:n_paths]
    push!(rows, (ticker = ticker, composer = "jumphmm_full",
                 cov_norm_mean = mean(cs_full), cov_norm_sd = std(cs_full)))

    # Residual-fit JumpHMM draw (if available)
    if marginals_resid !== nothing && haskey(marginals_resid, ticker)
        sim_r = simulate(marginals_resid[ticker], T;
                         n_paths = n_paths, seed = seed + i + 1_000_000)
        cs_r = [cov_norm(Float64.(sim_r.paths[r].observations), G_m) for r in 1:n_paths]
        push!(rows, (ticker = ticker, composer = "jumphmm_residual",
                     cov_norm_mean = mean(cs_r), cov_norm_sd = std(cs_r)))
    end

    # Gaussian draw (control — independent by construction)
    σ_εr = row.sigma_eps_real
    rng = MersenneTwister(seed + i + 2_000_000)
    cs_g = [cov_norm(σ_εr .* randn(rng, T), G_m) for r in 1:n_paths]
    push!(rows, (ticker = ticker, composer = "gaussian",
                 cov_norm_mean = mean(cs_g), cov_norm_sd = std(cs_g)))

    if i % 50 == 0
        @info "  Cov diagnostic ticker $i / $(nrow(calib))"
    end
end

diag = DataFrame(rows)
csv_path = joinpath(_PATH_TO_DATA, "cov-diagnostic.csv")
CSV.write(csv_path, diag)
@info "Cov diagnostic CSV → $csv_path"

# ── Plot ────────────────────────────────────────────────────────────────────
p = plot(; layout = (1, 1), size = (720, 420),
        xlabel = "Cov(ε̃, g_m) / (σ_ε σ_m)", ylabel = "count",
        title = "Empirical cross-term per ticker (mean over $n_paths draws)")
for (name, color) in (("jumphmm_full", :steelblue),
                       ("jumphmm_residual", :orange),
                       ("gaussian", :gray))
    sub = filter(r -> r.composer == name, diag)
    if nrow(sub) > 0
        histogram!(p, sub.cov_norm_mean; bins = 40, alpha = 0.55,
                   label = name, color = color, normalize = false)
    end
end
fig_path = joinpath(_PATH_TO_FIGS, "cov-diagnostic.pdf")
savefig(p, fig_path)
@info "Cov diagnostic figure → $fig_path"
