# =============================================================================
# 01c-Fit-GARCH.jl
#
# Fits a GARCH(1,1) model with Student-t standardized innovations to each
# non-market ticker's OLS residual series, then simulates n_paths synthetic
# innovation paths per ticker and caches them. The scoring loop in 03 then
# consumes the simulated paths directly without re-entering the ARCHModels
# simulate/fit call, keeping that dependency contained to this script.
#
# Used by the `compose_garch_t` baseline added to Composers.jl in response
# to peer-review point P1 (R3: the paper cites Bollerslev 1986 but never
# benchmarks against GARCH).
#
# Dependency:
#   This script requires ARCHModels.jl. If it is not already installed in the
#   project environment, uncomment the Pkg.add line below, or run in the REPL:
#       Pkg.activate("code"); Pkg.add("ARCHModels")
#
# Inputs:
#   data/universe.jld2          → tickers, growth_rates (T × N)
#   data/sim-calibration.jld2   → DataFrame with (α, β, R², σ_eps_real)
#
# Outputs:
#   data/garch-t-models.jld2    → Dict{String, Any} fitted GARCH-t models
#   data/garch-t-sims.jld2      → Dict{String, Vector{Vector{Float64}}}
#                                  one outer entry per ticker, inner length = n_paths
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

# Pkg.add("ARCHModels")  # uncomment on first run if needed
using ARCHModels

cfg           = load_config()
market_ticker = cfg["universe"]["market_ticker"]
n_paths       = Int(cfg["simulation"]["n_paths"])
seed          = Int(cfg["simulation"]["seed"])

models_cache = joinpath(_PATH_TO_DATA, "garch-t-models.jld2")
sims_cache   = joinpath(_PATH_TO_DATA, "garch-t-sims.jld2")
if isfile(models_cache) && isfile(sims_cache)
    @info "GARCH-t models and sims already cached — skipping."
    exit(0)
end

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading universe and calibration for GARCH-t fitting..."
ud = load(joinpath(_PATH_TO_DATA, "universe.jld2"))
cd = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))

tickers   = ud["tickers"]
G         = ud["growth_rates"]
calib     = cd["calibration"]

market_idx = findfirst(==(market_ticker), tickers)
G_m        = G[:, market_idx]
T          = length(G_m)
@info "GARCH-t setup" n_tickers = nrow(calib) T = T n_paths = n_paths

# ── 2. Fit GARCH(1,1)-t per ticker and pre-sample ───────────────────────────
Random.seed!(seed)

garch_models = Dict{String,Any}()
garch_sims   = Dict{String,Vector{Vector{Float64}}}()

for (i, row) in enumerate(eachrow(calib))
    ticker = row.ticker
    α, β   = row.alpha, row.beta
    j      = findfirst(==(ticker), tickers)

    # OLS residual series (already zero-mean by construction)
    e = G[:, j] .- α .- β .* G_m

    if i % 25 == 0 || i == 1
        @info "  fitting $i / $(nrow(calib)): $ticker"
    end

    # GARCH(1,1) with Student-t standardized innovations, no mean.
    spec  = GARCH{1,1}
    dist  = StdT
    model = fit(spec, e; dist = dist, meanspec = NoIntercept)

    # pre-sample n_paths synthetic residual series of length T
    sims_this = Vector{Vector{Float64}}(undef, n_paths)
    for r in 1:n_paths
        sims_this[r] = simulate(model, T).data
    end

    garch_models[ticker] = model
    garch_sims[ticker]   = sims_this
end

@info "Caching GARCH-t models to $models_cache"
jldsave(models_cache; models = garch_models)

@info "Caching GARCH-t sim paths to $sims_cache"
jldsave(sims_cache; sims = garch_sims)

@info "Done — $(length(garch_models)) GARCH-t fits and $(length(garch_sims)) sim sets cached."
