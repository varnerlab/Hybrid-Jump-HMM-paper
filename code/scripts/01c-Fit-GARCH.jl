# =============================================================================
# 01c-Fit-GARCH.jl
#
# Fits a GARCH(1,1) model with Student-t standardized innovations to each
# non-market ticker's OLS residual series and caches the fitted models.
# Simulation of synthetic innovation paths is done at scoring time in
# Pipeline.jl, so the seed used at scoring time controls the GARCH path
# realizations. (Earlier revisions of this script also pre-simulated and
# cached n_paths paths per ticker; that made the GARCH composer
# seed-invariant in seed-sweep runs and was removed.)
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
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

# Pkg.add("ARCHModels")  # uncomment on first run if needed
using ARCHModels

cfg           = load_config()
market_ticker = cfg["universe"]["market_ticker"]
seed          = Int(cfg["simulation"]["seed"])

models_cache = joinpath(_PATH_TO_DATA, "garch-t-models.jld2")
if isfile(models_cache)
    @info "GARCH-t models already cached — skipping."
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
@info "GARCH-t setup" n_tickers = nrow(calib) T = T

# ── 2. Fit GARCH(1,1)-t per ticker ──────────────────────────────────────────
Random.seed!(seed)

garch_models = Dict{String,Any}()
skipped      = NamedTuple[]

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
    # `fit` and `simulate` are qualified with `ARCHModels.` because
    # JumpHMM (imported in Include.jl) also exports both names; without
    # qualification Julia dispatches to JumpHMM's `fit`, which has no
    # method for GARCH types. Fits occasionally converge to a
    # nonstationary region (α+β ≥ 1); ARCHModels refuses to simulate
    # from those models, so we skip the ticker and record the reason.
    try
        model = ARCHModels.fit(GARCH{1,1}, e;
                               dist = StdT,
                               meanspec = NoIntercept{Float64})
        garch_models[ticker] = model
    catch err
        reason = sprint(showerror, err)
        @warn "GARCH skip" ticker=ticker reason=first(reason, 80)
        push!(skipped, (ticker = ticker, reason = reason))
    end
end

if !isempty(skipped)
    skipped_df = DataFrame(skipped)
    CSV.write(joinpath(_PATH_TO_DATA, "garch-t-skipped.csv"), skipped_df)
    @info "GARCH skipped tickers" n_skipped = nrow(skipped_df) n_fit = length(garch_models)
end

@info "Caching GARCH-t models to $models_cache"
jldsave(models_cache; models = garch_models)

@info "Done — $(length(garch_models)) GARCH-t fits cached."
