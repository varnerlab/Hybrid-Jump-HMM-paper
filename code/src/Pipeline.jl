# =============================================================================
# Pipeline.jl
# Universe loading, fitting, and artifact I/O helpers. Delegates to
# VLQuantitativeFinancePackage for data and JumpHMM for model fitting.
# =============================================================================

"""
    load_config(path = _PATH_TO_CONFIG) → Dict

Read the experiment configuration TOML.
"""
function load_config(path::AbstractString = _PATH_TO_CONFIG)
    isfile(path) || throw(ArgumentError("config file not found: $path"))
    return TOML.parsefile(path)
end

"""
    load_universe(min_obs) → (tickers, prices)

Load the JDIQ training universe via VLQuantitativeFinancePackage and filter to
tickers whose history matches the maximum trading-day count in the dataset
(AAPL is used as the reference). The `min_obs` argument is enforced as a
lower bound on the resulting common length: if AAPL has fewer than `min_obs`
days, the function errors. Returns a sorted vector of ticker symbols and a
`(T × N)` price matrix in the same order.
"""
function load_universe(min_obs::Int)
    raw = MyTrainingMarketDataSet()["dataset"]
    max_days = raw["AAPL"] |> nrow
    max_days >= min_obs ||
        error("AAPL has $max_days days; expected at least $min_obs")
    @info "Universe scan: AAPL has $max_days trading days"

    keep = Dict{String,DataFrame}()
    for (ticker, df) in raw
        if nrow(df) == max_days
            keep[ticker] = df
        end
    end

    tickers = sort(collect(keys(keep)))
    prices  = Matrix{Float64}(undef, max_days, length(tickers))
    for (j, t) in enumerate(tickers)
        prices[:, j] = keep[t].close
    end

    @info "Universe loaded" n_assets = length(tickers) n_obs = max_days
    return tickers, prices
end

"""
    growth_rate_matrix(prices; rf, dt) → G

Annualized excess log growth rates from a `(T × N)` price matrix, matching
the JumpHMM convention `G_t = (1/dt) log(P_t / P_{t-1}) - rf`. Delegates to
`JumpHMM.excess_growth_rates` so the units are identical to the per-asset
HMM marginals fit downstream.
"""
function growth_rate_matrix(prices::AbstractMatrix{<:Real};
                            rf::Float64 = 0.0,
                            dt::Float64 = 1.0 / 252.0)
    return JumpHMM.excess_growth_rates(prices; rf = rf, dt = dt)
end

"""
    fit_per_ticker_marginals(prices, tickers; cfg) → Dict{String,JumpHiddenMarkovModel}

Fit a JumpHMM marginal for each ticker. Caches the result to
`data/marginals.jld2` so subsequent calls hit the cache.
"""
function fit_per_ticker_marginals(prices::AbstractMatrix{<:Real},
                                  tickers::Vector{String};
                                  cfg::Dict)
    cache = joinpath(_PATH_TO_DATA, "marginals.jld2")
    if isfile(cache)
        @info "Loading cached marginals from $cache"
        return load(cache)["marginals"]
    end

    rf = Float64(cfg["hmm"]["risk_free_rate"])
    N  = Int(cfg["hmm"]["N"])
    ν  = Float64(cfg["hmm"]["nu"])
    dt = Float64(cfg["hmm"]["dt"])

    marginals = Dict{String,JumpHiddenMarkovModel}()
    for (j, t) in enumerate(tickers)
        @info "Fitting marginal $j / $(length(tickers)): $t"
        marginals[t] = fit(JumpHiddenMarkovModel, prices[:, j];
                           rf = rf, N = N, ν = ν, dt = dt)
    end

    @info "Caching marginals to $cache"
    jldsave(cache; marginals = marginals)
    return marginals
end

"""
    calibrate_sim(G, tickers, market_ticker) → DataFrame

OLS regression per ticker of `G[:, i]` on `G[:, market_idx]`. Returns a
DataFrame with columns `ticker, alpha, beta, r2_real, sigma_eps_real,
sigma_gen` (the last is filled in later from the HMM marginal).
"""
function calibrate_sim(G::AbstractMatrix{<:Real},
                       tickers::Vector{String},
                       market_ticker::String)
    market_idx = findfirst(==(market_ticker), tickers)
    market_idx === nothing && throw(ArgumentError("market ticker not in universe"))
    G_m   = G[:, market_idx]
    G_m̄   = mean(G_m)
    G_m_v = var(G_m)

    n = length(tickers)
    df = DataFrame(
        ticker         = String[],
        alpha          = Float64[],
        beta           = Float64[],
        r2_real        = Float64[],
        sigma_eps_real = Float64[],
    )

    for i in 1:n
        t = tickers[i]
        if t == market_ticker
            continue
        end
        G_i  = G[:, i]
        G_ī  = mean(G_i)
        β    = cov(G_i, G_m) / G_m_v
        α    = G_ī - β * G_m̄
        resid = G_i .- α .- β .* G_m
        σ_ε  = std(resid)
        SS_r = dot(resid, resid)
        SS_t = sum(abs2, G_i .- G_ī)
        r²   = SS_t > 0.0 ? 1.0 - SS_r / SS_t : 0.0
        push!(df, (t, α, β, r², σ_ε))
    end
    return df
end
