# =============================================================================
# SyntheticMarket.jl
# Market-path constructors for stress-test evaluation. The main experiment
# (scripts/03) uses the real SPY history as g_m; this module provides
# alternative paths that either scale the real market's volatility or generate
# a fresh synthetic path, so the hybrid composer's clipping branch can be
# exercised (peer-review P1: R1/R2 noted that clipping never triggered on the
# real SPY path and the HYBRID_CLIPPED branch was therefore unvalidated).
# =============================================================================

"""
    scale_market(gm_real, factor) → gm_scaled

Return a length-preserved market path with sample variance scaled by `factor²`
relative to `gm_real`. Preserves temporal ordering and mean-reverting
structure of the real path; only the variance is rescaled. Concretely:

    gm_scaled[t] = μ_real + factor · (gm_real[t] - μ_real)

where `μ_real = mean(gm_real)`. With `factor = 2`, clipping in
`compose_hybrid` triggers on high-β tickers whose idiosyncratic-variance
floor (`1 − f`) is exceeded under the inflated market component.
"""
function scale_market(gm_real::AbstractVector{<:Real}, factor::Real)
    @assert factor > 0.0 "scale factor must be positive"
    μ = mean(gm_real)
    return μ .+ factor .* (gm_real .- μ)
end

"""
    stress_factors(cfg = nothing) → Vector{Float64}

Return the list of market-variance scale factors for the stress-test sweep.
Defaults to `[1.0, 2.0, 3.0]`; a non-nothing `cfg` may override via
`cfg["stress"]["factors"]` in config.toml.
"""
function stress_factors(cfg = nothing)
    default = Float64[1.0, 2.0, 3.0]
    cfg === nothing && return default
    haskey(cfg, "stress") && haskey(cfg["stress"], "factors") ?
        Float64.(cfg["stress"]["factors"]) : default
end
