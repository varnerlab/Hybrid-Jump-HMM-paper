# =============================================================================
# Metrics.jl
# Recovery metrics (β, α, R²), marginal-fidelity metrics (KS, AD, W1, Hill,
# kurtosis), and cross-sectional metrics (rank-correlation Frobenius distance)
# for evaluating composers against real per-asset growth-rate histories.
# =============================================================================

# --- Recovery metrics --------------------------------------------------------

"""
    sim_recovery(g, gm) → (α̂, β̂, R²̂)

OLS regression of synthetic `g` on market path `gm`, returning the recovered
intercept, slope, and coefficient of determination.
"""
function sim_recovery(g::AbstractVector{<:Real}, gm::AbstractVector{<:Real})
    @assert length(g) == length(gm) "g and gm must have the same length"
    g_bar  = mean(g)
    m_bar  = mean(gm)
    m_var  = var(gm)
    m_var > 0.0 || throw(ArgumentError("market path has zero variance"))
    β̂      = cov(g, gm) / m_var
    α̂      = g_bar - β̂ * m_bar
    resid  = g .- α̂ .- β̂ .* gm
    SS_res = dot(resid, resid)
    SS_tot = sum(abs2, g .- g_bar)
    R²̂     = SS_tot > 0.0 ? 1.0 - SS_res / SS_tot : 0.0
    return α̂, β̂, R²̂
end

# --- Marginal-fidelity metrics ----------------------------------------------

"""
    ks_pvalue(g_synth, g_real) → p

Two-sample Kolmogorov–Smirnov p-value comparing synthetic and real growth-rate
samples.
"""
function ks_pvalue(g_synth::AbstractVector{<:Real}, g_real::AbstractVector{<:Real})
    return pvalue(ApproximateTwoSampleKSTest(collect(g_synth), collect(g_real)))
end

"""
    ad_pvalue(g_synth, g_real) → p

Two-sample Anderson–Darling p-value.
"""
function ad_pvalue(g_synth::AbstractVector{<:Real}, g_real::AbstractVector{<:Real})
    return pvalue(KSampleADTest(collect(g_synth), collect(g_real)))
end

"""
    wasserstein1(g_synth, g_real) → W1

Empirical 1-Wasserstein distance between two univariate samples (sort, then
mean absolute difference of order statistics).
"""
function wasserstein1(g_synth::AbstractVector{<:Real}, g_real::AbstractVector{<:Real})
    n = length(g_synth)
    m = length(g_real)
    s = sort(g_synth)
    r = sort(g_real)
    if n == m
        return mean(abs.(s .- r))
    end
    # piecewise-constant CDF integration via interpolated quantiles
    grid = range(0.0, 1.0; length = max(n, m))
    qs   = [quantile(s, q) for q in grid]
    qr   = [quantile(r, q) for q in grid]
    return mean(abs.(qs .- qr))
end

"""
    hill_index(x; tail_frac=0.05) → ξ

Hill estimator of the tail index for the upper tail. `tail_frac` selects the
fraction of largest order statistics used (0.05 = top 5%).
"""
function hill_index(x::AbstractVector{<:Real}; tail_frac::Float64 = 0.05)
    @assert 0.0 < tail_frac < 1.0 "tail_frac must lie in (0, 1)"
    sx = sort(x; rev = true)
    k  = max(2, floor(Int, tail_frac * length(sx)))
    @assert sx[k] > 0.0 "tail order statistic must be positive for Hill estimator"
    s  = 0.0
    for i in 1:(k - 1)
        s += log(sx[i] / sx[k])
    end
    return s / (k - 1)
end

"""
    excess_kurtosis(x) → κ - 3
"""
excess_kurtosis(x::AbstractVector{<:Real}) = kurtosis(x)

# --- Cross-sectional metrics -------------------------------------------------

"""
    rank_corr_frobenius(G_synth, G_real) → d

Frobenius distance between the Spearman rank-correlation matrices of the
synthetic and real `(T × N)` growth-rate matrices. Diagonals are zeroed before
the norm so the metric reflects off-diagonal dependence only.
"""
function rank_corr_frobenius(G_synth::AbstractMatrix{<:Real},
                             G_real::AbstractMatrix{<:Real})
    @assert size(G_synth, 2) == size(G_real, 2) "column counts must match"
    Cs = corspearman(G_synth)
    Cr = corspearman(G_real)
    D  = Cs .- Cr
    for k in 1:size(D, 1)
        D[k, k] = 0.0
    end
    return norm(D)
end

# --- Aggregate scorer --------------------------------------------------------

"""
    score_asset(g_synth, g_real, gm) → NamedTuple

Run all per-asset metrics on a single composed series against the real growth
rates and the market path. Returns a NamedTuple keyed by metric name, ready to
be pushed into a DataFrame row.
"""
function score_asset(g_synth::AbstractVector{<:Real},
                     g_real::AbstractVector{<:Real},
                     gm::AbstractVector{<:Real})
    α̂, β̂, R²̂ = sim_recovery(g_synth, gm)
    return (
        α_hat   = α̂,
        β_hat   = β̂,
        R²_hat  = R²̂,
        ks_p    = ks_pvalue(g_synth, g_real),
        ad_p    = ad_pvalue(g_synth, g_real),
        w1      = wasserstein1(g_synth, g_real),
        hill_up = hill_index(abs.(g_synth)),
        kurt    = excess_kurtosis(g_synth),
        var_g   = var(g_synth),
    )
end
