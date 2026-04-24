# =============================================================================
# Composers.jl
# Three composers for SIM construction over a per-asset generator draw ε̃.
# All three accept the same ε̃ vector so that the comparison across composers
# is paired (same idiosyncratic shock, different composition rule).
# =============================================================================

"""
    ConstructionFlag

Per-ticker tag identifying which branch of the hybrid construction was used.

* `HYBRID`           — variance-corrected, no clipping (`ρ ≤ 1 - f`)
* `HYBRID_CLIPPED`   — variance-corrected with `β` clipped to the floor
* `R2_PRESERVE`      — `R²`-preserving branch (tracker assets)
* `NAIVE`            — naive composition baseline
* `GAUSSIAN_SIM`     — Gaussian-residual SIM baseline
* `RESIDUAL_JUMPHMM` — JumpHMM fit to OLS residuals (not full returns), composed with `s=1`
* `BLOCK_BOOTSTRAP`  — Politis–Romano stationary bootstrap of real OLS residuals
* `GARCH_T`          — per-ticker GARCH(1,1)-t fit to OLS residuals
"""
@enum ConstructionFlag HYBRID HYBRID_CLIPPED R2_PRESERVE NAIVE GAUSSIAN_SIM RESIDUAL_JUMPHMM BLOCK_BOOTSTRAP GARCH_T

"""
    compose_naive(α, β, gm, ε̃) → g

Naive composition: `g = α + β * gm + ε̃`. No variance correction; the marginal
variance of `g` overshoots the generator variance by `β² * Var[gm]`.
"""
function compose_naive(α::Float64, β::Float64,
                       gm::AbstractVector{<:Real}, ε̃::AbstractVector{<:Real})
    return α .+ β .* gm .+ ε̃
end

"""
    compose_gaussian_sim(α, β, σ_ε_real, gm, rng) → g

Gaussian-residual SIM baseline. Replaces `ε̃` with a fresh
`N(0, σ_ε_real²)` draw, where `σ_ε_real` is the OLS residual standard
deviation from the real-data calibration. Recovers `(α, β, R²)` exactly but
loses any heavy-tail or regime structure carried by the generator.
"""
function compose_gaussian_sim(α::Float64, β::Float64, σ_ε_real::Float64,
                              gm::AbstractVector{<:Real}, rng::AbstractRNG)
    T = length(gm)
    return α .+ β .* gm .+ σ_ε_real .* randn(rng, T)
end

"""
    compose_hybrid(α, β, R²_real, gm, ε̃, σ²_m, σ²_gen;
                   f=0.10, R²_threshold=0.80) → (g, β_eff, flag)

Hybrid SIM composition. Branch selection is driven by `R²_real`:

* If `R²_real ≥ R²_threshold`, use the `R²`-preserving branch
  (Eq. 8 of the paper): set `σ²_ε_target = β² σ²_m (1 - R²_real)/R²_real` and
  rescale `ε̃` to that target variance.
* Otherwise, use the variance-preserving branch (Eqs. 4–7): if
  `ρ = β² σ²_m / σ²_gen ≤ 1 - f`, set `s² = 1 - ρ` and `β_eff = β`;
  if `ρ > 1 - f`, clip `β` to
  `sign(β) sqrt((1-f) σ²_gen / σ²_m)` and set `s² = f`.

Returns the composed series `g`, the effective `β_eff`, and a
`ConstructionFlag` tagging the branch.
"""
function compose_hybrid(α::Float64, β::Float64, R²_real::Float64,
                        gm::AbstractVector{<:Real}, ε̃::AbstractVector{<:Real},
                        σ²_m::Float64, σ²_gen::Float64;
                        f::Float64 = 0.10,
                        R²_threshold::Float64 = 0.80)

    @assert 0.0 < f < 1.0           "idiosyncratic floor must lie in (0, 1)"
    @assert 0.0 < R²_threshold ≤ 1.0 "R² threshold must lie in (0, 1]"
    @assert σ²_m > 0.0              "market variance must be positive"

    if R²_real ≥ R²_threshold
        # R²-preserving branch
        σ²_ε_target = (R²_real ≥ 1.0 - 1e-12) ? 0.0 :
                      β^2 * σ²_m * (1.0 - R²_real) / R²_real
        scale = (σ²_gen > 0.0 && σ²_ε_target > 0.0) ?
                sqrt(σ²_ε_target / σ²_gen) : 0.0
        ε     = scale .* ε̃
        β_eff = β
        flag  = R2_PRESERVE
    else
        # variance-preserving branch with optional β clipping
        ρ = β^2 * σ²_m / max(σ²_gen, 1e-30)
        if ρ > 1.0 - f
            β_eff = sign(β) * sqrt((1.0 - f) * σ²_gen / σ²_m)
            s²    = f
            flag  = HYBRID_CLIPPED
        else
            β_eff = β
            s²    = 1.0 - ρ
            flag  = HYBRID
        end
        ε = sqrt(s²) .* ε̃
    end

    g = α .+ β_eff .* gm .+ ε
    return g, β_eff, flag
end

"""
    compose_residual_jumphmm(α, β, gm, ε̃_resid) → g

Alternative composition baseline: the per-asset generator is fit to the OLS
residual series `e_i(t) = r_i(t) - α_i - β_i r_m(t)` rather than to the full
return series. A draw `ε̃_resid` from that residual-fit generator is already
idiosyncratic by construction, so composition proceeds without variance
correction: `g = α + β · gm + ε̃_resid`.

Functionally identical to `compose_naive`; the distinction is upstream (where
the generator was fit). We keep it as a separate entry point so the caller's
intent and the construction flag are unambiguous in the output scoreboard.
"""
function compose_residual_jumphmm(α::Float64, β::Float64,
                                   gm::AbstractVector{<:Real},
                                   ε̃_resid::AbstractVector{<:Real})
    return α .+ β .* gm .+ ε̃_resid
end

"""
    compose_block_bootstrap(α, β, gm, residuals_real, block_length, rng) → g

Politis–Romano stationary bootstrap composition baseline. Draws a bootstrap
replicate of length `length(gm)` from the real OLS residual series
`residuals_real`, using random block lengths with mean `block_length`, then
composes `g = α + β · gm + ε̂` with no variance correction. Preserves any
serial correlation present in the real residuals up to the block scale.
"""
function compose_block_bootstrap(α::Float64, β::Float64,
                                  gm::AbstractVector{<:Real},
                                  residuals_real::AbstractVector{<:Real},
                                  block_length::Real,
                                  rng::AbstractRNG)
    T = length(gm)
    n = length(residuals_real)
    @assert n > 0 "residual series must be non-empty"
    @assert block_length > 0 "mean block length must be positive"
    p = 1.0 / block_length
    ε̂ = Vector{Float64}(undef, T)
    idx = rand(rng, 1:n)
    for t in 1:T
        if t > 1 && rand(rng) < p
            idx = rand(rng, 1:n)
        end
        ε̂[t] = residuals_real[idx]
        idx = idx == n ? 1 : idx + 1
    end
    return α .+ β .* gm .+ ε̂
end

"""
    compose_garch_t(α, β, gm, ε̃_garch) → g

GARCH(1,1)-t composition baseline. The innovation `ε̃_garch` is assumed to be
a simulated path from a per-ticker GARCH(1,1) with Student-t standardized
innovations, fit to the real OLS residual series. Composes with no variance
correction: `g = α + β · gm + ε̃_garch`.

Fitting is performed externally (see `scripts/01c-Fit-GARCH.jl`) against
`ARCHModels.jl`; this function only consumes the sampled innovations so that
`Composers.jl` remains free of third-party dependencies.
"""
function compose_garch_t(α::Float64, β::Float64,
                          gm::AbstractVector{<:Real},
                          ε̃_garch::AbstractVector{<:Real})
    return α .+ β .* gm .+ ε̃_garch
end

"""
    apply_copula_reorder!(ε, U)

Rank-reorder a vector of innovations `ε` so its empirical rank order matches
the ranks of a uniform draw `U` of the same length. Used to inject
cross-sectional dependence into the residuals before composition. Modifies
nothing in place; returns the reordered vector.
"""
function apply_copula_reorder(ε::AbstractVector{<:Real}, U::AbstractVector{<:Real})
    @assert length(ε) == length(U) "ε and U must have the same length"
    sorted_eps    = sort(collect(ε))
    copula_ranks  = ordinalrank(U)
    return sorted_eps[copula_ranks]
end
