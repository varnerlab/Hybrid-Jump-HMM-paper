# =============================================================================
# Composers.jl
# Three composers for SIM construction over a per-asset generator draw ε̃.
# All three accept the same ε̃ vector so that the comparison across composers
# is paired (same idiosyncratic shock, different composition rule).
# =============================================================================

"""
    ConstructionFlag

Per-ticker tag identifying which branch of the hybrid construction was used.

* `HYBRID`         — variance-corrected, no clipping (`ρ ≤ 1 - f`)
* `HYBRID_CLIPPED` — variance-corrected with `β` clipped to the floor
* `R2_PRESERVE`    — `R²`-preserving branch (tracker assets)
* `NAIVE`          — naive composition baseline
* `GAUSSIAN_SIM`   — Gaussian-residual SIM baseline
"""
@enum ConstructionFlag HYBRID HYBRID_CLIPPED R2_PRESERVE NAIVE GAUSSIAN_SIM

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
