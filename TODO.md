# TODO — Peer-review revision run sequence

State as of commit `b8b701f`. Scripts are in `code/scripts/`; all take their
config from `code/config.toml`. Run from the repo root.

## One-time setup

- [ ] Install GARCH dependency:
      `julia --project=code -e 'using Pkg; Pkg.add("ARCHModels")'`
      (Required only for `01c-Fit-GARCH.jl`; skip if you drop the GARCH
      baseline.)

## Required — produces the new headline numbers

- [ ] `julia --project=code code/scripts/01b-Fit-Residual-Marginals.jl`
      Hours. Fits 423 JumpHMM models on OLS residuals. Writes
      `data/marginals-residuals.jld2`.
- [ ] `julia --project=code code/scripts/01c-Fit-GARCH.jl`
      ~30 min. Needs ARCHModels.jl. Writes `data/garch-t-models.jld2`
      and `data/garch-t-sims.jld2`.
- [ ] `julia --project=code code/scripts/03-Compose-And-Evaluate.jl`
      ~30 min. Rewrites `data/results.jld2` with all six composers.
- [ ] `julia --project=code code/scripts/04-Tables.jl`
      Fast. Rebuilds Tables 1–3 with the new composer rows.
- [ ] `julia --project=code code/scripts/05-Figures.jl`
      Fast. Existing Figs 1–3 (unchanged composer set).

This minimum is enough to rewrite Abstract + Table 1 for the revision.

## Optional — each answers a specific reviewer ask

- [ ] `julia --project=code code/scripts/03b-Stress-Eval.jl`
      ~40 min. Clipping branch validation at 2× and 3× SPY variance.
      Writes `data/results-stress.jld2`. Addresses R1/R2 P1 item 4.
- [ ] `julia --project=code code/scripts/06-VaR-Backtest.jl`
      ~10 min. Per-ticker 95% / 99% VaR + Kupiec coverage test across
      all composers. Writes `data/var-backtest.jld2`. Addresses R3's
      downstream-validation demand.
- [ ] `julia --project=code code/scripts/07-Cov-Diagnostic.jl`
      ~5 min. Independence-assumption diagnostic: histogram of
      `Cov(ε̃, g_m) / (σ_ε σ_m)` per composer. Writes
      `data/cov-diagnostic.csv` and `figs/cov-diagnostic.pdf`.
      Addresses R2 P2.
- [ ] `julia --project=code code/scripts/08-Synthetic-Tracker-Eval.jl`
      ~5 min. `r2-preserve` branch correctness at
      R² ∈ {0.80, 0.85, 0.90, 0.95, 0.99} × β ∈ {0.80, 1.00, 1.20}.
      Writes `data/synth-tracker.jld2`. Supports the "n=2 is structural"
      reframe.
- [ ] `julia --project=code code/scripts/09-Extra-Figures.jl`
      Fast. Produces Figs 4–10, each conditional on its inputs existing.
      Run this last so it picks up whatever has been generated.

## Overnight-grade — not required for arXiv; useful for Quantitative Finance submission

- [ ] `julia --project=code code/scripts/03c-Seed-Sweep.jl`
      ~4 h. 8 seeds × 6 composers; produces uncertainty CIs.
      Writes `data/results-seed-*.jld2`. Addresses R2 P1 item 3.
- [ ] `julia --project=code code/scripts/03d-Sensitivity-Sweep.jl`
      ~1.5 h. Hybrid-only, 5 × 5 grid over
      (R²_threshold ∈ {0.70, 0.75, 0.80, 0.85, 0.90},
       f ∈ {0.05, 0.10, 0.15, 0.20, 0.30}). Writes
      `data/results-thresh-*-f-*.jld2`. Addresses R2 P2.

## Minimum-viable path to a revised arXiv draft

```
01b  →  01c  →  03  →  06  →  09
```

Table 1 with 6 composers + the downstream VaR validation + the new figure
set. Clipping stress (03b), seed sweep (03c), sensitivity (03d) can be
deferred to the QF resubmission.

## Then

- [ ] Recompile `paper/Paper_v1.tex` (expect warnings for references to
      Tables 4–5 and Figures 4–10 until those environments are added).
- [ ] Ping Claude to run Task 12 (paper revisions: Abstract, Tables 1 and 3
      expansion, new Tables 4–5 and Figure environments in `Paper_v1.tex`,
      References bib TODO resolution).
