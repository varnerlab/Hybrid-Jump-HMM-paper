# Code: Variance-Corrected Single-Index Composition

Reproducibility code for the paper *Variance-Corrected Single-Index Composition
for Heavy-Tailed Return Generators*.

## Layout

```
code/
├── Project.toml         Julia environment
├── Include.jl           paths, package imports, source loading
├── config.toml          experiment configuration
├── src/
│   ├── Composers.jl     three composers: naive, Gaussian SIM, hybrid (paired ε̃)
│   ├── Metrics.jl       KS, AD, Wasserstein-1, Hill tail index, β recovery
│   └── Pipeline.jl      universe loading, fitting, artifact I/O
├── scripts/
│   ├── 01-Fit-Marginals.jl       per-ticker JumpHMM fits over the IS window
│   ├── 02-Calibrate-SIM.jl       per-ticker OLS against SPY
│   ├── 03-Compose-And-Evaluate.jl  paired draws through the three composers
│   ├── 04-Tables.jl              paper tables
│   └── 05-Figures.jl             paper figures
├── data/                cached JLD2 artifacts (gitignored)
└── figs/                paper figures (gitignored)
```

## Dependencies

The hybrid composer is implemented in
[`JumpHMM.jl`](https://github.com/varnerlab/JumpHMM.jl) as
`HybridSingleIndexModel`. The paper code uses `JumpHMM.jl` for per-ticker
marginal fits and re-implements the three composers locally so that all three
share the same per-ticker generator draw `ε̃`, yielding a paired comparison.

The 424-ticker universe and price history are loaded via
`MyTrainingMarketDataSet()` from
[`VLQuantitativeFinancePackage.jl`](https://github.com/varnerlab/VLQuantitativeFinancePackage.jl).

## Running

```julia
include("Include.jl")
include("scripts/01-Fit-Marginals.jl")
include("scripts/02-Calibrate-SIM.jl")
include("scripts/03-Compose-And-Evaluate.jl")
include("scripts/04-Tables.jl")
include("scripts/05-Figures.jl")
```

Scripts are idempotent: re-running a script that finds its cached artifact in
`data/` will skip the expensive recompute and use the cache.
