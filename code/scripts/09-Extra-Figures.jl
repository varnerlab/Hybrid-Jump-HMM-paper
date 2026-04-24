# =============================================================================
# 09-Extra-Figures.jl
#
# Produces the new figures introduced in the peer-review revision. Each block
# is conditional on its input cache existing, so the script runs cleanly at
# any stage of the experimental pipeline.
#
# Outputs (paper/figs/):
#   fig4_baseline_panel.pdf    KS pass and Hill index by composer (all 6)
#   fig5_sensitivity.pdf       R² threshold × floor f heatmap (hybrid KS pass)
#   fig6_stress.pdf            clipping incidence + KS pass vs gm_factor
#   fig7_r2_distribution.pdf   R² histogram on 423-ticker universe
#   fig8_var_backtest.pdf      VaR 99%/95% exceedance rates by composer
#   fig9_synth_tracker.pdf     synthetic-tracker R² recovery
#   fig10_cov_diagnostic.pdf   Cov(ε̃, g_m) / (σ_ε σ_m) histograms
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

const _PAPER_ROOT    = abspath(joinpath(_ROOT, "..", "paper"))
const _PATH_TO_PFIGS = joinpath(_PAPER_ROOT, "figs")
isdir(_PATH_TO_PFIGS) || mkpath(_PATH_TO_PFIGS)

default(
    fontfamily     = "Computer Modern",
    titlefontsize  = 13,
    guidefontsize  = 12,
    tickfontsize   = 10,
    legendfontsize = 10,
    foreground_color_legend = nothing,
    background_color_legend = :white,
    grid           = true, gridalpha = 0.25,
    framestyle     = :box,
)

# Okabe-Ito palette (colourblind-safe)
const OI_BLUE       = RGB(0 / 255,   114 / 255, 178 / 255)
const OI_ORANGE     = RGB(230 / 255, 159 / 255,   0 / 255)
const OI_VERMILLION = RGB(213 / 255,  94 / 255,   0 / 255)
const OI_GREEN      = RGB(0 / 255,   158 / 255, 115 / 255)
const OI_SKY        = RGB(86 / 255,  180 / 255, 233 / 255)
const OI_YELLOW     = RGB(240 / 255, 228 / 255,  66 / 255)
const OI_PURPLE     = RGB(204 / 255, 121 / 255, 167 / 255)

composer_order = ["naive", "gaussian", "hybrid", "residual_jumphmm",
                  "block_bootstrap", "garch_t"]
composer_label = Dict("naive" => "Naive",
                       "gaussian" => "Gaussian SIM",
                       "hybrid" => "Hybrid",
                       "residual_jumphmm" => "JumpHMM-on-resid",
                       "block_bootstrap" => "Bootstrap",
                       "garch_t" => "GARCH(1,1)-t")
composer_color = Dict("naive" => OI_BLUE,
                       "gaussian" => OI_PURPLE,
                       "hybrid" => OI_VERMILLION,
                       "residual_jumphmm" => OI_GREEN,
                       "block_bootstrap" => OI_SKY,
                       "garch_t" => OI_ORANGE)

# ── Figure 4: all-composer headline panel ────────────────────────────────────
results_file = joinpath(_PATH_TO_DATA, "results.jld2")
if isfile(results_file)
    @info "Building Figure 4: baseline panel"
    r = load(results_file)["results"]
    present = [c for c in composer_order if c in unique(r.composer)]

    agg = combine(groupby(r, :composer),
        :ks_p    => (p -> 100 * mean(p .> 0.05)) => :ks_pct,
        :hill_up => median => :hill_med,
        :kurt    => median => :kurt_med,
        :w1      => median => :w1_med,
    )
    agg = agg[sortperm([findfirst(==(c), composer_order) for c in agg.composer]), :]

    # KS pass rate bars
    labels = [composer_label[c] for c in agg.composer]
    colors = [composer_color[c] for c in agg.composer]

    p4a = bar(labels, agg.ks_pct;
              ylabel = "KS pass (%)", title = "Marginal fidelity", legend = false,
              color = colors, ylims = (0, 110), xrotation = 15, fillalpha = 0.85)
    p4b = bar(labels, agg.hill_med;
              ylabel = "Hill (upper 5% tail)", title = "Tail index",
              legend = false, color = colors, xrotation = 15, fillalpha = 0.85)

    fig4 = plot(p4a, p4b; layout = (1, 2), size = (1200, 460),
                left_margin = 7Plots.mm, bottom_margin = 9Plots.mm)
    savefig(fig4, joinpath(_PATH_TO_PFIGS, "fig4_baseline_panel.pdf"))
    @info "Wrote fig4_baseline_panel.pdf"
else
    @info "results.jld2 missing — skipping Figure 4."
end

# ── Figure 5: sensitivity heatmap ────────────────────────────────────────────
R²_grid = [0.70, 0.75, 0.80, 0.85, 0.90]
f_grid  = [0.05, 0.10, 0.15, 0.20, 0.30]
have_sensitivity = all(isfile(joinpath(_PATH_TO_DATA, "results-thresh-$(t)-f-$(f).jld2"))
                       for t in R²_grid, f in f_grid)
if have_sensitivity
    @info "Building Figure 5: sensitivity heatmap"
    KS = zeros(length(f_grid), length(R²_grid))
    for (j, t) in enumerate(R²_grid), (i, f) in enumerate(f_grid)
        d = load(joinpath(_PATH_TO_DATA, "results-thresh-$(t)-f-$(f).jld2"))["results"]
        hy = filter(row -> row.composer == "hybrid", d)
        KS[i, j] = 100 * mean(hy.ks_p .> 0.05)
    end
    p5 = heatmap(string.(R²_grid), string.(f_grid), KS;
                 xlabel = "\$R^2_{\\mathrm{preserve}}\$",
                 ylabel = "idiosyncratic floor \$f\$",
                 title = "Hybrid KS pass rate (%)",
                 color = :viridis, colorbar_title = "KS pass (%)",
                 size = (640, 480), left_margin = 6Plots.mm)
    savefig(p5, joinpath(_PATH_TO_PFIGS, "fig5_sensitivity.pdf"))
    @info "Wrote fig5_sensitivity.pdf"
else
    @info "Sensitivity grid incomplete — skipping Figure 5."
end

# ── Figure 6: clipping stress ────────────────────────────────────────────────
stress_file = joinpath(_PATH_TO_DATA, "results-stress.jld2")
if isfile(stress_file)
    @info "Building Figure 6: clipping stress"
    rs = load(stress_file)["results"]
    hy = filter(row -> row.composer == "hybrid", rs)

    agg_s = combine(groupby(hy, :factor),
        :ks_p => (p -> 100 * mean(p .> 0.05)) => :ks_pct,
        :flag => (flags -> 100 * mean(string.(flags) .== "HYBRID_CLIPPED")) => :clip_pct,
    )
    sort!(agg_s, :factor)

    p6a = plot(agg_s.factor, agg_s.ks_pct;
               xlabel = "market-variance scale factor",
               ylabel = "Hybrid KS pass (%)",
               title = "Marginal fidelity under stress",
               lw = 2.5, color = OI_VERMILLION, marker = :circle, legend = false,
               ylims = (0, 110))
    p6b = plot(agg_s.factor, agg_s.clip_pct;
               xlabel = "market-variance scale factor",
               ylabel = "tickers clipped (%)",
               title = "Clipping incidence",
               lw = 2.5, color = OI_YELLOW, marker = :diamond, legend = false)

    fig6 = plot(p6a, p6b; layout = (1, 2), size = (1100, 420),
                left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
    savefig(fig6, joinpath(_PATH_TO_PFIGS, "fig6_stress.pdf"))
    @info "Wrote fig6_stress.pdf"
else
    @info "results-stress.jld2 missing — skipping Figure 6."
end

# ── Figure 7: R² distribution (supports the structural n=2 argument) ────────
@info "Building Figure 7: R² distribution"
cal = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))["calibration"]
p7 = histogram(cal.r2_real; bins = 50,
               xlabel = "\$R^2_{\\mathrm{real}}\$ (OLS vs SPY)",
               ylabel = "ticker count",
               title = "Calibrated \$R^2\$ across 423 non-SPY tickers",
               color = OI_BLUE, alpha = 0.6, legend = false,
               size = (760, 440), left_margin = 6Plots.mm)
vline!(p7, [0.80]; lw = 2, ls = :dash, color = :black,
       label = "\$R^2_{\\mathrm{preserve}} = 0.80\$")
sorted = sort(cal, :r2_real, rev = true)
for i in 1:2
    annotate!(p7, sorted.r2_real[i] + 0.005, 5 + 4 * i,
              text(sorted.ticker[i], OI_GREEN, 9, :left))
end
savefig(p7, joinpath(_PATH_TO_PFIGS, "fig7_r2_distribution.pdf"))
@info "Wrote fig7_r2_distribution.pdf"

# ── Figure 8: VaR backtest ───────────────────────────────────────────────────
var_file = joinpath(_PATH_TO_DATA, "var-backtest.jld2")
if isfile(var_file)
    @info "Building Figure 8: VaR exceedance"
    v = load(var_file)["results"]
    present = [c for c in composer_order if c in unique(v.composer)]

    agg_v = combine(groupby(v, [:composer, :alpha_level]),
        :rate => mean => :mean_rate,
        :rate => std  => :sd_rate,
        :kupiec_p => (p -> 100 * mean(p .> 0.05)) => :kupiec_pass,
    )

    labels = [composer_label[c] for c in present]
    colors = [composer_color[c] for c in present]

    # 99% VaR panel
    agg99 = filter(row -> row.alpha_level == 0.99, agg_v)
    agg99 = agg99[sortperm([findfirst(==(c), composer_order) for c in agg99.composer]), :]

    p8a = bar(labels, agg99.mean_rate * 100;
              ylabel = "exceedance rate (%)",
              title = "99% VaR: mean exceedance", legend = false,
              color = colors, xrotation = 15, fillalpha = 0.85,
              yerror = agg99.sd_rate * 100)
    hline!(p8a, [1.0]; lw = 2, ls = :dash, color = :black)

    p8b = bar(labels, agg99.kupiec_pass;
              ylabel = "Kupiec pass (%)",
              title = "99% VaR: coverage test",
              legend = false, color = colors, xrotation = 15, fillalpha = 0.85,
              ylims = (0, 110))

    fig8 = plot(p8a, p8b; layout = (1, 2), size = (1200, 460),
                left_margin = 7Plots.mm, bottom_margin = 9Plots.mm)
    savefig(fig8, joinpath(_PATH_TO_PFIGS, "fig8_var_backtest.pdf"))
    @info "Wrote fig8_var_backtest.pdf"
else
    @info "var-backtest.jld2 missing — skipping Figure 8."
end

# ── Figure 9: synthetic tracker recovery ─────────────────────────────────────
synth_file = joinpath(_PATH_TO_DATA, "synth-tracker.jld2")
if isfile(synth_file)
    @info "Building Figure 9: synthetic tracker"
    s = load(synth_file)["summary"]

    p9a = plot(s.r2_true, s.median_r2_hat; group = s.beta_true,
               xlabel = "target \$R^2\$", ylabel = "recovered \$\\hat{R}^2\$",
               title = "R² recovery on synthetic trackers",
               lw = 2, marker = :circle, legend = :bottomright,
               size = (560, 440))
    plot!(p9a, 0.75:0.01:1.0, 0.75:0.01:1.0;
          ls = :dash, color = :black, label = "identity")

    p9b = plot(s.r2_true, s.median_beta_hat; group = s.beta_true,
               xlabel = "target \$R^2\$", ylabel = "recovered \$\\hat{\\beta}\$",
               title = "β recovery on synthetic trackers",
               lw = 2, marker = :diamond, legend = :bottomright,
               size = (560, 440))

    fig9 = plot(p9a, p9b; layout = (1, 2), size = (1150, 460),
                left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
    savefig(fig9, joinpath(_PATH_TO_PFIGS, "fig9_synth_tracker.pdf"))
    @info "Wrote fig9_synth_tracker.pdf"
else
    @info "synth-tracker.jld2 missing — skipping Figure 9."
end

# ── Figure 10: Cov diagnostic ────────────────────────────────────────────────
cov_file = joinpath(_PATH_TO_DATA, "cov-diagnostic.csv")
if isfile(cov_file)
    @info "Building Figure 10: Cov diagnostic"
    dcov = CSV.read(cov_file, DataFrame)
    p10 = plot(; xlabel = "Cov(\\ε̃, gm) / (σ_ε σ_m)", ylabel = "ticker count",
               title = "Empirical cross-term (mean over draws)",
               size = (760, 440), left_margin = 6Plots.mm,
               legend = :topright)
    for name in ("jumphmm_full", "jumphmm_residual", "gaussian")
        sub = filter(row -> row.composer == name, dcov)
        isempty(sub) && continue
        histogram!(p10, sub.cov_norm_mean; bins = 35, alpha = 0.5,
                   label = name, normalize = false)
    end
    savefig(p10, joinpath(_PATH_TO_PFIGS, "fig10_cov_diagnostic.pdf"))
    @info "Wrote fig10_cov_diagnostic.pdf"
else
    @info "cov-diagnostic.csv missing — skipping Figure 10."
end

@info "All available new figures written to $_PATH_TO_PFIGS"
