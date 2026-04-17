# =============================================================================
# 05-Figures.jl
#
# Produces the paper figures from data/results.jld2.
#
# Outputs (paper/figs/):
#   fig1_preservation.pdf   KS pass rate and variance ratio vs β
#   fig2_tails.pdf          excess kurtosis (clipped) and Hill tail index vs β
#   fig3_branch_map.pdf     (β, R²_real) scatter with construction flag
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

const _PAPER_ROOT    = abspath(joinpath(_ROOT, "..", "paper"))
const _PATH_TO_PFIGS = joinpath(_PAPER_ROOT, "figs")
isdir(_PATH_TO_PFIGS) || mkpath(_PATH_TO_PFIGS)

# paper-friendly defaults
default(
    fontfamily    = "Computer Modern",
    titlefontsize = 13,
    guidefontsize = 12,
    tickfontsize  = 10,
    legendfontsize = 10,
    foreground_color_legend = nothing,
    background_color_legend = :white,
    grid           = true,
    gridalpha      = 0.25,
    framestyle     = :box,
)

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading results + calibration..."
r   = load(joinpath(_PATH_TO_DATA, "results.jld2"))["results"]
cal = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))["calibration"]
uni = load(joinpath(_PATH_TO_DATA, "universe.jld2"))

G    = uni["growth_rates"]
tks  = uni["tickers"]
σ²_m = var(G[:, findfirst(==("SPY"), tks)])

cmap_β  = Dict(zip(cal.ticker, cal.beta))
cmap_R² = Dict(zip(cal.ticker, cal.r2_real))
σ²_real = Dict(tks[j] => var(G[:, j]) for j in 1:length(tks))

r.beta_cal = [cmap_β[t]  for t in r.ticker]
r.r2_cal   = [cmap_R²[t] for t in r.ticker]

# ── 2. Per-ticker per-composer summaries ────────────────────────────────────
summary = combine(groupby(r, [:ticker, :composer]),
    :β_hat   => median => :β_hat_med,
    :R²_hat  => median => :R²_hat_med,
    :ks_p    => (p -> mean(p .> 0.05)) => :ks_pass,
    :w1      => median => :w1_med,
    :kurt    => median => :kurt_med,
    :hill_up => median => :hill_med,
    :var_g   => median => :var_med,
)
summary.beta_cal = [cmap_β[t]  for t in summary.ticker]
summary.r2_cal   = [cmap_R²[t] for t in summary.ticker]
summary.var_rel  = summary.var_med ./ [σ²_real[t] for t in summary.ticker]

function composer_frame(df::DataFrame, name::String)
    return filter(row -> row.composer == name, df)
end

# Okabe-Ito colourblind-safe palette
const OI_BLUE       = RGB(0 / 255,   114 / 255, 178 / 255)  # deep blue
const OI_ORANGE     = RGB(230 / 255, 159 / 255,   0 / 255)  # orange
const OI_VERMILLION = RGB(213 / 255,  94 / 255,   0 / 255)  # red-orange (hero)
const OI_GREEN      = RGB(0 / 255,   158 / 255, 115 / 255)  # bluish green
const OI_SKY        = RGB(86 / 255,  180 / 255, 233 / 255)  # sky blue
const OI_YELLOW     = RGB(240 / 255, 228 / 255,  66 / 255)  # yellow
const OI_PURPLE     = RGB(204 / 255, 121 / 255, 167 / 255)  # reddish purple

color_map  = Dict("naive" => OI_BLUE, "gaussian" => OI_PURPLE, "hybrid" => OI_VERMILLION)
marker_map = Dict("naive" => :circle, "gaussian" => :square,   "hybrid" => :diamond)
label_map  = Dict("naive" => "Naive", "gaussian" => "Gaussian SIM", "hybrid" => "Hybrid")

function scatter_composers!(plt, df::DataFrame, ycol::Symbol;
                            alpha::Float64 = 0.4, markersize::Real = 3.5,
                            labels::Bool = true)
    for c in ("naive", "gaussian", "hybrid")
        sc = composer_frame(df, c)
        scatter!(plt, sc.beta_cal, sc[!, ycol];
                 label  = labels ? label_map[c] : nothing,
                 color  = color_map[c],
                 marker = marker_map[c],
                 markersize = markersize,
                 alpha  = alpha,
                 markerstrokecolor = color_map[c],
                 markerstrokewidth = 0.0)
    end
    return plt
end

function binned_median_line!(plt, df::DataFrame, ycol::Symbol;
                             n_bins::Int = 12, lw::Real = 2.5)
    for c in ("naive", "gaussian", "hybrid")
        sc = composer_frame(df, c)
        edges = quantile(sc.beta_cal, range(0.0, 1.0; length = n_bins + 1))
        xs = Float64[]; ys = Float64[]
        for k in 1:n_bins
            lo, hi = edges[k], edges[k + 1]
            mask = (sc.beta_cal .>= lo) .& (sc.beta_cal .<= hi)
            any(mask) || continue
            push!(xs, (lo + hi) / 2)
            push!(ys, median(sc[!, ycol][mask]))
        end
        plot!(plt, xs, ys;
              label = nothing,
              color = color_map[c],
              lw    = lw)
    end
    return plt
end

# ── 3. Figure 1: KS pass rate + variance ratio ──────────────────────────────
@info "Building Figure 1: preservation headline..."

p1a = plot(title = "KS pass rate per ticker",
           xlabel = "calibrated \$\\beta\$",
           ylabel = "KS pass rate (\$\\alpha=0.05\$)",
           ylims = (-0.02, 1.05),
           legend = :outerright)
scatter_composers!(p1a, summary, :ks_pass; markersize = 3.5, alpha = 0.45)
binned_median_line!(p1a, summary, :ks_pass)

p1b = plot(title = "Variance ratio \$\\mathrm{Var}(g)\\,/\\,\\sigma^2_{\\mathrm{gen}}\$",
           xlabel = "calibrated \$\\beta\$",
           ylabel = "\$\\mathrm{Var}(g)\\,/\\,\\sigma^2_{\\mathrm{gen}}\$",
           legend = false)
scatter_composers!(p1b, summary, :var_rel; markersize = 3.5, alpha = 0.45, labels = false)
βs_dense = range(0.0, maximum(summary.beta_cal) * 1.02; length = 200)
σ²_gen_med = median(values(σ²_real))
naive_ref = 1.0 .+ βs_dense.^2 .* σ²_m / σ²_gen_med
plot!(p1b, βs_dense, naive_ref;
      label = nothing, color = OI_BLUE, ls = :dot, lw = 2)
hline!(p1b, [1.0]; label = nothing, color = OI_VERMILLION, ls = :dash, lw = 2)
annotate!(p1b, βs_dense[end-10], naive_ref[end-10] * 1.02,
          text("naive theory \$1+\\rho\$", OI_BLUE, 9, :right))
annotate!(p1b, 0.05, 1.03, text("hybrid target", OI_VERMILLION, 9, :left))

fig1 = plot(p1a, p1b;
            layout = (1, 2), size = (1300, 450),
            left_margin = 6Plots.mm, right_margin = 3Plots.mm,
            bottom_margin = 5Plots.mm, top_margin = 3Plots.mm)
savefig(fig1, joinpath(_PATH_TO_PFIGS, "fig1_preservation.pdf"))
@info "Wrote paper/figs/fig1_preservation.pdf"

# ── 4. Figure 2: kurtosis (clipped) and Hill index ──────────────────────────
#
# Caption caveat (for paper): the "real data" reference line on the left panel
# is the median empirical excess kurtosis across the universe. The hybrid
# median (~7) sits below it (~13) because the JumpHMM marginal is fit to
# preserve the generator's *own* heavy-tailed marginal, which tracks each
# asset's distributional shape but does not replicate the extreme empirical
# 4th moment exactly. The point of the panel is that hybrid tracks the naive
# composition (both inherit the generator tails) while Gaussian SIM collapses
# to κ ≈ 1 regardless of β.
@info "Building Figure 2: tails and kurtosis..."

# real-data reference median
summary.kurt_real = [kurtosis(G[:, findfirst(==(t), tks)]) for t in summary.ticker]
real_kurt_med = median(summary.kurt_real)

p2a = plot(title = "Excess kurtosis vs \$\\beta\$",
           xlabel = "calibrated \$\\beta\$",
           ylabel = "excess kurtosis",
           ylims = (-2, 20),
           legend = :topright)
scatter_composers!(p2a, summary, :kurt_med; markersize = 3.5, alpha = 0.5)
binned_median_line!(p2a, summary, :kurt_med)
hline!(p2a, [real_kurt_med];
       label = "real data (median)",
       color = :black, lw = 2, ls = :dash)

p2b = plot(title = "Hill tail index vs \$\\beta\$",
           xlabel = "calibrated \$\\beta\$",
           ylabel = "Hill index, upper 5 pct tail",
           legend = false)
scatter_composers!(p2b, summary, :hill_med; markersize = 3.5, alpha = 0.5, labels = false)
binned_median_line!(p2b, summary, :hill_med)

fig2 = plot(p2a, p2b;
            layout = (1, 2), size = (1200, 450),
            left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
            top_margin = 3Plots.mm)
savefig(fig2, joinpath(_PATH_TO_PFIGS, "fig2_tails.pdf"))
@info "Wrote paper/figs/fig2_tails.pdf"

# ── 5. Figure 3: branch map in (β, R²) space ────────────────────────────────
@info "Building Figure 3: branch map..."
hyb_summary = composer_frame(summary, "hybrid")
flag_by_ticker = combine(groupby(filter(row -> row.composer == "hybrid", r), :ticker),
    :flag => first => :flag)
flag_lookup = Dict(zip(flag_by_ticker.ticker, flag_by_ticker.flag))
hyb_summary.flag = [flag_lookup[t] for t in hyb_summary.ticker]

p3 = plot(title = "Branch selection in \$(\\beta,\\, R^2_{\\mathrm{real}})\$ space",
          xlabel = "calibrated \$\\beta\$",
          ylabel = "calibrated \$R^2_{\\mathrm{real}}\$",
          ylims = (-0.02, 1.05),
          legend = :topleft, size = (900, 500),
          left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)

flag_color  = Dict("HYBRID" => OI_VERMILLION, "HYBRID_CLIPPED" => OI_YELLOW, "R2_PRESERVE" => OI_GREEN)
flag_label  = Dict("HYBRID" => "hybrid (variance-preserving)",
                   "HYBRID_CLIPPED" => "hybrid-clipped",
                   "R2_PRESERVE" => "\$R^2\$-preserving")
flag_marker = Dict("HYBRID" => :circle, "HYBRID_CLIPPED" => :xcross, "R2_PRESERVE" => :star5)
flag_size   = Dict("HYBRID" => 4, "HYBRID_CLIPPED" => 6, "R2_PRESERVE" => 8)

for flg in ("HYBRID", "HYBRID_CLIPPED", "R2_PRESERVE")
    sub = filter(row -> row.flag == flg, hyb_summary)
    isempty(sub) && continue
    scatter!(p3, sub.beta_cal, sub.r2_cal;
             label  = "$(flag_label[flg]) (\$n=$(nrow(sub))\$)",
             color  = flag_color[flg],
             marker = flag_marker[flg],
             markersize = flag_size[flg],
             alpha = 0.55,
             markerstrokecolor = flag_color[flg],
             markerstrokewidth = 0.0)
end
hline!(p3, [0.80]; label = "\$R^2_{\\mathrm{preserve}} = 0.80\$",
       lw = 2, ls = :dash, color = :black)

# label the two trackers
for row in eachrow(filter(r -> r.flag == "R2_PRESERVE", hyb_summary))
    annotate!(p3, row.beta_cal + 0.02, row.r2_cal,
              text(row.ticker, OI_GREEN, 10, :left))
end

savefig(p3, joinpath(_PATH_TO_PFIGS, "fig3_branch_map.pdf"))
@info "Wrote paper/figs/fig3_branch_map.pdf"

@info "All figures written to $_PATH_TO_PFIGS"
