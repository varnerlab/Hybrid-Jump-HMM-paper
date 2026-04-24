# =============================================================================
# 04-Tables.jl
#
# Reads data/results.jld2 and produces the paper tables.
#
# Outputs:
#   paper/sections/tables/table1_aggregate.tex       aggregate scorecard by composer
#   paper/sections/tables/table2_by_branch.tex       hybrid per-construction-flag
#   paper/sections/tables/table3_by_beta_bucket.tex  per-β-quartile breakdown
#
# The LaTeX snippets are hand-written as `tabular` environments ready to be
# wrapped in `\begin{table}` by the paper. Column headers are explicit so the
# tables render cleanly regardless of PrettyTables version.
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

const _PAPER_ROOT      = abspath(joinpath(_ROOT, "..", "paper"))
const _PATH_TO_TABLES  = joinpath(_PAPER_ROOT, "sections", "tables")
isdir(_PATH_TO_TABLES) || mkpath(_PATH_TO_TABLES)

# ── 1. Load artifacts ───────────────────────────────────────────────────────
@info "Loading results + calibration..."
r   = load(joinpath(_PATH_TO_DATA, "results.jld2"))["results"]
cal = load(joinpath(_PATH_TO_DATA, "sim-calibration.jld2"))["calibration"]

cmap = Dict(zip(cal.ticker, zip(cal.beta, cal.r2_real)))
r.beta_cal = [cmap[t][1] for t in r.ticker]
r.r2_cal   = [cmap[t][2] for t in r.ticker]
r.dβ       = abs.(r.β_hat .- r.beta_cal)
r.dR²      = abs.(r.R²_hat .- r.r2_cal)

composer_order_full = ["naive", "gaussian", "hybrid",
                        "residual_jumphmm", "block_bootstrap", "garch_t"]
composer_display    = Dict("naive"            => "Naive",
                            "gaussian"         => "Gaussian SIM",
                            "hybrid"           => "Hybrid",
                            "residual_jumphmm" => "JumpHMM-on-residuals",
                            "block_bootstrap"  => "Block bootstrap",
                            "garch_t"          => "GARCH(1,1)-\$t\$")

# Drop any composers not present in this results DataFrame, so the script
# runs cleanly whether or not 01b / 01c have been executed.
present = Set(unique(r.composer))
composer_order = [c for c in composer_order_full if c in present]
@info "Composers in results" all = composer_order_full present = composer_order

# ── helpers for LaTeX output ────────────────────────────────────────────────
fmt3(x) = @sprintf("%.3f", x)
fmt1(x) = @sprintf("%.1f", x)

function latex_table(io::IO, header::Vector{String}, rows::Vector{Vector{String}};
                     colspec::String = "l" * "r"^(length(header) - 1))
    println(io, "\\begin{tabular}{", colspec, "}")
    println(io, "\\toprule")
    println(io, join(header, " & "), " \\\\")
    println(io, "\\midrule")
    for row in rows
        println(io, join(row, " & "), " \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

# ── 2. Table 1: aggregate scorecard ─────────────────────────────────────────
tbl1 = combine(groupby(r, :composer),
    :dβ      => median => :dβ,
    :dR²     => median => :dR²,
    :ks_p    => (p -> 100 * mean(p .> 0.05)) => :ks_pct,
    :ad_p    => (p -> 100 * mean(p .> 0.05)) => :ad_pct,
    :w1      => median => :w1,
    :kurt    => median => :kurt,
    :hill_up => median => :hill,
)
tbl1 = tbl1[[findfirst(==(c), tbl1.composer) for c in composer_order], :]

println("\n=== Table 1: Aggregate scorecard by composer ===")
show(tbl1, allcols = true, allrows = true); println()

header1 = ["Composer", "\$|\\hat\\beta-\\beta|\$", "\$|\\hat R^2 - R^2_{\\real}|\$",
           "KS pass (\\%)", "AD pass (\\%)", "\$W_1\$", "\$\\kappa\$", "Hill"]
rows1 = [
    [composer_display[tbl1.composer[i]],
     fmt3(tbl1.dβ[i]), fmt3(tbl1.dR²[i]),
     fmt1(tbl1.ks_pct[i]), fmt1(tbl1.ad_pct[i]),
     fmt3(tbl1.w1[i]), fmt3(tbl1.kurt[i]), fmt3(tbl1.hill[i])]
    for i in 1:nrow(tbl1)
]
open(joinpath(_PATH_TO_TABLES, "table1_aggregate.tex"), "w") do io
    latex_table(io, header1, rows1)
end
@info "Wrote paper/sections/tables/table1_aggregate.tex"

# ── 3. Table 2: hybrid per construction flag ────────────────────────────────
hyb = filter(row -> row.composer == "hybrid", r)
tbl2 = combine(groupby(hyb, :flag),
    :ticker  => (t -> length(unique(t))) => :n_tickers,
    :dβ      => median => :dβ,
    :dR²     => median => :dR²,
    :ks_p    => (p -> 100 * mean(p .> 0.05)) => :ks_pct,
    :w1      => median => :w1,
    :kurt    => median => :kurt,
)
flag_display = Dict("HYBRID" => "\\texttt{hybrid}",
                    "HYBRID_CLIPPED" => "\\texttt{hybrid-clipped}",
                    "R2_PRESERVE" => "\\texttt{r2-preserve}")

println("\n=== Table 2: Hybrid composer, by construction flag ===")
show(tbl2, allcols = true, allrows = true); println()

header2 = ["Branch", "\$N_{\\mathrm{tickers}}\$", "\$|\\hat\\beta-\\beta|\$",
           "\$|\\hat R^2 - R^2_{\\real}|\$", "KS pass (\\%)", "\$W_1\$", "\$\\kappa\$"]
rows2 = [
    [get(flag_display, tbl2.flag[i], tbl2.flag[i]),
     string(tbl2.n_tickers[i]),
     fmt3(tbl2.dβ[i]), fmt3(tbl2.dR²[i]),
     fmt1(tbl2.ks_pct[i]),
     fmt3(tbl2.w1[i]), fmt3(tbl2.kurt[i])]
    for i in 1:nrow(tbl2)
]
open(joinpath(_PATH_TO_TABLES, "table2_by_branch.tex"), "w") do io
    latex_table(io, header2, rows2)
end
@info "Wrote paper/sections/tables/table2_by_branch.tex"

# ── 4. Table 3: by-β-quartile breakdown ─────────────────────────────────────
β_vals = cal.beta
qs     = quantile(β_vals, [0.25, 0.5, 0.75])
bucket_for(β) = β ≤ qs[1] ? "Q1 (low \$\\beta\$)" :
                β ≤ qs[2] ? "Q2"                 :
                β ≤ qs[3] ? "Q3"                 : "Q4 (high \$\\beta\$)"
r.beta_bucket = [bucket_for(b) for b in r.beta_cal]

tbl3 = combine(groupby(r, [:beta_bucket, :composer]),
    :dβ      => median => :dβ,
    :ks_p    => (p -> 100 * mean(p .> 0.05)) => :ks_pct,
    :w1      => median => :w1,
    :kurt    => median => :kurt,
    :var_g   => median => :var,
)
# Sort: bucket asc, then composer in canonical order
bucket_order = ["Q1 (low \$\\beta\$)", "Q2", "Q3", "Q4 (high \$\\beta\$)"]
tbl3.bucket_rank   = [findfirst(==(b), bucket_order)   for b in tbl3.beta_bucket]
tbl3.composer_rank = [findfirst(==(c), composer_order) for c in tbl3.composer]
sort!(tbl3, [:bucket_rank, :composer_rank])

println("\n=== Table 3: By β-quartile × composer ===")
show(tbl3[:, Not([:bucket_rank, :composer_rank])], allcols = true, allrows = true); println()

header3 = ["\$\\beta\$ bucket", "Composer", "\$|\\hat\\beta-\\beta|\$",
           "KS pass (\\%)", "\$W_1\$", "\$\\kappa\$", "Var(\$g\$)"]
rows3 = let rs = Vector{Vector{String}}(), prev_bucket = ""
    for i in 1:nrow(tbl3)
        b = tbl3.beta_bucket[i]
        bucket_cell = b == prev_bucket ? "" : b
        prev_bucket = b
        push!(rs, [
            bucket_cell,
            composer_display[tbl3.composer[i]],
            fmt3(tbl3.dβ[i]),
            fmt1(tbl3.ks_pct[i]),
            fmt3(tbl3.w1[i]),
            fmt3(tbl3.kurt[i]),
            fmt1(tbl3.var[i]),
        ])
    end
    rs
end
open(joinpath(_PATH_TO_TABLES, "table3_by_beta_bucket.tex"), "w") do io
    latex_table(io, header3, rows3; colspec = "llrrrrr")
end
@info "Wrote paper/sections/tables/table3_by_beta_bucket.tex"

@info "All tables written to $_PATH_TO_TABLES"
