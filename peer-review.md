# Simulated Peer Review: *Variance-Corrected Single-Index Composition for Heavy-Tailed Return Generators*

**Manuscript:** `paper/Paper_v1.tex` (commit `dcacfdb`, 2026-04-24).
**Assumed target venue:** mid-to-high-tier quantitative-finance / empirical-finance venue (Journal of Financial Econometrics, Quantitative Finance, Journal of Empirical Finance) or an arXiv preprint for the ML-adjacent audience implied by the NeurIPS template.
**Reviewer-persona note:** the skill prompt describes metabolic-engineering reviewers. I have substituted domain-qualified personas (quantitative finance, factor models, competing synthetic-return generators) while preserving the three review levels (moderate / hard / very hard).

---

## Reviewer 1 (Moderate) — Quantitative-finance methods researcher

### Summary
The paper introduces a closed-form variance correction for composing per-asset heavy-tailed return generators with a single-index market factor. The correction rescales the idiosyncratic draw by $s_i = \sqrt{1 - \rho_i}$ so that the composed series inherits the generator's marginal variance exactly while OLS recovers the calibrated $\beta$. The method is clean, the experimental protocol is sensible, and the result is a useful addition to the practitioner's toolkit; I recommend minor revision.

### Strengths
1. **Transparent algebra.** The variance-correction derivation (Eqs.~3–5) is a one-line calculation and the resulting scalar is closed-form. The fact that the SPY-on-itself identity falls out as the $R^2 \to 1$ degenerate limit (method.tex:223–227) without special-casing is a small but satisfying property.
2. **Auditability.** The persisted construction flag (`hybrid`, `hybrid-clipped`, `r2-preserve`) plus the preserved pair $(\beta_i, \beta_i^{\eff})$ under clipping (method.tex:178–179) is exactly what a downstream user needs to trust the output.
3. **Paired-innovation protocol.** Using the same $\teps_i^{(r)}$ across all three composers for each (ticker, replication) pair (results.tex:40–46) is a good experimental design choice — the metric deltas isolate the composition rule from generator stochasticity.
4. **Clear separation of failure modes.** Figure~2 (tails + Hill) and Figure~1 (KS + variance ratio) decompose the two baselines' failure modes cleanly: naive inflates variance but keeps tails; Gaussian-SIM recovers $R^2$ but flattens tails.
5. **Honest scope statement.** The Discussion (discussion.tex:36–64) is explicit about what the evaluation does and does not cover (deferred copula reorder, fixed SPY path, high-turnover evaluation).

### Weaknesses
1. **Clipping branch is untested.** Method-section text devotes a full subsection to `hybrid-clipped` (method.tex:148–179), but results.tex:143–147 reports that clipping "did not trigger anywhere in the universe." The entire clipping branch is therefore unvalidated on the experimental universe. This is fine if framed as forward-looking, but the paper currently presents clipping as a first-class design element. *Fix:* add a synthetic stress experiment with an elevated-volatility $g_m$ (e.g., $2\times$ or $3\times$ SPY sample variance) so at least Q4 tickers trigger clipping, and report the clipped $\beta_i^{\eff}$ recovery.
2. **R²-preserve branch has $n = 2$.** Only QQQ and SPYG cross the $R^2 = 0.80$ threshold (results.tex:135–142). Pitching r2-preserve as a major branch with two data points is thin. *Fix:* broaden the universe to include a dozen or two index ETFs (sector SPDRs, iShares country ETFs, bond ETFs) where $R^2$ preservation is the natural target. Report per-ticker KS and $|\hat R^2 - R^2_{\real}|$ for that expanded tracker subset.
3. **No uncertainty quantification on aggregate pass rates.** The headline number "KS pass 50.4% vs 4.5% naive vs 0.6% Gaussian" (abstract.tex:29–31, results.tex:94–99) is a point estimate over 42,300 observations but with a single RNG seed (results.tex:76). A previous draft claimed "standard errors on aggregate pass rates below 0.3%" (git-history removed in 2ed6645), but no uncertainty figure is in the current version. *Fix:* re-run with 8–16 seeds and report aggregate-pass-rate $\pm$ seed-to-seed standard deviation, or bootstrap CIs on the within-seed ticker distribution.
4. **Variance-ratio denominator uses real variance as proxy.** results.tex:70–72 states that the ratio $\mathrm{Var}(g_i^{(r)}) / \sigma_{\gen,i}^2$ is computed "using the real-data variance as a proxy for the generator's targeted variance." This is a subtle choice: if the JumpHMM generator's sample variance differs materially from the real-data variance (it usually does, since JumpHMM targets the distribution, not a variance constraint), the reported ratio mixes composition-rule error with generator-miscalibration error. *Fix:* report both (a) $\mathrm{Var}(g)/\sigma_{\gen}^2$ using the actual generator sample variance averaged over replications, and (b) $\mathrm{Var}(g)/\mathrm{Var}(r_{i,\real})$ using the real-data variance. Clarify which one Figure~1 (right) uses.
5. **Typo/copyedit in intro.** introduction.tex:17 reads "Theses two requirements" — should be "These".

### Questions for Authors
1. In the abstract you report a median absolute $\beta$ error of $0.018$. All three composers recover $\beta$ to within $0.018$–$0.022$ (Table~1). How much of that $0.018$ is composition-rule error versus finite-$T$ sampling noise in OLS on $T = 2767$ days? A calibration line showing $|\hat\beta - \beta|$ vs $1/\sqrt{T}$ under each composer would be informative.
2. The floor $f = 0.10$ is treated as a tunable knob (method.tex:175). What is the sensitivity of the reported aggregate metrics to $f \in \{0.05, 0.10, 0.20\}$ on the subset of tickers where clipping would trigger under an elevated-volatility $g_m$?
3. The Hill tail index is reported at 0.37 (hybrid/naive) vs 0.22 (Gaussian) — please clarify which Hill estimator convention (tail index $\theta$, or reciprocal shape $1/\alpha$) and what $k$ / fraction of the upper tail is used.

### Requested Experiments / Analyses
1. **Clipping stress test.** Run the 424-asset universe with $g_m$ constructed by scaling the SPY growth-rate path to $2\times$ and $3\times$ its sample variance. Report per-ticker $|\hat\beta^{\eff} - \beta^{\eff}_{\text{expected}}|$ for clipped tickers and the consequences for the marginal.
2. **Seed-to-seed variability.** Aggregate metrics with 8–16 RNG seeds. Report means $\pm$ SD (or 90% CI).
3. **Expanded tracker subset.** At least 15–20 additional index/sector ETFs in the r2-preserve bucket to give that branch statistical weight.

### Minor Comments
- Abstract ("the hybrid construction recovered the factor loading to a median absolute error of $0.018$") — median of what (across tickers? across replications? across both?). One extra clause would clarify.
- Algorithm~1 lists a `fallback` flag (method.tex:245) but the algorithm never assigns it and no text describes when `fallback` is chosen. Either define or remove.
- Figure~2 (tails): $y$-axis clip $[-2, 20]$ is noted in caption but not justified; any points clipped out?
- Table~3: the "Var($g$)" column — is that the median of per-ticker sample variances across the bucket, averaged over replications?

### Recommendation
**Minor revision.** The core contribution is sound and well presented. The three requested experiments are incremental and directly address the weaknesses above without requiring a redesign.

---

## Reviewer 2 (Hard) — Factor-model / empirical-finance expert

### Summary
The paper proposes a rescaling scalar $s_i^2 = 1 - \rho_i$ applied to a per-asset generator draw to preserve marginal variance under SIM composition. The derivation is correct under the stated independence assumption but the evaluation leaves important statistical claims underpowered, several design choices are unjustified, and alternative approaches are not seriously benchmarked. I recommend major revision.

### Strengths
1. The preservation-criteria decomposition (method.tex:74–104) into SIM recovery, marginal fidelity, and cross-sectional dependence is the right conceptual frame.
2. Sign-preserving clip (method.tex:164) correctly handles negative-$\beta$ assets; this is the kind of edge case that silently breaks lesser implementations.
3. The $R^2$-target identity (Eq.~8) and its $R^2 = 1$ limit are elegant.

### Weaknesses
1. **The "fit to residuals" alternative is not benchmarked.** The paper's setup is: fit the per-asset generator to the asset's *full* return history, then apply a variance correction because the generator output double-counts the market component. The natural alternative — fit the generator to the OLS residual $e_i(t) = r_i(t) - \hat\alpha_i - \hat\beta_i r_m(t)$ directly, then use $s_i = 1$ — is mentioned nowhere in the paper. Without this baseline the evaluation omits the most obvious competing recipe. *Fix:* add a fourth composer, **"residual-JumpHMM"**, that fits JumpHMM to $e_i$ rather than $r_i$. Report all the same metrics. If hybrid beats residual-JumpHMM on some metric, that is the actual contribution; if not, the paper's scope narrows.
2. **Independence assumption between $\teps_i$ and $g_m$.** method.tex:57–63 asserts $\Cov[\teps_i, g_m] \approx 0$ "by construction" because JumpHMM is sampled independently of the market path. That is correct for the *sampled* $\teps_i$, but the generator itself was fit to $r_i$, which has sample covariance with $r_m$ of magnitude $\beta_i \sigma_m^2$. Any temporal structure in $\teps_i$ that reflects this dependence (regime states, jump timing) is then injected alongside $g_m$. The paper should report the empirical $\widehat{\Cov}(\teps_i, g_m)$ distribution across the 423 tickers, not just assert the population-level approximation holds.
3. **Single seed; no reported variability.** results.tex:76 reports seed 1234. Aggregate numbers in Table~1 and Table~3 are presented to three significant figures (50.4%, 47.6%, 0.365), implying a precision that a single-seed run does not support. With 42,300 observations per composer this is probably small-SE, but "probably small-SE" is not a substitute for reporting SE. *Fix:* report confidence intervals, either via seed replication or bootstrap across tickers within a seed.
4. **$R^2_{\mathrm{preserve}} = 0.80$ is data-peeked.** The text explicitly says the threshold was chosen by looking at where the universe's cumulative $R^2$ distribution "sorts itself cleanly" (method.tex:197–200). That is a legitimate engineering choice but it is *not* a prospectively defensible threshold. *Fix:* present the branch assignment under a grid of thresholds (0.70, 0.75, 0.80, 0.85, 0.90) and show that the main metric ordering (hybrid > naive > Gaussian SIM on KS, etc.) is insensitive to the threshold within a reasonable range. If it is sensitive, acknowledge and discuss.
5. **Rescale factor > 1 on the $R^2$-preserve branch.** method.tex:227–231 notes that the rescale factor in Eq.~(10) can exceed one when the synthetic market has lower variance than the calibration market, "stretching the generator's tails to fit a target larger than the original draw." This is flagged as an auditable case but is neither evaluated nor bounded. On which fraction of plausible synthetic-market variances does this trigger? What is the distortion to tail statistics when it does?
6. **Floor $f = 0.10$ is arbitrary.** method.tex:175 states $f = 0.10$ "guarantees that at least 10% of each asset's variance comes from idiosyncratic noise." That is a definition of $f = 0.10$, not a justification. Why 0.10 and not 0.05 or 0.20? The answer likely involves a tradeoff between clipping aggressiveness and marginal fidelity at high $\rho_i$; the paper should plot that tradeoff explicitly.
7. **Cross-sectional fidelity is deferred but still listed as a preservation criterion.** method.tex:97–104 lists cross-sectional dependence as the third preservation criterion, then results.tex:56–60 turns off the rank reorder and defers this question to the companion paper. This creates a narrative mismatch: the reader is told three properties matter; two are evaluated. Either (a) evaluate the third with at least a sanity check (e.g., rank-correlation matrix recovery under empirical-copula reorder, even if abbreviated), or (b) remove cross-sectional dependence from the preservation-criteria list and present it only as a downstream concern.

### Questions for Authors
1. Why is the residual-JumpHMM baseline excluded? A reader's first instinct is to avoid double-counting at the generator-fitting stage rather than correct for it post-hoc.
2. The aggregate KS pass rate is 50.4%. In an unconditional KS test at $\alpha = 0.05$ with $T \approx 2767$ against a perfectly matched distribution, pass rate should be 95%. A pass rate of 50% means the hybrid marginal rejects the real marginal at the 5% level for roughly half the tickers. Is the argument that 50% is "good enough" for downstream use, or that the residual 50% failure reflects JumpHMM miscalibration rather than composition-rule error? How do you separate these?
3. On the $R^2$-preserve branch, why rescale $\teps_i$ by a scalar rather than draw a Gaussian with exactly the target residual variance? If the whole point of that branch is to target $R^2$ (not the marginal), the JumpHMM draw loses its motivation.

### Requested Experiments / Analyses
1. **Residual-JumpHMM baseline.** Fit JumpHMM to OLS residuals, compose with $s = 1$, report all the same metrics. This is non-negotiable for publication.
2. **Threshold sensitivity.** Report Table~1 under $R^2_{\mathrm{preserve}} \in \{0.70, 0.75, 0.80, 0.85, 0.90\}$.
3. **Floor-$f$ tradeoff curve.** Under the elevated-volatility $g_m$ of Reviewer 1's request, sweep $f \in \{0.05, 0.10, 0.15, 0.20, 0.30\}$ and plot KS pass rate against clipping incidence.
4. **Empirical cross-term diagnostic.** Histogram of $\widehat{\Cov}(\teps_i, g_m)/(\sigma_{\gen,i}\sigma_m)$ across the 423 tickers, per composer, to show the independence assumption is empirically satisfied.

### Minor Comments
- results.tex:82–83: "$423 \times 100 \times 3 = 126{,}900$ scored series" — arithmetic check: correct.
- The text calls $s_i$ a "scalar" throughout but it is time-invariant per ticker per replication; worth saying so once.
- Algorithm~1 line 15 ("Optional") uses English italic rather than a standard algorithmic convention; minor.
- The variance ratio numbers in Table~3 (column "Var($g$)") have no target column; the reader has to compute what "unit preservation" looks like in that column.

### Recommendation
**Major revision.** The residual-JumpHMM baseline omission is the decisive issue; without it the paper's claim of a novel tradeoff-resolving recipe is not defended against the most obvious alternative. The threshold sensitivity, floor sensitivity, and uncertainty quantification are all addressable without restructuring but are not optional for a methods paper in this venue.

---

## Reviewer 3 (Very Hard) — Author of competing synthetic-return generators

### Summary
The paper presents a one-line algebraic rescaling and claims it "closes a gap" in SIM composition. The rescaling is correct, but the framing as a novel contribution substantially overclaims relative to the decades-long literature on factor-residual modeling, and the evaluation is structured in a way that the contribution cannot fail: the two chosen baselines are both strawmen. I recommend reject (at current scope) or major revision with a substantially broader benchmark and sharpened novelty claims.

### Strengths
1. The paper is clearly written and the algorithm is reproducible.
2. The $R^2$-preserving branch's SPY-on-itself degenerate limit is a nice property.
3. The construction flag auditability is good engineering practice.

### Weaknesses
1. **The baselines are strawmen.** "Naive" composition (plug the generator draw in as the residual) is not a recipe any careful practitioner uses; it is the obvious wrong answer. "Gaussian-residual SIM" is not an alternative *to* the hybrid, it is a different method targeting a different quantity ($R^2$ rather than marginal). With these two baselines the hybrid has no real competitor. *Fix:* benchmark against at least three legitimate alternatives from the literature: (a) GARCH-$t$ residuals fit per ticker to OLS residuals (Bollerslev; already in your references but not used as a baseline), (b) block-bootstrap residuals from joint history (Politis–Romano; already cited but not used), (c) multivariate $t$-copula on OLS residuals fit jointly. All three address the same marginal-vs-factor problem; the paper's 424-asset evaluation would be genuinely informative against them.
2. **The contribution is mathematically trivial.** The identity $\Var[\beta g_m + s\teps] = \beta^2\sigma_m^2 + s^2\sigma_{\gen}^2 = \sigma_{\gen}^2$ solving to $s^2 = 1 - \rho$ is one line of undergraduate probability. The paper even acknowledges this (discussion.tex:29–30: "No part of the construction is new in a deep sense"). A methods paper at a serious venue needs a non-trivial contribution. The candidates are (i) the branching logic (clipping + R²-preserve), (ii) the SPY-limit elegance, (iii) the empirical demonstration. (i) is engineering, (ii) is a one-paragraph remark, (iii) is weak given the strawman baselines.
3. **Novelty relative to the companion paper is unclear.** discussion.tex:46–48 defers cross-sectional fidelity to "the companion JumpHMM.jl paper" [Alswaidan-Varner 2026]. Readers cannot tell what this paper adds on top of that companion. If the companion paper already uses JumpHMM + SIM composition, the present paper's incremental contribution needs to be stated precisely with a one-paragraph diff, otherwise the two papers will read as the same work split for LPU.
4. **The Discussion's new bias/autocorrelation caveat undercuts the paper's production-readiness claim.** The abstract and introduction sell this as a "practitioner-ready recipe" for workflows where "both tail shape and market dependence matter." The Discussion (discussion.tex:112–128) now admits that the iid residual structure creates a large-magnitude diversification-return artifact under daily rebalancing, and that the fix requires AR(1)/GARCH/block-bootstrap — i.e., *not* the recipe this paper presents. What downstream workflow, concretely, can use the paper's output as a drop-in without running into the caveat? If the answer is "buy-and-hold and low-turnover strategies," the introduction should say so.
5. **Data-peeking on the threshold.** Already raised by R2. Worth emphasizing that the paper describes picking $R^2_{\mathrm{preserve}} = 0.80$ by inspection of the universe's own $R^2$ distribution (method.tex:197–200). In a different universe (sector-heavy, single-country, fixed-income), this threshold will be different. The algorithm is advertised as metadata-free; in practice it has a hand-tuned hyperparameter fit to this universe.
6. **Only 2 trackers in r2-preserve.** Generalizing a branch that governs two data points into a method-section subsection is not proportionate.
7. **No downstream validation.** The paper proposes composed synthetic paths for "portfolio construction, risk attribution, and stress testing" (introduction.tex:4–5). None of these downstream tasks is evaluated. A VaR backtest, a portfolio-optimization out-of-sample comparison, or a stress-test scenario accuracy would turn "50.4% KS pass rate" into a consequential number. Without that, the paper is a method-fit exercise.
8. **Kurtosis claim is weaker than presented.** The hybrid kurtosis is $7.2$ vs real-data $\sim 13$ (results.tex:90–91). That is 55% of the real value. The paper's framing ("hybrid and naive track the generator's kurtosis") is accurate but deflects from the fact that neither hybrid nor naive recovers real-data kurtosis; both fall meaningfully short because the JumpHMM marginal itself does (Figure~2 caption, Paper_v1.tex:202–205).

### Questions for Authors
1. What does your method add over fitting a per-asset GARCH-$t$ to OLS residuals and composing with $s = 1$? Please benchmark.
2. What precisely does this paper claim that the Alswaidan-Varner (2026) companion paper does not?
3. At what target fraction of real-data kurtosis (80%? 90%? 100%?) would you call the method successful? The current 55% recovery is not a headline result.
4. If a user wants to evaluate a daily-rebalanced allocator on the composed paths, is the recommendation to use this paper's method and apply a post-hoc bias correction externally, or to use a different generator entirely? Either answer narrows the paper's claimed use case.

### Requested Experiments / Analyses
1. **Real-alternative benchmarks.** GARCH-$t$-on-residuals, block-bootstrap-on-residuals, $t$-copula on residuals. Drop the Gaussian-SIM strawman or retain it only as a reference.
2. **One downstream-use validation.** A historical VaR backtest on the composed paths vs real paths, or a historical portfolio-optimization out-of-sample comparison. Pick one; it should tell the reader which real problem the method solves.
3. **Precise novelty diff against the companion paper.** A paragraph in the Introduction that names the exact claims in Alswaidan-Varner 2026 and the exact claims here, with no overlap. If overlap is unavoidable, consider merging.

### Minor Comments
- The abstract ends with "All artifacts, scripts, and the fitted marginal cache are released with the paper" — verify that the cache URL resolves at submission time.
- introduction.tex:17 "Theses" → "These".
- The unnumbered footnote to JumpHMM.jl (method.tex:47) duplicates the code-availability URL in the main .tex. One or the other.
- Table~2: the "r2-preserve" branch row reports 99% KS pass rate over $n = 2$ tickers; that number has no statistical weight. Either drop the row or annotate the $n$.

### Recommendation
**Reject, with encouragement to resubmit after (a) adding real alternative baselines, (b) establishing a non-trivial, clearly-scoped novelty claim distinct from the companion paper, and (c) providing at least one downstream-use validation.** The mathematical core is correct but the present scope does not clear the bar for an independent methods contribution.

---

## Summary of Actionable Items

Consolidated and prioritized across the three reviews.

### Priority 1 — blocking for publication
1. **Add a "residual-JumpHMM" baseline.** Fit JumpHMM to OLS residuals, compose with $s = 1$, report all metrics. *(R2, R3)*
2. **Add a legitimate competing baseline from the cited literature.** GARCH-$t$-on-residuals (Bollerslev) and/or block-bootstrap-on-residuals (Politis–Romano). Drop "Gaussian SIM" as the primary baseline; keep as a reference if desired. *(R3)*
3. **Uncertainty quantification on aggregate metrics.** Re-run with 8–16 seeds OR bootstrap CIs on ticker distribution. Report ±SD or 90% CI on all three-sig-fig numbers. *(R1, R2)*
4. **Validate the clipping branch.** Run a synthetic stress test with elevated-volatility $g_m$ (e.g., $2\times$ and $3\times$ SPY sample variance) that actually triggers clipping on Q4 tickers. Report clipped-$\beta^{\eff}$ recovery and marginal consequences. *(R1, R2)*
5. **Sharpen novelty claim vs. Alswaidan-Varner (2026) companion paper.** One paragraph in Introduction naming the exact non-overlapping contribution. *(R3)*

### Priority 2 — substantially strengthens the paper
6. **Expand the r2-preserve tracker subset.** 15–20 index/sector ETFs, bond ETFs, country ETFs. Report per-ticker metrics. *(R1, R3)*
7. **Threshold sensitivity for $R^2_{\mathrm{preserve}}$.** Grid of $\{0.70, 0.75, 0.80, 0.85, 0.90\}$; show main ordering is insensitive. *(R2, R3)*
8. **Floor-$f$ sensitivity / tradeoff curve.** Sweep $f$ on the clipping-triggered subset; plot KS pass rate vs clipping incidence. *(R2)*
9. **One downstream-use validation.** Historical VaR backtest OR out-of-sample portfolio optimization comparison. *(R3)*
10. **Empirical $\widehat{\Cov}(\teps_i, g_m)$ histogram.** Verify the independence assumption holds empirically across the 423 tickers. *(R2)*

### Priority 3 — clarifications and minor
11. Clarify variance-ratio denominator: generator sample variance vs real-data variance. Report both if needed. *(R1)*
12. Clarify Hill-estimator convention ($\theta$ vs $1/\alpha$, choice of $k$). *(R1)*
13. Define or remove the `fallback` construction flag in Algorithm~1. *(R1)*
14. Fix typo "Theses" → "These" at introduction.tex:17. *(R1, R3)*
15. Deduplicate the JumpHMM.jl URL (footnote in method vs code-availability block). *(R3)*
16. Annotate $n = 2$ on the r2-preserve row of Table~2. *(R3)*
17. Add a sentence on computational cost (fit time, generation time for the 424-asset × 100-replication × 2767-day run).
18. Restate scope in the Introduction after the new Discussion caveat: name which downstream use cases the method supports as-is and which require the residual-autocorrelation extension. *(R3)*

### Aggregate recommendation
Two of three reviewers recommend revision (minor and major); one recommends reject-with-resubmit. The mathematical core is not in dispute; the decisive issues are (a) baseline scope (R2, R3), (b) statistical rigor on reported numbers (R2), and (c) novelty framing relative to the companion paper (R3). Addressing Priority 1 items 1–5 should convert this to a revise-and-resubmit with a realistic path to acceptance.
