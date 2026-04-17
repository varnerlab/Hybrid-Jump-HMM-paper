# Journal submission prospects — hybrid-jump-hmm paper

Date: 2026-04-17

## Short verdict

The paper is journal-submittable — not arXiv-only — but PLoS ONE is not
the right home. A practitioner-oriented quantitative-finance methods
journal is a much better fit.

## Why not PLoS ONE

PLoS ONE is "sound science" in principle and broad-scope in scope, but
in practice ~90% of what they publish is life-sciences. Their editorial
pool is thin on quantitative finance, review quality is unpredictable,
and the APC is not cheap. The paper would get a journal stamp but not a
readership that cares about SIM composition or synthetic return
generation.

## Better-fit venues (in rough order of match)

1. **Journal of Financial Data Science** — short practitioner-oriented
   methods papers are the core format; reproducibility and released
   artifacts are a plus rather than a curiosity. Most natural home.
2. **Journal of Computational Finance** — fits the algorithmic /
   methodology framing, tolerates elementary math when paired with
   clean empirics.
3. **Journal of Risk** — works if the paper is reframed around
   stress-test use (the clipping story plus the "scenario markets more
   volatile than calibration" point in the discussion).
4. **Quantitative Finance** — possible, but reviewers there often want
   either more theory or broader empirics than the current draft
   delivers.
5. **Journal of Statistical Software / SoftwareX** — would work if the
   paper is repackaged as a software paper around `JumpHMM.jl` plus a
   sibling composer package. Not the current framing but an option.

Not in scope: high-impact general-finance journals (JFE, RFS, JF) — the
contribution is too narrow and the math too elementary for those
venues.

## Likely referee pushbacks

Worth pre-empting in the current draft:

1. **The math is elementary.** A one-line variance correction, which
   the introduction already concedes is not new in a deep sense. The
   defense is practical impact (50.4% vs 4.5% vs 0.6% KS pass rates is
   a decisive improvement), the auditable per-asset construction flag,
   and the principled treatment of the degenerate tracker limit.
2. **Single universe, single market, single window.** A reviewer will
   likely ask for at least an out-of-sample robustness check or a
   second asset class (international equities, crypto, commodities, or
   a scenario-market run with $g_m$ more volatile than the calibration
   market — exactly the regime the discussion flags as where clipping
   matters).
3. **Heavy reliance on the companion JumpHMM.jl paper.** The paper
   needs to land standalone, which the current introduction mostly
   achieves but the Results section leans on.

## Suggested sequencing

1. Post to arXiv first under **q-fin.CP** (computational finance) — the
   right primary category; **q-fin.RM** (risk management) is a
   reasonable cross-listing. Establishes priority and invites informal
   feedback before a formal review.
2. Submit to **Journal of Financial Data Science** as the first
   journal target. Short-paper format fits the scope, turnaround is
   reasonable, and the audience is right.
3. If rejected there, Journal of Computational Finance is the natural
   second target.
4. Keep PLoS ONE as a distant fallback only if the practitioner
   journals fail and a journal stamp is needed for tenure/grant
   reasons.

## What would strengthen a journal submission

If the author wants to harden the paper before submission rather than
submit as-is:

- Add a second universe or a scenario-market run to address the
  single-dataset concern (and to let clipping actually trigger, which
  it does not in the current run).
- Extend the lagged-factor SIM sketch in the discussion into a short
  additional subsection with one worked example.
- Consider a cross-sectional dependence experiment using the
  empirical-copula reorder, even a short one, to close the third
  preservation criterion with data rather than leaving it deferred to
  the companion paper.

None of these are required to submit; they are what would move the
paper from "accept with minor revisions" to "accept" at a methods
journal.
