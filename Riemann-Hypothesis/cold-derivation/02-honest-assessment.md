---
title: "Riemann Cold Derivation — Honest Assessment"
date: 2026-04-28
companion: 01-riemann-cold.md
purpose: Separate what was derived from what was asserted. No varnish.
---

# Honest Assessment of the Cold Derivation

This is the companion to `01-riemann-cold.md`. The derivation is what we wish were true; this is what we know is true about the derivation.

## What the derivation actually is

A **physical reading** of the Riemann functional equation in RS2 language. Not a proof of RH. Not a quantitative prediction. A structural correspondence with a physical interpretation of "½" — and that is the entire RS2 contribution.

If we strip out RS2 vocabulary, what remains is:
- The functional equation $\xi(s) = \xi(1-s)$ (Riemann, 1859, derived).
- The Berry–Keating operator $H = -i(xd/dx + 1/2)$ as a candidate Hilbert–Pólya operator (Berry & Keating, 1999, postulated).
- The critical line σ = ½ (open).

What RS2 adds, in honest accounting:
- A **physical interpretation** of why the operator has +½ rather than, say, +1 or 0. The interpretation: birotation ≡ metaplectic double cover, and metaplectic representations carry weight ½.

That is the whole new content. Everything else in §1–§9 is reorganization or notation.

## Load-bearing assertions

These are claims the derivation makes that are *not* derived inside RS2:

### LB1: RS2 birotation IS the metaplectic double cover (§5.2)

This is the single load-bearing identification. The derivation argues for it on three grounds:

1. Both are double covers of the same projective base.
2. Both halve natural angular periods.
3. RS2's spin-½ assignment to magnetic rotation matches metaplectic spin-½ representations.

These are *structural similarities*, not a proof of identity. Mathematicians do not use the word "birotation"; Peret does not use the word "metaplectic." We are asserting the bridge.

If LB1 is wrong, the +½ in (6.1) loses its physical origin and the cold derivation collapses to standard Berry–Keating.

LB1 is, however, **falsifiable** in a sharp sense. The metaplectic group is *unique* — there is exactly one connected double cover of SL(2,ℝ). If RS2 birotation has the same algebraic data (involution, double cover structure, half-period behavior under the relevant action), then they must be the same up to isomorphism. The derivation's burden is to show RS2 birotation has this data; we have argued plausibly but not proved it.

### LB2: The σ-coordinate map (4.1) is canonical

We wrote $\sigma = \tfrac{1}{2} + \tfrac{r}{2}$ and asserted this is the natural map from RS2 displacement to the Riemann coordinate. The justification (the half-shift is forced by the metaplectic weight ½) is *retrofitted*: we identified the half-shift by what makes the involution work, then claimed the identification was forced.

A more honest order: we *observed* that the coordinate map putting unit speed at σ = ½ reproduces the functional equation involution; we did not derive the map from RS2 axioms alone.

The map is consistent and uniquely determined by the constraint "match the involution s ↔ 1−s onto the RS2 reciprocal r ↔ −r with unit speed at the fixed point." But that constraint *is* the conclusion we want; using it to derive the map is circular if we want to derive the conclusion from the map.

The non-circular claim: the map (4.1) is the unique linear map ℝ → ℝ that satisfies (i) involution-preserving and (ii) unit speed at the fixed point. The "½" appears only because we *demanded* unit speed at the fixed point of the standard Riemann involution. An alternative coordinate could be devised, but not without a different physical-frame anchoring.

### LB3: The Hilbert–Pólya operator construction is open

§6 sketches a structural argument for "why σ = ½" but does not close. The Hilbert–Pólya program is widely understood to be the canonical bottleneck for proving RH. RS2 does not solve it. The +½ in (6.1) is now physically interpreted, but the construction of boundary conditions / domain such that (6.1) has eigenvalues at the imaginary parts of the zeros is not given.

To be explicit: this derivation, even granted LB1 and LB2, does **not** prove RH.

## Soft claims

These are claims the derivation makes that are weaker than they read on first encounter:

### SC1: §6.3 unit-boundary heuristic is heuristic only

"Zeros ought to be at the natural datum because zeros are resonances and resonance-level destructive interference requires scaling balance" is a physical motivation, not a constraint. The destructive-interference argument has the structure of a wishful explanation: it predicts what we want and does not falsify the alternative (zeros could exist as quartets, the interference could partially cancel, etc.).

### SC2: §7.1 trivial-zero quaternion reading is qualitative

The trivial zeros at −2, −4, −6, … fitting the −1 quaternion direction × magnetic 2D doubling is a numerical coincidence reframed as a structural identity (compare honest-flag #3 in `00-FIRST-PRINCIPLES.md`). It does not add quantitative content; it could be reverse-engineered from any coarse RS2 mapping.

### SC3: §7.2 GUE statistics is consistent, not predicted

The argument "RS2 birotation breaks time-reversal symmetry, hence GUE rather than GOE" is consistent with Montgomery–Odlyzko numerics but does not predict the GUE statistics from RS2 — the statistics were known empirically before this argument, and RS2 is being adapted to fit.

### SC4: §7.3 4n² prediction is speculative

The suggestion that low-lying zeros show 4n²-related structure is the kind of post-hoc numerology RS2 has historically been criticized for (see Peret's own RS2-105 notes on quantum π = 4 vs. analog π and the reference-system relativity). It is offered as a falsifiable check, not as a load-bearing prediction. Numerical investigation should proceed cautiously.

## What would falsify the cold derivation

The derivation collapses or weakens substantially if any of:

| Falsifier | Effect |
|---|---|
| RS2 birotation does *not* have metaplectic algebraic structure | LB1 fails; +½ loses physical origin; reduces to Berry–Keating |
| Spin-½ in RS2-107 has a different group-theoretic origin than metaplectic ½ | LB1 weakened; alternative origin needed |
| Numerical investigation of zero spacings shows no 4n²-related structure under RS2 displacement binning | SC4 falsified — derivation survives but loses one falsifiable check |
| Hilbert–Pólya operator constructed for some non-birotational H | +½ would no longer be uniquely birotational; RS2 contribution dilutes |
| Counterexample to RH found (zero off σ = ½) | Whole project ends |

## Comparison: what RS2 adds vs. what's already in the literature

| Component | Source | RS2 contribution |
|---|---|---|
| Functional equation $\xi(s) = \xi(1-s)$ | Riemann 1859 | none — quoted |
| Theta modular weight $\sqrt{t}$ in θ(1/t) | Jacobi 1828 | physical reading |
| Mellin transform construction of ξ | Riemann 1859 | none — quoted |
| Hilbert–Pólya conjecture | Hilbert/Pólya 1910s–50s | none — quoted |
| Berry–Keating operator $H = -i(xd/dx + ½)$ | Berry & Keating 1999 | **physical origin of +½** |
| Metaplectic double cover Mp(2,ℝ) | Weil 1964, Shale 1962 | identified with **birotation** |
| GUE statistics of zeros | Montgomery 1973, Odlyzko 1980s | physical reason (time-reversal breaking) |

The RS2 contribution is concentrated in **two cells of this table**: physical origin of +½, and identification of birotation with the metaplectic cover. Everything else is reorganization and notation.

## What this means for Hard Problems #7

A respectable Paper #7 would:

1. State the physical reading clearly.
2. Acknowledge it is a reading, not a proof.
3. Cite Berry–Keating, Hilbert–Pólya, Montgomery–Odlyzko honestly as the antecedents.
4. Make the LB1 identification (birotation ≡ metaplectic) the **central thesis** of the paper. Defend it on the three structural grounds in §5.2, and at greater rigor than §5.2 reaches.
5. Suggest §7-style falsifiable predictions, marked clearly as speculative.
6. Not claim RH is proved.

Honest viability rating: **B/B+**. Conceptual paper, comparable in scope to Paper #1 (Navier-Stokes qualitative). Not the series climax — that role belongs to Paper #6 (Master Formula, 2.3 ppm) which has quantitative content. Paper #7 is the *philosophical* climax: it explains why the standard model's most stubborn open problem has the structure it has, in RS2 physical terms.

The publication path is harder than for Papers #2–#6 because the contribution is interpretive rather than quantitative. Suitable venues: Foundations of Physics, possibly Annals of Physics with a long honest-assessment section, or as a chapter in a longer monograph rather than a stand-alone Annals submission.

## Recommended next moves

1. **Open prior-art directory** and write `03-prior-art-comparison.md`. Compare against the Jan 13 reasoning chain captured in `PRIOR-ART-NOTES.md`. Where prior-art and cold derivation agree, note the agreement; where they diverge, decide which is right and why.

2. **Numerical sanity check**. Run the .pkl test suites against the cold derivation's predictions in §7. If §7.3 (4n² structure) shows nothing, mark it falsified and remove from the paper draft.

3. **Sharpen LB1**. Either rigorize the identification of birotation with the metaplectic cover (cite or derive from RS2-107, 108, 109 the algebraic data needed) or back off to "structural analogy" language.

4. **Decide on paper format**. Full Annals-style derivation, or shorter Foundations-of-Physics interpretive paper. Joe's call.

The cold derivation is complete as a structural reframing. It is honest. It is not a proof. It supplies the +½ with physical origin in RS2 birotation, which is the entire claim to novelty.
