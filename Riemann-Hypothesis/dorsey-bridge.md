---
title: "Dorsey Bridge — QA, RS2, and the §7 Numerics on a Single Page"
author: Joe Van Horn (with Claude Opus 4.7)
date: 2026-04-28 (late evening)
status: connective tissue — pulls together the four-corner mapping between Dorsey's *Prime Numbers Encode a Wavefield* (2023), the QA Comprehensive paper (Vanhorn 2025), the RS2 cold derivation (`cold-derivation/01-riemann-cold.md`), and the §7 numerical sanity check (`section-7-tests/results.md`). Designed to drop in as a Paper #7 appendix or to stand alone as the orientation document for someone reading the project for the first time.
---

# Dorsey Bridge

## 1. What this document is

The Riemann-zero work in this directory threads four sources that talk past each other in the original literature:

- **Damon Dorsey** (2023), *Prime Numbers Encode a Wavefield* (Zenodo, DOI 10.5281/zenodo.17269878; YouTube channel "The Prime Scalar Field") — geometric framework: prime numbers as waves on a 3D sphere, bounded yet infinite, **non-touching**.
- **Joe Van Horn** (2025), *Qualia Algebra Comprehensive* — cites Dorsey as the visualization of [0,0,0,0] / Potential Space; positions Dorsey as one of five independently-developed frameworks (QA, RS2, Knot, Prime Field, Pendulum) that converge on `n = 3 dimensions + harmonic organization + discrete spectrum + observer-critical` predictions.
- **RS2 cold derivation** (this directory, 2026-04-28) — derives σ = ½ as a metaplectic / birotational half-quantum from the canonical commutator [x, p] = i, identifies the Berry–Keating operator H = −i(x d/dx + ½) as the natural Hilbert–Pólya candidate, eight readings of ½ unified.
- **§7 numerical sanity check** (`section-7-tests/`, 2026-04-28) — first 500 non-trivial zeros, GUE pair correlation 89% better fit than uniform, nearest-neighbor spacing 3.8× better fit than GOE.

This document explains how the four pieces are different faces of one object, and why the convergence is load-bearing rather than coincidental.

## 2. The single fact, four ways

The fact: **σ = ½ is the unique critical-strip value at which all prime-waves share one amplitude scale**, so they can populate a bounded structure without colliding.

| Vantage | How σ = ½ shows up |
|---|---|
| **Number theory** (Riemann–von Mangoldt explicit formula) | Each non-trivial zero ρ = ½ + iγ contributes a wave $x^\rho/\rho = (\sqrt{x}/\rho) \exp(i\gamma \log x)$ to the prime-counting function. The √x amplitude is exactly σ = ½. |
| **Dorsey geometry** (3D-sphere wavefield) | Bounded sphere supports infinite non-touching prime-waves. "Non-touching" requires a single shared amplitude scaling across all waves. The unique σ that supplies this is ½. |
| **RS2 cold derivation** (§5.4 (h), eight readings) | The √x amplitude is the spectral-acoustic face of the canonical commutator [x, p] = i. Same +½ that lifts the dilation generator −i x d/dx to the Berry–Keating operator. |
| **§7 numerical statistics** (GUE level repulsion) | Zeros physically repel: R₂(r) → 0 as r → 0 in the pair correlation, observed at N = 500. "Non-touching" is the empirical signature of GUE level repulsion. |

The four columns are the **same fact**: the critical line is the unique amplitude-locus where the prime-wave manifold is bounded, the spectrum is discrete, and the levels repel.

## 3. The QA bridge — why Dorsey lands at [0,0,0,0]

In QA Comprehensive (lines 481, 569–577), Joe places Dorsey's prime sphere at the **Potential Space [0,0,0,0]** vertex of the four-space structure:

> "Dorsey's Prime Scalar Field: 3D sphere with infinite non-touching prime waves — bounded finite structure with infinite variation."
>
> "The bounded-yet-infinite structure resonates with [0,0,0,0] as: bounded (finite structure / sphere), infinite (infinite variation / non-touching prime waves), pre-manifest (exists prior to observation/collapse), potential (all possibilities present but none actualized)."

In QA's four-space scheme:
- **[0,0,0,0]** Potential Space — pre-manifest, formless void
- **[1,0,0,0]** Witness Space — pure observer
- **[1, x, y, z]_i** Personal Space — individual interface
- **[1,1,1,1]** Consensus Space — physical realm

Dorsey's prime sphere is the **rendering** of the bounded-infinite structure at [0,0,0,0]. It's not directly the physical / Consensus realm; it's the spectral substrate from which physical structure draws. This matches the role of ζ(s) in physics: not a force or field, but a generating function whose zeros encode the harmonic skeleton of arithmetic.

## 4. The RS2 bridge — why Dorsey lands at σ = ½

The cold derivation §5.4 unifies eight readings of ½ under the canonical commutator [x, p] = i. Reading (h) — the √x prime-wave amplitude — is Dorsey's geometric reading translated into the RS2 frame:

> "if even one zero sat off the critical line at σ = ½ ± ε, its wave would carry amplitude $x^{1/2 \pm \epsilon}$ and would either grow or shrink relative to the rest, breaking the bounded wavefield."

This is exactly Dorsey's "non-touching" requirement, stated in number-theoretic terms. The RS2 contribution is to identify *why* this requirement points to a specific operator: the Berry–Keating H = −i(x d/dx + ½), whose half-integer weight matches the metaplectic / birotational structure of the EE-quaternion algebra (RS2-107..109, Peret) and the spin-½ structure of the Cayley–Dickson doubling chain (Nehru, *Quaternion Organon* §8).

## 5. The §7 numerical bridge — empirical content of "non-touching"

Dorsey's "non-touching prime waves" is a geometric statement. The §7 numerical sanity check turns it into an empirical one. From `section-7-tests/results.md`:

- **Pair correlation** R₂(r) = 1 − (sin πr / πr)² (Montgomery's conjecture) fits observed pair correlation **89% better than uniform**. The signature feature: R₂(r) → 0 as r → 0. Levels repel; they do not touch.
- **Nearest-neighbor spacing distribution** fits the GUE Wigner surmise P(s) = (32/π²) s² exp(−4s²/π) **3.8× better than GOE** and **23× better than Poisson**. Mean unfolded spacing 0.9997 (essentially 1, perfect unfolding). The Wigner surmise has a strict zero at s = 0; level repulsion is structural.

The convergence: the geometric "non-touching" requirement of Dorsey's sphere is the **empirical** content of the GUE level-repulsion that random matrix theory predicts for a time-reversal-broken Hermitian operator. Both express the same fact about the zero set, from different angles. The §7 numerical results provide the bridge.

## 6. Why the convergence is load-bearing

Four independent vantages reach σ = ½ for *different* reasons:

1. **Riemann's** original analytic-continuation construction lands σ = ½ as the fixed line of the functional equation s ↔ 1 − s.
2. **Dorsey's** sphere lands σ = ½ as the unique amplitude-locus where prime-waves can share one scaling.
3. **RS2 cold derivation** lands σ = ½ as the metaplectic / birotational / single-loop-inductance / Cayley–Dickson half-step quantum from the canonical commutator [x, p] = i.
4. **§7 numerical statistics** lands σ = ½ as the empirical fact that observed zeros exhibit GUE level repulsion at N = 500.

When four independent passes — analytic, geometric, algebraic, statistical — converge on the same answer for different reasons, the answer graduates from "one framework says so" to "the structure forces this." This is the same epistemic move described in `feedback-cold-rederivation-methodology.md` for the Jan 13 / Apr 28 cold-vs-prior-art convergence: convergence between independent passes is itself a load-bearing result.

## 7. Open extensions

Three threads remain to close, in order of weight:

1. **Hilbert–Pólya construction** (`02-honest-assessment.md` §9 honest gap): build the specific Hilbert space + measure + boundary conditions that make H_RS = xp + px have spectrum exactly equal to {2γ : ζ(½ + iγ) = 0}. Open since 1910s. The RS2 reading + Dorsey rendering + §7 statistics together suggest the Hilbert space wants to be a function space on the unit-quaternion sphere S³ (the EE-quaternion atomic-zone manifold), but the precise measure and boundary conditions remain to be written down.

2. **Higher-N numerical extension** (§7 follow-up): extend zero count from 500 to 10⁵ using Odlyzko's precomputed tables. GUE convergence tightens visibly with N. Form factor K(τ) = |Σ exp(2πi γ_n τ)|² / N is a different statistical lens that should also match RMT.

3. **Cross-paper coupling** (§7.3 4n² prediction): the speculative test failed to confirm structure under naive binning. A more principled binning rule from RS2-105 (quantum π = 4) and the cold derivation Hard Problem #5 might still find the predicted 4n² envelope. Currently flagged speculative.

## 8. References

- Dorsey, D. (2023). *Prime Numbers Encode a Wavefield* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17269878
- Dorsey, D. (2025). *FINALLY! The hidden patterns in PRIME NUMBERS show themselves — we were only looking in one dimension*. YouTube channel "The Prime Scalar Field". https://www.youtube.com/watch?v=Y9f-Gq42Pxg
- Vanhorn, J. (2025). *Qualia Algebra: Comprehensive*. Zenodo. (See lines 481, 569–577 for Dorsey integration with [0,0,0,0] Potential Space.)
- Vanhorn, J. (2026). *Riemann Hypothesis — Cold Re-Derivation from RS2 First Principles*. `cold-derivation/01-riemann-cold.md`. Eight readings of ½ unified, §5.4.
- Vanhorn, J. (2026). *Section 7 Numerical Sanity Check Results*. `section-7-tests/results.md`. First 500 zeros, GUE pair correlation + nearest-neighbor confirmed.
- Peret, B. (2011–2016). *Reciprocal System v2 Papers* (RS2-101..109). Distillation: `~/RS-Framework-Bridge/RS2-foundations/00-FIRST-PRINCIPLES.md`.
- Nehru, K. V. K. (1997). *Some Thoughts on Spin*. In *Quaternion Organon* (transcribed: `~/RS-Framework-Bridge/RS2-foundations/nehru-thoughts-on-spin.md`).
- Berry, M. V., & Keating, J. P. (1999). The Riemann zeros and eigenvalue asymptotics. *SIAM Review*, 41(2), 236–266.
- Montgomery, H. L. (1973). The pair correlation of zeros of the zeta function. *Analytic Number Theory*, AMS Proc. Symp. Pure Math. 24, 181–193.
- Odlyzko, A. M. (1987–). Numerical verification of Riemann zeros and pair correlation, multiple papers.
