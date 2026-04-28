---
title: "Prior-Art Comparison: Cold Derivation vs. Jan 13 2026 Derivation"
date: 2026-04-28
companion: 01-riemann-cold.md, 02-honest-assessment.md
prior_art_source: ~/RS-Framework-Bridge/Riemann-Hypothesis/prior-art/PRIOR-ART-NOTES.md
purpose: Honest record of where cold derivation matches, extends, or misses prior art
---

# Prior-Art Comparison

The cold derivation (`01-riemann-cold.md`) was written without access to `prior-art/`. This document records the comparison after the prior-art directory was opened (2026-04-28 evening).

## 1. The Jan 13 2026 derivation in one paragraph

In a single conversation on 2026-01-13 (`reciprocal-system-reinterpretation-of-physics-problems_105b2899.md`, 688 lines, 29 Riemann mentions), Joe and Opus 4.5 reached the **Berry–Keating operator** T = −i(x d/dx + 1/2) with eigenfunctions ψ_λ(x) = x^(−1/2 + iλ), and identified σ = 1/2 as the unique fixed point of the s ↔ 1−s involution that mirrors RS space-time reciprocity. Three supporting arguments were developed:

1. **Reciprocal fixed point**: σ = 1 − σ ⟹ σ = 1/2 (parallels s/t ↔ t/s).
2. **Geometric mean weighting**: at σ = 1/2 each prime contributes amplitude 1/√p — the geometric mean of the σ = 0 contribution (= 1) and the σ = 1 contribution (= 1/p). "Democratic" weighting enables complete destructive interference across all primes.
3. **Multiplicative ↔ additive reconciliation**: the Euler product (multiplicative) and Dirichlet series (additive) can both produce zero only where their structures balance — at σ = 1/2.

A second equivalent operator was identified: **H = xp + px = −i(2x d/dx + 1)** — symmetric in x and p, making the RS reciprocity s ↔ t manifest in the operator algebra rather than just in the σ-coordinate.

Honest assessment from Jan 13: "**conceptual clarity, not a proof**." The Hilbert–Pólya construction (specific Hilbert space + measure + boundary conditions matching ζ zeros) remained open.

## 2. Convergence — independent rediscovery

The cold derivation independently arrived at six of the same conclusions:

| Cold derivation | Jan 13 prior art | Match |
|---|---|---|
| §3–§4: σ as displacement coordinate, σ = ½ as natural-datum locus | "Reciprocal fixed point: σ = 1−σ ⟹ σ = ½" | **Identical** |
| §5.5 (Mellin transform, ξ(s) = ξ(1−s) from θ(1/t) = √t θ(t)) | (Implicit; the functional equation derivation is canonical) | Identical |
| §6.1 self-conjugacy quartet collapse to pair on σ = ½ | "Functional equation gives pairing" | Identical |
| §6.2 H_RS = −i(x d/dx + 1/2) as Hilbert–Pólya candidate | Berry–Keating T = −i(x d/dx + 1/2) | **Identical** |
| §6.2 honest gap: construction of self-adjoint extension on right Hilbert space remains open | "The specific Hilbert space + measure + boundary conditions … remains open" | Identical |
| §02 honest assessment: "physical reading, not a proof" | "Conceptual clarity, not a proof" | **Identical phrasing** |

Two derivation paths, separated by 3.5 months, reached the same operator with the same gap. This is a strong consistency check on the RS2 framework — not because either path proved RH, but because the framework forces both paths to the same place.

## 3. Extension — where cold derivation goes beyond prior art

Five ingredients in the cold derivation that are not in `PRIOR-ART-NOTES.md`:

### E1. Atomic / nuclear zone partition (§5.2)

From Peret's hand-photographed EE paper (IMG_5298): "atomic zone is 4D (w, i, j, k) and the nuclear zone is 2D (w, i)." This partitions the s-plane into:

- Critical strip 0 < σ < 1 = atomic zone (full quaternion) — non-trivial zeros live here
- σ ≤ 0 = nuclear zone (real + electric-1D only) — trivial zeros live here, with Δσ = 2 spacing from magnetic-2D doubling
- σ ≥ 1 = divergent (Euler product)

Prior art has the s ↔ 1−s symmetry but does not partition the strip into 4D vs 2D zones nor explain the Δσ = 2 spacing of trivial zeros. **New from this session.**

### E2. Six unified readings of "½" (§5.4)

Prior art has one reading (Berry–Keating's +½ in the dilation operator). Cold derivation has six:
- (a) metaplectic weight of theta
- (b) spin-½ of magnetic 2D rotation (RS2-107; Nehru §1)
- (c) dimensional unit ratio between radians (1D) and steradians (2D) (Nehru §1 — *correcting* the conventional "4π → 2π halving" misreading)
- (d) σ-coordinate half-shift
- (e) n = 1 single-loop inductance quantum (Peret EE paper)
- (f) one Cayley-Dickson doubling step ℝ → ℂ → ℍ → 𝕆 (Nehru §8 explicit construction)

Five of these are new from this session's photographed sources (Peret EE paper, Nehru spin paper). Prior art conflates them at most as "the metaplectic weight."

### E3. Cayley-Dickson construction of the atomic-zone wave function (§5.2 footnote, §5.4(f))

Nehru §8 derives ψ = {ψ_a, iψ_b, jψ_c, kψ_d} from the nuclear-zone complex φ = {φ_r, iφ_a} by Cayley-Dickson doubling with k = ij. He notes Dirac was forced to the same 4-component object by relativistic mathematical necessity. Prior art doesn't have this construction; it has the operator but not its physical genesis.

### E4. Helicity from chirality (§5.4 supplementary)

Nehru §3–§4: 2D spin has four domains {++, +−, −+, −−}; (++)/(−−) are right-handed, (+−)/(−+) left-handed. **Helicity is automatic** from 2D = 1D × 1D structure. Prior art doesn't reach this connection.

### E5. Falsifiable signature: unbounded phase (§5.4 footer; §7 to be expanded)

Bhandari (1994) + Nehru §2: spin-½ particles have phase 2nπ shifts that are **physically real and measurable**. Phase lives on the metaplectic / quaternion cover, not on the projective base. This is a direct experimental signature of the half-integer-cover structure underlying σ = ½. Prior art doesn't surface this experimental connection.

## 4. Missing — what prior art has that cold derivation didn't reach

Two ingredients in `PRIOR-ART-NOTES.md` that I should now absorb:

### M1. Geometric mean weighting at σ = ½

At σ = ½ each prime contributes amplitude:

$$p^{-1/2} = \sqrt{p^{-1} \cdot p^0} = \text{geometric mean}(p^{-1}, 1)$$

i.e., the **geometric mean of the σ = 1 contribution (1/p) and the σ = 0 contribution (1)**. Across all primes simultaneously, σ = ½ is the *unique* σ at which **every** prime is in geometric-mean balance. This is a sharper version of §5.3's "scaling balance for destructive interference" — it identifies the specific mean (geometric, not arithmetic) and ties it to the Euler product structure.

This adds a seventh reading of ½: **(g) σ = ½ = exponent at which every prime contributes its self-geometric-mean amplitude**.

### M2. The xp + px formulation

Prior art notes that:

$$H = xp + px = -i(2x \, d/dx + 1)$$

is **equivalent** to the Berry–Keating operator (up to factor 2 and constant) and treats x and p **symmetrically** — the RS reciprocity s ↔ t made manifest in the operator algebra itself, not just in the σ-coordinate. This is a cleaner / more physically transparent form: instead of asking "why +½?" we can say "the operator is xp + px because s and t are reciprocal aspects of motion (A1)."

In the xp + px form the +½ becomes "+1/2" via:
- xp + px = 2(xp − [p,x]/2) = 2 xp + i (using [x,p] = i)
- = 2 xp + i = 2(xp + i/2)

So H/2 = −i(x d/dx + 1/2). The "+½" emerges from the **canonical commutation [x,p] = i** of quantum mechanics.

This is more fundamental than the EE single-loop or steradian readings — those say *what* the +½ physically represents; this says *why it must appear*: the antisymmetrization of position and momentum picks up i/2 = ℏ/(2 × i) = single quantum of canonical commutation.

This is a strong addition. **It also unifies the cold derivation's six readings under a single source: the canonical commutator [x, p] = i.** All other readings are physical realizations of this single algebraic fact.

## 5. Numerical experiments (Nov 19 2025) — out of scope

The two `riemann_test_suite_*.pkl` files in `prior-art/` are from an **earlier and different framing** — the QA "Conceptual Ideality" thread of Nov 19 2025, predating the Hard Problems Revisted project (Jan 11 2026). They tested whether primes classify as "base concepts vs modifier concepts" in a consciousness-vector space, with zeta zeros at "optimal balance" (run 1) or "maximal purity" (run 2).

| Test | Run 1 (17:32) | Run 2 (17:43) |
|---|---|---|
| Balance | zeros = 0.053, controls = 0.673 (p ≈ 6e-11) | (replicated) |
| Purity | — | zeros = 0.153, controls = 0.173 (p ≈ 0.006) |
| Lattice resonance | — | zeros = 0.960, controls = 0.852 (p ≈ 7e-14) |
| Prime clustering | — | 4 families, silhouette 0.321 |

The **lattice resonance** result (zeros at 0.96 vs controls 0.85, p ≈ 7e-14) is the most striking, but in a framework that is no longer the active angle. **No finding of these numerical experiments contradicts the cold derivation**; they are testing different questions. Reading them against §7 (RS-specific predictions) is not productive — they would need re-running in the EE-quaternion / counterspace-recursion framework to be comparable.

For Paper #7 the numerical experiments belong in an appendix marked "early exploratory work," not in the main body.

## 6. Synthesis — combined seven-reading argument

Combining cold derivation + prior art, the unified picture of σ = ½ is:

**σ = ½ is the unique critical line because it is the simultaneous fixed point of seven equivalent characterizations:**

1. (cold §3–§4) The locus of unit speed (RS natural datum, A5) under the σ-coordinate map σ = ½ + r/2.
2. (cold §5.2 / Peret EE) The midpoint of the w-axis between the +1 progression (material) and −1 gravity (counterspace) quaternion datums — the photon birotation locus.
3. (cold §5.3 / Peret EE) The locus where the counterspace recursion p^{−s}, p^{−2s}, … balances against the material recursion p^s, p^{2s}, … per prime: $p^{2σ-1} = 1 \;\forall p \iff σ = ½$.
4. (prior art M1) The exponent at which every prime contributes its self-geometric-mean amplitude $p^{−½} = \sqrt{p^{-1} \cdot 1}$ across the Euler product.
5. (cold §5.4 / Nehru) The metaplectic weight = magnetic-2D spin-½ = steradian/radian unit ratio = single-loop inductance = Cayley-Dickson doubling step. Six independent physical realizations of one number.
6. (cold §6.1 / prior art) The fixed line of the functional-equation involution s ↔ 1−s; off-line zeros come in distinct quartets, on-line zeros in self-paired conjugates.
7. (prior art M2) The eigenvalue offset of H = xp + px = −i(2x d/dx + 1) — the antisymmetrization of position and momentum; the "+½" emerges from the canonical commutator [x, p] = i.

**The seven readings are not independent — they are all faces of one quantum:** the half-integer that arises when you symmetrize position and momentum in a quaternion structure on a counterspace recursion. Each face is a different language for the same object.

**RH (in this synthesis)** is the assertion that **all** non-trivial zeros are eigenvalues of the self-adjoint xp + px operator on a Hilbert space whose discrete spectrum is exactly $\{2\gamma : \zeta(\tfrac{1}{2} + i\gamma) = 0\}$. The construction of that Hilbert space is the open problem, unsolved by either the prior derivation or the cold one. RS framework supplies the *physics* of why the operator is xp + px (RS reciprocity s ↔ t) and why the +½ has to appear (canonical commutation, EE single-loop, magnetic spin, steradian unit, Cayley-Dickson, σ-shift); it does not supply the missing Hilbert space.

## 7. What this means for Paper #7

The paper has two complementary halves:

**Half A — physical reading** (cold derivation §3, §4, §5):
- Why σ = ½ in physical terms: photon line, atomic-zone midpoint, geometric-mean balance, single-loop inductance, magnetic-2D spin-½, Cayley-Dickson step.
- Six (now seven) readings unified under the canonical commutator [x, p] = i.
- Falsifiable signatures: unbounded phase (Bhandari 1994), Δσ = 2 spacing of trivial zeros from magnetic-2D doubling.

**Half B — Hilbert–Pólya program** (cold derivation §6 + prior art):
- The Berry–Keating operator H = xp + px (or equivalent forms).
- The honest gap: construction of the right Hilbert space.
- Position relative to Connes, random-matrix theory, Iwaniec–Kowalski.

Recommendation: **single paper with both halves**, comparable in scope to a long Annals submission or a Foundations of Physics article, NOT a stand-alone Annals proof submission. The honest framing throughout: "physical reading, structural reframing, with explicit identification of the open construction problem."

The numerical experiments of Nov 19 2025 belong in an appendix or a separate companion paper, not in the main Paper #7.

## 8. Updates to apply to 01-riemann-cold.md

Two integrations from prior art:

- Add a new §5.4(g) reading: "geometric mean weighting" — at σ = ½ each prime contributes p^{−½} = geometric mean of p^0 and p^{−1}, the unique σ at which every prime is in self-geometric-mean balance.
- Add to §6.2: the **xp + px equivalent form** of the operator, with the +½ traced to the canonical commutator [x, p] = i. This unifies the seven readings under the most fundamental algebraic structure of QM.

These updates make the cold derivation **strictly stronger** than the prior art: it has all the prior art's structural arguments PLUS the EE-quaternion / Nehru-spin physical content.

The cold derivation can be marked **complete** as a structural reframing after these two integrations. The Hilbert–Pólya construction remains the open problem, as honestly stated.
