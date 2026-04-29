---
title: "FSC Cold Derivation — Pre-Seal Honest Assessment"
problem: HP#4 (Fine Structure Constant)
status: pre-seal — written before opening `../prior-art/`
date: 2026-04-29
companion: 01-fsc-cold.md
purpose: lock the cold-pass epistemic state in writing before any contact with the Jan 13/14 derivations, so the prior-art comparison in 03 can detect drift in either direction
---

# Honest Assessment — FSC Cold Pass

This document records, *before* opening `../prior-art/`, what the cold pass in `01-fsc-cold.md` actually established. Any drift between this assessment and the prior-art comparison in `03-prior-art-comparison.md` is a signal — either the cold pass overclaimed (and the prior pass had it more honestly), or the prior pass overclaimed (and the cold pass has caught drift, as in HP#3 Hierarchy).

## Taxonomy

The cold pass produces claims at four epistemic levels. Each claim in `01-fsc-cold.md` belongs to exactly one level.

### Level 1 — Numerically verified

These are direct calculations whose verification requires no RS2 framework, only standard mathematics and a calculator.

- **Claim**: 4π³ + π² + π = 137.036303776... (to 9 decimals)
  - **Status**: numerically verified
  - **Where in §**: §1.1, §6
  - **What's at risk**: nothing — this is arithmetic

- **Claim**: This value differs from CODATA 1/α = 137.035999084 by 0.000305, a relative 2.22 × 10⁻⁶ (2.22 ppm)
  - **Status**: numerically verified
  - **Where in §**: §1.2, §6
  - **What's at risk**: nothing — arithmetic against published CODATA

- **Claim**: V(S¹) = 2π, V(S²) = 4π, V(S³) = 2π², V(S⁵) = π³ via standard sphere-volume formulas
  - **Status**: numerically verified
  - **Where in §**: §3.1, §3.2
  - **What's at risk**: nothing — standard mathematics

- **Claim**: 4π³ = V(S³) · V(S¹), π² = V(S³)/2, π = V(S¹)/2
  - **Status**: numerically verified (algebraic identity given the volumes above)
  - **Where in §**: §3.3, §4
  - **What's at risk**: nothing — algebraic identity

### Level 2 — Structurally established within RS2

These are identifications that follow from canonical RS2 axioms (LB1-LB6 in §2). They depend on the framework being correct as published; if RS2 is wrong about, e.g., the atomic zone being 4D quaternion, these claims fall.

- **Claim**: The atomic-zone unit-rotation manifold is S³, the unit sphere in ℍ ≅ ℝ⁴
  - **Status**: structural — follows from LB2 (atomic zone = quaternion ℍ) + standard topology (unit sphere in ℝ⁴ is S³)
  - **Where in §**: §3.2
  - **Sources**: Peret EE → RS dictionary §5; Nehru "Thoughts on Spin" §8
  - **What's at risk**: if Peret's atomic-zone-as-quaternion claim is incorrect, this falls

- **Claim**: The nuclear-zone unit-phase manifold is S¹, the unit circle in ℂ ≅ ℝ²
  - **Status**: structural — follows from LB2 (nuclear zone = complex ℂ) + standard topology
  - **Where in §**: §3.2
  - **Same sources**: Peret EE §5; Nehru §7
  - **What's at risk**: same as above

- **Claim**: The factor (2π rad)/(4π sr) = ½ is the spin-½ unit conversion (rad → sr), not a halving of any single quantity
  - **Status**: structural — directly quoted from Nehru §1 (LB3)
  - **Where in §**: §2 (LB3), §3.1, §4.2 reading (a)
  - **Source**: Nehru "Some Thoughts on Spin," *Reciprocity* XXVI/3 §1
  - **What's at risk**: if the canonical Nehru reading of spin-½ is wrong, this falls; otherwise robust

- **Claim**: V(S²) = 4π is the steradian coverage that appears in F = q₁q₂/(4πε₀r²)
  - **Status**: structural — standard physics identification
  - **Where in §**: §3.1 (V(S²)), §6 (consistency)
  - **What's at risk**: nothing in RS2 specifically — this is canonical Coulomb law

### Level 3 — Plausible structural reading, post-hoc fit

These claims are the heart of the cold pass and where epistemic care matters most. Each is a *reading* of a feature of the closed form that fits within RS2 but is not derived from RS2 axioms alone — the closed form was visible above the seal, and the cold pass identified the volume content of each term post-hoc.

- **Claim**: The leading term 4π³ is the joint phase volume of atomic×nuclear-zone EM coupling, V(S³) × V(S¹)
  - **Status**: post-hoc structural reading
  - **Where in §**: §3.3, §4.1
  - **What's derived**: V(S³)×V(S¹) = 4π³ (algebra)
  - **What's NOT derived**: that the EM coupling event's phase-space measure is the *product* of atomic and nuclear zone volumes (versus, e.g., a sum, integral over a fiber bundle, or some weighted combination). This is asserted by analogy to product-space measure, with the volumes matching as the verification.
  - **What's at risk**: a different RS2-natural decomposition (e.g., V(SU(2)) acting on V(S²)/2 via projective fibration) might also give 4π³ but not via this product
  - **Honesty check**: this claim is *suggested* by RS2 + the closed form; it is not *forced* by RS2 alone

- **Claim**: The middle term π² is the half-volume V(S³)/2, the atomic-zone single-side contribution
  - **Status**: post-hoc structural reading
  - **Where in §**: §4.2
  - **What's derived**: V(S³)/2 = π² (algebra)
  - **What's NOT derived**: that the "single-side without partner" contribution should be exactly half the manifold volume (versus, e.g., the full volume, or the volume divided by the four-domain projection ratio 4 → 2 = ½ which gives the same answer here but for different reasons)
  - **Honesty check**: two distinct readings (hemisphere; four-domain projection) both give the right number ½. Convergence on the half-coefficient is suggestive, but the underlying mechanism that *forces* this term to be half-volume is not derived.

- **Claim**: The trailing term π is the half-volume V(S¹)/2, the nuclear-zone single-side contribution
  - **Status**: post-hoc structural reading
  - **Where in §**: §4.3
  - **What's derived**: V(S¹)/2 = π (algebra)
  - **What's NOT derived**: same caveat as the middle term — the half-cycle reading is suggestive but not forced
  - **Honesty check**: by far the simplest of the three terms, but inherits the same "why halved" question

- **Claim**: The form 1/α = T₁ + T₂ + T₃ exhausts the channels
  - **Status**: post-hoc structural reading, weakest claim in this level
  - **Where in §**: §4.4
  - **What's NOT derived**: that there are exactly three channels (joint, atomic-only, nuclear-only) and not, e.g., a "neither" channel contributing 1, or higher cross-terms beyond product/half. The closed form has three terms; the cold pass argues the natural classification gives three; consistency, not derivation.

### Level 4 — Open / unresolved / speculative

These are flagged in `01-fsc-cold.md` but not claimed.

- **Open**: the 2.22 ppm residual against CODATA. Four candidate origins listed in §6 (higher-order corrections, quantum-π/analog-π transition, cross-talk between terms, closed form is approximate). None endorsed.
- **Open**: the "/2" coefficients' first-principles origin. Two readings given (hemisphere, four-domain projection); cold pass does not commit to which is operative or whether they're the same thing in different language.
- **Open**: the connection of the +1 unit in the factored form π(4π² + π + 1) to the natural-progression-rate contribution. Suggested in §9 as a cross-check forward to HP#5 (sin²θ_W = π/(4π+1) has the same +1 closure), not derived.
- **Speculative**: Bhandari 2nπ phase memory as a falsifiable signature of the SU(2) cover. Flagged S5 in §7 with "much weaker prediction" caveat.

## What this cold pass does NOT claim

Listed in negative form to be explicit:

(a) That the cold pass has *derived* α from RS2 first principles. It has not. It has identified RS2-natural manifolds whose volumes match the three terms of a closed form known to within 2 ppm of CODATA.

(b) That the closed form 1/α = 4π³ + π² + π is exact. It is not — there is a 2.22 ppm residual that the cold pass does not explain.

(c) That the SU(2)×U(1) decomposition reaches the Standard Model electroweak gauge group structure. It does not — the SU(2)×U(1) here is RS2-native (atomic zone × nuclear zone), not the SU(2)_W × U(1)_Y of electroweak before symmetry breaking. The numerical agreement is suggestive of a deep connection, but the cold pass does not claim to derive electroweak unification.

(d) That the "halving" coefficients are first-principles. They follow from the rad/sr conversion (LB3) and/or four-domain projection (LB2 + Nehru §3) at the structural level, but the *necessity* of half-volume for single-side contributions is not derived.

(e) That the cold pass has predicted α. The closed form's value was visible in `MEMORY.md` above the seal; only the *derivation chain* / *structural reading* was reproduced from sealed conditions. The numerical claim (137.0363 vs CODATA 137.0360) is a verification of the closed form, not a prediction.

## Comparison to prior cold passes (within methodology)

This cold pass sits in the structural-reading class:

- **Riemann (Phase 5.1)**: structural reading of σ = ½ as fixed point of the canonical commutator [x,p]=i; identified Berry-Keating operator H = xp + px. Similar epistemic standing — structural, with numerical verification at N = 500 zeros (§7 tests). FSC analog: this cold pass identified the volume content; numerical verification is the 2.2 ppm match.
- **Yang-Mills (Phase 5.2)**: numerical-prediction class — Δ = ln(2π) × p ≈ 1711 MeV at 1.1% from quenched lattice. Stronger than FSC because the prediction is closed-form and the comparison is within experimental + theoretical uncertainty.
- **Hierarchy (Phase 5.3)**: structural reading + drift detection. Identified F_EM/F_grav = α · (M_Pl/m_p)² decomposition + RS Identity III. Detected prior-pass overclaim that "(M_Pl/m_p) is derived from RS geometry." Similar standing to FSC: structural reading, with one numerical figure (1.5% / 0.06% calibrated) and an honest open gap.

FSC's epistemic standing is closest to Hierarchy: structural reading, numerical match within a known-limited precision, honest open gap on the residual. Both produce volume- or scale-based readings of well-known physical constants.

The 2.22 ppm precision is *higher* than Hierarchy's 1.5% because the FSC closed form is more precise — but the cold-pass *contribution* (structural reading) is no stronger than Hierarchy's. What the cold pass adds is the geometric content (volumes on RS2-natural spheres), not new precision.

## Pre-commitment to drift detection

Three things to watch for when opening `../prior-art/`:

(1) **Different volume reading?** The cold pass settled on V(S³)·V(S¹) for the leading term. If the prior pass derived 4π³ via a different geometric path (e.g., as 8π³/2 = V(S²)·V(S³)/4, or via a Cayley-Dickson-doubling chain ratio), that's a structural fork, and the comparison should note it.

(2) **Tighter half-coefficient derivation?** The cold pass left the /2 and /4 origin as a post-hoc structural reading with two non-uniquely-distinguished candidate origins. If the prior pass derived these from a single forced argument (e.g., from the Nehru rad/sr conversion alone, or from the four-domain projection alone), that would be a structurally tighter reading worth absorbing.

(3) **Prior-pass overclaim?** As in Hierarchy. If the prior pass claimed to *derive* α from RS2 first principles (and not just produce a structural reading), the cold pass should flag this as drift. The honesty-check is: did the prior pass commit to an end-to-end derivation chain whose steps are all level-1 or level-2 claims, or did it (like the cold pass) rest on level-3 post-hoc structural readings? The taxonomy in this document gives the test.

## What would change this assessment

If, after opening `../prior-art/`, the cold pass finds:

- A prior derivation that *forces* the sum-of-three structure from RS2 axioms alone → upgrade the §4.4 claim to Level 2.
- A prior derivation that closes the 2 ppm gap with a clean RS2 mechanism → upgrade Level 4 "open" item to Level 2 or 3.
- A different volume reading that's stronger → demote the cold-pass V(S³)·V(S¹) reading to "one of N candidate readings."
- A prior overclaim on derivation status → upgrade the drift-detection note in `03-prior-art-comparison.md`.

This document is committed in writing before any of those tests. It is the cold-pass anchor.

— end of pre-seal assessment —
