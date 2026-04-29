---
title: "The Inter-Regional Ratio — Nehru 1985"
author_original: Prof. K.V.K. Nehru, Ph.D.
source_publication: Reciprocity Vol. 14 No. 2, 1985 (ISUS)
recovered_from: Google Drive (file ID 19cN7PuPZ…) — folder "Nehru papers"
distillation_by: Claude Opus 4.7 (cold read 2026-04-28 late evening)
purpose: Quantify cross-region transmission factors for the Yang-Mills cold derivation
---

# The Inter-Regional Ratio — Nehru 1985

The Inter-Regional Ratio (IRR) is Larson's quantitative factor for transmission of motion across region boundaries. Nehru clarifies the calculation by reframing "orientation" → "degrees of freedom" / "intrinsic existential possibilities."

This paper is the **quantitative companion** to the atomic vs nuclear zone partition we had only sketched in Peret's EE → RS dictionary (IMG_5298). Nehru gives the actual numbers.

---

## §1. The basic counting

For a unit of motion to be transmitted across a region boundary, the probability is **1/f**, where f = total possibilities for the motion within the source region.

### 1.1 Rotational degrees of freedom in 3D vector space

A 1D rotation in 3D vector (extension) space has 2 directions (+1, -1) per dimension, for 3 dimensions:

$$f_{1D\text{-rot in 3D-vec}} = 2 \times 2 \times 2 = 8$$

A 2D rotation is a coupled pair of 1D rotations. The coupling causes degeneracy — reversing both component rotations leaves the 2D rotation unchanged: (+1)(+1) = (-1)(-1) = +1. So the 8 → 4 reduction:

$$f_{2D\text{-rot in 3D-vec}} = 8 / 2 = 4$$

Generally: `f = p^d / n` for a d-dimensional motion in n-dimensional vector space with p possibilities per dimension.

### 1.2 Scalar dimensions — independent

Scalar (space-time) dimensions, unlike vector dimensions, are **independent**. So degrees of freedom *add* rather than multiply:

$$f_{\text{scalar 3D}} = 2 + 2 + 2 = 6 \text{ (per Larson, NFoS p. 84)}$$

### 1.3 Atomic structure

The atom comprises **two 2D magnetic rotations + one 1D electric rotation**. Each rotation lives in 3D vector space:

$$f_{\text{atom}} = (2^2/3 \text{-style factor})(2^2/3)(2^1/3\text{-electric-1D}) \cdot (\text{3D vector space})$$

Working it out per Nehru's formula (8 for 1D-rot, 4 for 2D-rot in 3D-vec):

$$f_{\text{atom rotation}} = 4 \times 4 \times 8 = \boxed{128}$$

This is the **total rotational degrees of freedom** for the atomic structure in 3D vector time.

---

## §2. Vibrational degrees of freedom

The photon basic vibration adds vibrational degrees of freedom on top of the rotational base.

A 1D vibration has only 1 possibility per dimension (forward and backward = one oscillation), so $1^3 = 1$ per scalar dimension. But there are 3 scalar dimensions and the vibration occupies only one: 3 free choices, so $f_{\text{1D-vib in 3 scalar dims}} = 3$.

For the rotational unit, this contributes 1/3 added degree of freedom.

For 2D rotation founded on 1D vibration: $p = 3, n = 2$ → $3^2 = 9$ vibrational possibilities. So one rotational system gets 1/9th additional degree of freedom from vibration.

Atomic structure has **two** 2D rotational systems, so vibrational contribution is **2/9th**:

$$f_{\text{atom total}} = 128 + 128 \times \frac{2}{9} = 128 + 28.444 = \boxed{156.4444}$$

Subatom (one 2D rotational system only):

$$f_{\text{subatom total}} = 128 + 128 \times \frac{1}{9} = 128 + 14.222 = \boxed{142.2222}$$

These are the quantitative IRR values.

---

## §3. Implications for the Yang-Mills cold derivation

### 3.1 The "128" factor as the atomic-zone rotational quantum

128 = 2 × 4 × 4 × 8 / (degeneracy factors) = the total existential possibilities of an atomic structure within 3D vector time.

In RS spherical-harmonic decomposition (Joe's `larson_128_report.txt` from Drive), the factor 128 emerges at quantization ratio 1.0000 with error 0.0000 across multiple (ℓ, m) modes — confirming 128 as a fundamental atomic-zone rotational quantum.

For YM mass-gap purposes: the YM excitation lives in the atomic zone (4D quaternion), so the relevant degree-of-freedom count is 128 + structure.

### 3.2 Pure-YM (no electric) IRR adjustment

A pure-YM excitation has only the magnetic-2D rotational content — no electric-1D. Removing the 1D electric rotation factor (8 in the chain 4 × 4 × 8 = 128) gives:

$$f_{\text{pure-magnetic}} = 4 \times 4 = 16$$

Plus 2/9 vibrational = 16 + 16 × 2/9 = 19.555...

Or 4 × 4 × 1 (electric-1D contributes 1 in pure magnetic case, since only the trivial direction remains) = 16, IRR 16 × 11/9 = 19.555...

This is the *transmission factor* for a pure-magnetic-2D YM excitation across the atomic-zone boundary. It does NOT directly give the mass; the mass is set by the t³/s³ inductance quantum (= primary mass p, calibrated to 1 amu). But it constrains the *eigenmode density* of YM excitations within the atomic zone.

### 3.3 Higher glueballs and the IRR ladder

The lightest glueball is the lowest-frequency mode of a pure-magnetic-2D excitation. The next glueballs (2++ tensor at ~2.3 GeV, 0-+ pseudoscalar at ~2.5 GeV in lattice predictions) might map onto the next available rotational eigenmodes within the atomic-zone constraints.

If `m_n / m_amu = ln(2π × n)` for n = 1, 2, 3, ..., we'd predict:
- n = 1: 1.711 GeV (matches lightest scalar glueball, f₀(1710))
- n = 2: ln(4π) × m_amu = 2.531 × 931.494 MeV = 2.358 GeV (close to 2++ tensor glueball)
- n = 3: ln(6π) × m_amu = 2.937 × 931.494 MeV = 2.736 GeV (close to 0-+ pseudoscalar)

This is a **falsifiable prediction**: the YM excitation spectrum should follow `ln(2πn) × m_amu` for low n, with the IRR factor 128 (or 16 for pure-magnetic) limiting the mode count.

(Speculative; needs careful analysis. Lattice predictions for higher glueballs have larger error bars.)

---

## §4. Open questions

1. Why 1D-rotation in scalar dimensions counts as 6 (additive) but 1D-rotation in vector dimensions counts as 8 (multiplicative)? Nehru cites Larson NFoS p. 84 but the asymmetry between scalar and vector counting needs sharpening.
2. The 1/9 vibrational factor depends on a 2D-rotation founded on 1D-vibration with 3 scalar dimension choices — but why exactly 2/9 for atomic (two rotational systems) vs 1/9 for subatomic (one)?
3. Is the "16 × 11/9 = 19.55..." pure-magnetic IRR a real physical quantity, or does removing the electric component require a separate calculation?
