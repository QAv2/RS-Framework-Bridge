---
title: "RS2 First Principles — synthesis of foundation papers 101–109"
author: distillation by Claude (cold read of Bruce Peret's RS2-101..109)
read_date: 2026-04-28
purpose: Foundation document for cold re-derivation of the Riemann Hypothesis from RS2 first principles.
---

# RS2 First Principles

This document synthesizes Peret's RS2-101 through RS2-109 into a single first-principles working set. Per-paper distillations live alongside this file (`RS2-101-creating-a-theory.md` ... `RS2-109-dimensional-thinking.md`).

---

## 1. Minimal axiom set

The minimal axiom set for RS2, drawn from across the nine papers:

| # | Axiom | Source paper |
|---|---|---|
| A1 | The universe is composed of one component, *motion*, in three dimensions, in discrete units, with two reciprocal aspects, *space* and *time*. | RS2-103 (= RS2 Postulate I, modified from Larson's 1979 Postulate I) |
| A2 | The universe conforms to the relations of ordinary mathematics, its primary magnitudes are absolute, and its geometry is *projective*. | RS2-103 (= RS2 Postulate II; Euclidean → projective; commutative dropped) |
| A3 | Minimum scalar magnitude is unity (= 1). There is no zero, no negative, no fractional in the natural framework. Maximum is unlimited but finite. | RS2-104 (with prep in 101) |
| A4 | Scalar motion is the projectively invariant *cross-ratio* of two scalar orientations, with the specific aspects being space and time. Two forms: velocity s/t and energy t/s. | RS2-104 |
| A5 | The natural datum of measurement is *unit speed* (the speed of light in conventional units). All motion is measured as displacement from unity. | RS2-104, 106 |
| A6 | Each scalar dimension has a tristate of motion: unity (default progression), speed (s/t < 1, material), energy (t/s < 1, equivalently s/t > 1, cosmic). | RS2-106 |
| A7 | The number of stable scalar dimensions is uniquely 3, by Nehru's binary-combination argument: n(n−1)/2 = n has solutions n ∈ {0, 3, ∞}. | RS2-106 |
| A8 | Both linear (yang) translation and angular (yin) rotation are *primary* motions. Vibration is not primary; it is shear strain between counter-directed motions. | RS2-107, 109 |
| A9 | Direction can change only at unit boundaries — the discrete-unit "links of a chain" can flex only at junctures. | RS2-105 |

Notes:
- A1 and A2 are the *postulates* proper. A3–A9 are explicit consequences or supplementary commitments stated as foundational across the papers (they are not "derived" within RS2 from A1 and A2 alone — they are independent commitments).
- Larson's pre-1959 third postulate (scalar nature of space-time) was dropped as "necessary consequence" but the proof is not given; it remains a definitional commitment.
- A7 is presented as a derivation, but Nehru's algebra `n(n-1)/2 = n` is offered as a stability heuristic, not a rigorous theorem.

---

## 2. Construction chain

The dependency graph from primitives to the most complex objects (each construction depends on those above it):

```
A1, A2 (motion as one component, projective geometry)
   │
   ├── Discrete unit (A3)
   │      │
   │      └── Unit boundary as the only direction-reversal locus (A9)
   │
   ├── Scalar / quantity / magnitude (104)
   │      │
   │      └── Scalar ratio (A:B with three orientations: A<B, A=B, A>B)
   │             │
   │             └── Cross-ratio = ratio of ratios = projective invariant (104)
   │                    │
   │                    └── Scalar motion (A4) — cross-ratio with space/time aspects
   │                           │
   │                           ├── Two forms: speed s/t (yang/material) ↔ energy t/s (yin/cosmic)
   │                           │
   │                           └── Scalar dimension (cross-ratio with unity datum)
   │                                  │
   │                                  └── Three-dimensional motion (A7, Nehru)
   │
   ├── Natural reference system (unit speed = c) — progression (A5)
   │      │
   │      ├── Material sector (s/t, clock time, observable)
   │      └── Cosmic sector (t/s, clock space, inverse)
   │
   ├── Tristate dimension of motion (A6: unity / speed / energy)
   │      │
   │      ├── Six units of motion across 3 dimensions
   │      ├── Three speed ranges 1-x, 2-x, 3-x
   │      └── Displacement = (count − 1) from unit-speed datum
   │             │
   │             └── A-B-C displacement notation
   │                    │
   │                    ├── Magnetic rotation (A, B): t/s × t/s = t²/s²  — 2D
   │                    └── Electric rotation (C):    s/t              — 1D
   │
   ├── Quantum π = 4 (105) — discrete-unit perimeter / diameter
   │      │
   │      ├── Circumference 8r ↔ electric rotation (1D, perimeter)
   │      └── Area 4r² ↔ magnetic rotation (2D, area = Periodic Table 4n²)
   │
   ├── Mass and gravity (107)
   │      │
   │      ├── Mass = magnetic × electric rotation = (t/s)² · (s/t) = t³/s³ · t
   │      ├── Mass / gravity full reciprocity (4-fold inversion: t↔s, out↔in, angular↔linear, volume↔point)
   │      ├── Step measure ↔ growth measure conversion: Δs = ln(Δt)
   │      └── Gravity threshold: gravity = floor( ln(Δt) ); massless when Δt ∈ {0, 1, 2}
   │
   ├── Lorentz unit-circle (108)
   │      │
   │      ├── 1 = γ² + v²  (Pythagorean / unit-circle equation, c = 1)
   │      ├── Three special points on unit circle: (+1,0) progression / (0,±γ) photon-birotation / (−1,0) gravity
   │      ├── ± double-root involution (negative branch normally suppressed)
   │      └── c-progression as fulcrum between motion-in-space and motion-in-time
   │
   └── Quaternion / division-algebra extension (108, 109)
          │
          ├── Hamilton: i² = j² = k² = i·j·k = -1; i·j = k
          ├── Four units of motion in one scalar dimension: {+1, i, i·j, -1}
          ├── Birotation: (k)(-k) = +1 ⟺ (e^(ix) + e^(-ix))/2 = cos(x)
          └── Variable units per dimension hierarchy:
                 1D real      (progression)
                 2D complex   (single rotating system: e⁻, e⁺)
                 4D quaternion (double rotating system: photon, neutrinos, atoms)
                 8D octonion  (life unit / α-helix)
```

The **photon** is the central derived object in RS2: a quaternion birotation (i·j paired with -k) that traces a cosine waveform at unit speed. The **proton/atom** is a dual quaternion (two double rotating systems). Charges are electric (1D) or magnetic (2D) rotational vibrations.

---

## 3. Glossary of named RS2 objects

| Term | Definition | Paper |
|---|---|---|
| **Aspect** | A reciprocal facet of motion: space (s) and time (t) are the two aspects. | 102 |
| **Affine stratum** | Geometric stratum of the biologic level; scalar recursion; Euclid's axioms partially fail. | 103 |
| **Atom** | A "double rotating system" (two quaternion-like structures), notated A-B-C. A unit-sized chunk of 3D time arranged on a 3D spatial grid. | 102, 106, 107 |
| **Atomic Z (atomic number)** | Encoded in the A-B-C displacement triple; not a 1D integer in RS2. | 106, 109 |
| **Birotation** | Two oppositely directed rotations whose sum is a cosine (Euler). The photon. In RS2 quaternion form: i·j paired with -k. | 107, 108, 109 |
| **Charge (electric)** | 1D rotational vibration. Occupies a scalar dimension. | 107 |
| **Charge (magnetic)** | 2D rotational vibration. | 107 |
| **Clock** | The unit-speed datum; the "tick" by which all motion is measured. Two forms: clock time (material), clock space (cosmic). | 106 |
| **Cosmic sector** | Reciprocal sector: 3D time + clock space + inverse (energy t/s) frame. Contains "cosmic matter" (= antimatter by inversion of scalar direction). | 102, 106 |
| **Cross-ratio** | Ratio of ratios; the only projective invariant across all geometric strata. Underlies scalar motion. | 104 |
| **Direction reversal** | Larson's 1D construction: a "diameter" on which to spin a rotation. Retired in RS2 (rotation is primary). | 107, 108 |
| **Dimension (scalar)** | A cross-ratio where one scalar orientation is fixed at unity and the other varies. | 104 |
| **Dimension of motion** | Tristate object (unity / speed / energy). RS2 has exactly three. | 106 |
| **Direction reversal (Lorentz)** | In the unit-circle picture: motion from (+1, 0) to (0, ±γ) — onto the imaginary axis. | 108 |
| **Discrete unit** | The minimum quantity of anything (= 1). No zero, no fractional, no negative. | 101, 104 |
| **Displacement** | Difference between count and unity, measured from the unit-speed datum. Temporal "n" or spatial "(n)". | 106 |
| **Division algebra** | The 1, 2, 4, 8 algebras (real, complex, quaternion, octonion). Each maps to an RS2 unit-of-motion structure. | 109 |
| **Equivalent space / equivalent time** | Spatial expression of temporal motion (and vice versa). Effectively the "imaginary axis" of RS2. | 106, 107, 109 |
| **Force field** | Apparent interaction caused by scalar-induced location changes. Lines of force = t/s². No real interaction; it is independent motions appearing to interact. | 104 |
| **Gravity** | Inward, linear velocity in space, expressed at a point — the full reciprocal of mass. Equals "the speed of light running backwards" at the (-1, 0) coordinate of the Lorentz unit circle. | 107, 108 |
| **Growth measure** | Logarithmic / integral counting from unity to infinity (1 → ∞). Paired with energy unit. | 107 |
| **Imaginary aspect** | Equivalent space (or time) — the v² / orbital component. | 109 |
| **Levels of existence** | Inanimate (Euclidean stratum), Biologic (affine stratum), Ethical (further projective stratification). RS2 derives these from geometric stratification (replacing Larson's 1995 metaphysical postulates). | 102, 103 |
| **Life unit** | Motion *between* space and time (rather than in either); octonion structure. | 103, 109 |
| **Lorentz factor** | Reinterpreted as unit-circle equation 1 = γ² + v²; both ± roots kept. | 108 |
| **Magnetic rotation** | 2D rotation, t²/s² dimensions, area-like (4r² = Periodic Table 4n²). Spin-½ in conventional physics (period 4π). | 105, 106, 107 |
| **Magnitude** | "Greatness of size or extent." Counting-number value, ≥ 1. | 104 |
| **Manifest realm** | The dual material+cosmic sector simultaneity. | 102 |
| **Mass** | Outward angular velocity in time over a volume. RS-formula: magnetic × electric rotation = t³/s³ · t. | 107 |
| **Material sector** | Our observable sector: 3D space + clock time + speed (s/t). | 102, 106 |
| **Motion** | The single component of the universe; ratio of space to time. | 101, 102, 103 |
| **Natural datum** | Unit speed (c). The reference point from which displacements are measured. | 104, 106 |
| **Natural reference system / Progression** | Default outward expansion of the universe at unit speed. The "tick of the clock." Hubble expansion. | 106 |
| **Octonion** | 8D non-associative division algebra. Maps to "life unit" / α-helix structure. | 109 |
| **Orientation (scalar)** | The three categorical relations of two magnitudes: A<B (low), A=B (unit), A>B (high). | 104 |
| **Pixelation / quantization** | The Universe as a 3D grid of unit cubes, the geometric consequence of A3 in Euclidean projection. | 105 |
| **Photon** | Birotation; Euler cosine in 1D Larson, quaternion (i·j coupled to -k) in RS2. Carried at unit speed by the progression. | 105, 107, 108, 109 |
| **Projective geometry** | The base stratum of RS2 geometry; affine and Euclidean are sub-cases. | 103, 104 |
| **Quantum π** | Discrete-unit perimeter/diameter ratio = 4.0 (vs. analog 3.14159 and transitional 4 → 3.14159). | 105 |
| **Quaternion** | 4D division algebra {1, i, j, k}. Models photon, particles, atoms in RS2. | 108, 109 |
| **Rotational base** | Larson's construction: gravitational opposition to the progression at (-1, 0) on the Lorentz unit circle. Retired/replaced in RS2 by every-location-rotational quaternion structure. | 107, 108 |
| **Scalar** | Quantity possessing only magnitude. No "of what." | 104 |
| **Scalar boundary** | Unity (A=B); the natural datum separating high and low orientations. | 104 |
| **Scalar motion** | The projectively invariant cross-ratio with space/time aspects. | 104 |
| **Sector** | One of two simultaneous frames (material / cosmic). | 102, 106 |
| **Speed range (1-x, 2-x, 3-x)** | Three regions of net motion: low / intermediate / ultra-high. | 106 |
| **Speed / energy** | Two reciprocal forms of motion: s/t and t/s. | 102, 104, 106 |
| **Step measure** | Linear / counting measurement, range 0 → 1. Paired with speed unit. | 107 |
| **Tao of motion** | Yang (linear) + yin (angular) as joint primaries. | 107 |
| **Tristate dimension** | Unity / speed / energy. | 106 |
| **Unit boundary** | Inter-unit juncture where direction reversal can occur. | 105 |
| **Unit speed** | The natural datum; speed of light. | 104, 106 |

---

## 4. Honest assessment — non-rigorous, hand-wavy, or contradictory items

These flags are for Joe's transparency. The papers contain real mathematical structure but also substantial gaps where claims are asserted, not derived.

1. **Dropping of the third postulate is asserted, not proved (102).** Larson states the scalar nature of space-time is "a necessary consequence" of the previous postulates but the proof is not given. Peret carries this forward without re-examining.

2. **Nehru's three-dimension argument (106) is heuristic.** The equation `n(n-1)/2 = n` yielding {0, 3, ∞} is an algebraic stability count for "binary combinations within a closed group," but the framing assumes the group structure {1, i, j, k} *as the model* — it is not a theorem-grade derivation that physical space *must* be 3D. The 0 and ∞ solutions are excluded by hand.

3. **Quantum π = 4 (105) is reference-system relative, not a universal constant.** The paper acknowledges three values (analog 3.14159, quantized 4.0, transitional). The claim that "Larson's 4n² Periodic Table formula = πr² with π = 4" is a numerical coincidence reframed as a structural identity. The paper does not justify why the discrete grid is square (rather than hexagonal, e.g., which would give a different "π").

4. **Force fields as "no actual interaction" (104) is metaphysically strong.** "There is actually no interaction at all — they are independent motions that, when observed, appear to be interacting." This is asserted but not formally connected to the equations of electromagnetism or gravitation.

5. **Mass/gravity reciprocity (107) is constructed by selecting four inversions (t↔s, out↔in, angular↔linear, volume↔point).** Why these four (and not others) constitute "the" reciprocal is not derived. The fourth (volume ↔ point) is geometric and only justified post-hoc.

6. **Massless threshold ln(Δt) (107) uses a logarithm of a discrete count.** The growth-measure / step-measure conversion `Δs = ln(Δt)` is given by citation (Larson, *Basic Properties of Matter*, Eq. 1-1) but not derived in the RS2 papers. The threshold ln(3) ≈ 1.1 is offered as the explanation for proton being the first massive particle — but the rounding `gravity = floor(ln(Δt))` is asserted, and the small overshoot (1.1, not exactly 1) is not addressed.

7. **Lorentz factor as unit circle (108) drops β, conflates v and v/c, and treats c = 1 as a unit conversion.** This is consistent within the system but the rewriting `1/γ = sqrt(c² - v²)` instead of conventional `γ = 1/sqrt(1 - β²)` requires the reader to track that "γ" has been redefined on the y-axis as `sqrt(1 - v²)` (the *reciprocal* of the conventional γ). This is a notational hazard.

8. **Negative root of Lorentz factor (108).** The argument that suppressing the negative root creates a 0/0 problem is illustrated by the joke 2 = 1 paradox — but the paradox uses cancellation of (a² - ab) = 0, which is a different issue from `± sqrt(1 - v²)`. The two are conflated rhetorically.

9. **Quaternion birotation (108, 109).** The mapping of the four units of motion {+1, i, i·j, -1} to {progression, electric, magnetic, gravity} is offered without showing that this *uniquely* identifies them — many other quaternion structures could be assigned. The argument is "this works and is consistent."

10. **Octonion = life unit (109)** is a mapping by analogy. The paper does not derive that the 8D octonion produces "α-helix" structures; it asserts the identification.

11. **Two postulates are silent on what selects unit speed.** A1 introduces "discrete units" but does not specify *which* magnitude gets the label "unit." A5 declares unit speed the natural datum. These are not the same statement, and the second is not derived from the first.

12. **Projective stratification produces "seven possible levels" (103).** No explanation of why seven; the three documented levels (inanimate, biologic, ethical) are described, the remaining four are not.

13. **Larson's "ultra-high speed" 3-x range (106) is replaced in RS2 by the quaternion -1 unit (109).** The two pictures are not formally reconciled; the speed-range table in 106 and the quaternion table in 109 carve up the same domain differently.

14. **Charge "vibration created by a photon captured in a rotation" (107).** This is a verbal model of how charge arises but the precise capture mechanism is not specified.

15. **Spin-½ as "cone sweeping a sphere" (107).** This is a visualization, not a derivation. The 720° period is asserted to match conventional spin-½, but the connection to fermion statistics, exchange symmetry, etc., is not made.

16. **No explicit derivation of the Riemann critical line σ = ½.** None of the nine papers contains the constant 1/2 used in the critical-line sense, the Riemann zeta function, or any direct connection to prime distribution / non-trivial zero structure. The closest structural analogues are:
    - Spin-½ from the magnetic 4π rotation (107).
    - The ± involution of the Lorentz factor's two roots (108).
    - The reciprocal involution s/t ↔ t/s with c as fulcrum, fixed only at |s/t| = 1 (108, 109).
    - The middle quaternion unit i·j as the "halfway" point in {+1, i, i·j, -1} (109).
    - Unit-circle structure with three distinguished points (108).
    None of these is presented as a Riemann-style fixed-point or critical-line argument; the RS2 corpus does not connect to the zeta function in 101–109.

---

## 5. Riemann-relevance summary

For Joe's cold derivation, the structural seeds available *in the foundation papers* (excluding the prior-art directory):

- **Reciprocal involution s ↔ 1/s** with unique fixed point at s = 1 (i.e., unit speed). This is the *only* fixed point of the involution in RS2's natural-framework numbers (≥ 1, finite). Note: the "boundary at unity" is at *one*, not at *one-half* — RS2 does not place its critical point at ½. Any σ = ½ derivation needs a *second* construction layer.
- **Lorentz unit circle** with three special points (+1, 0), (0, ±γ), (-1, 0). The middle point — the photon birotation — is the only point at v = 0 (i.e., the imaginary axis only). This is the structural analogue of "real part = 0" — but Riemann's critical line is real part = ½, not 0.
- **Spin-½ as 4π / 2 period** (107) — a 1/2 appears as the *result* of identifying a 4π-period rotation with what conventional physics measures as 2π-period spin. The factor of ½ is structural to the magnetic rotation. This is the most direct "½" appearance in RS2.
- **Three-points-of-unit-circle ↔ 3D rotational stability** (108, 109).
- **Cross-ratio as projective invariant** (104) — historically, cross-ratios are central to the projective-geometric formulation of the upper half-plane and modular forms (which is the natural home of zeta-like objects). RS2 does not invoke this connection but the cross-ratio primitive is present.
- **n(n-1)/2 = n stability count** (106) — this is structurally an *involution* (n ↔ n-1) plus a halving. Joe should consider whether this halving is the source of the ½.

The most promising bridge for the cold derivation: identify the magnetic-rotation 4π → 2π halving (RS2-107 spin-½) and/or the Nehru `n(n-1)/2 = n` halving as the source of σ = ½, with the projective cross-ratio (104) and Lorentz unit circle (108) providing the involution structure.
