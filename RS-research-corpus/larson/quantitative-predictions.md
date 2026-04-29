# Larson's Quantitative Predictions — A Reference Catalog

This file lists every numerical prediction Larson (or his immediate followers Peret/Satz/Hatch) makes, with comparison to observation and the source chapter. The purpose: give Joe a reference for what "RS quantitative result" looks like — the standard for Hard Problems Paper #2 to match.

## Natural unit values (no observational comparison — definitional)

These are RS-derived unit values, not predictions per se. Used to convert RS dimensional quantities to SI.

| Quantity | Larson value | Modern (CODATA) | Agreement |
|---|---|---|---|
| Speed of light | 2.997930 × 10¹⁰ cm/s | 2.99792458 × 10¹⁰ | <10 ppm |
| Atomic mass unit (g) | 1.65979 × 10⁻²⁴ | 1.66054 × 10⁻²⁴ | 0.05% |
| Gravitational G (cgs) | 6.67537 × 10⁻⁸ | 6.67430 × 10⁻⁸ | 0.02% |
| Boltzmann k (erg/K) | 1.38044 × 10⁻¹⁶ | 1.38065 × 10⁻¹⁶ | 0.02% |
| Electric potential nat. | 9.31146 × 10⁸ V | n/a (definitional) | matches 1 amu = 931.494 MeV/c² to 0.04% |

Source: NBM Ch 13, BPM Ch 5, BPM Ch 20.

## Atomic mass predictions (BPM Ch 24)

Larson's worked example for lead:

| Element | Z | G | Theoretical (m_v + 2Z) | Observed | Error |
|---|---|---|---|---|---|
| Lead (²⁰⁷Pb) | 82 | 43 | 207 | 207.2 | 0.1% |

Larson Eq 24-1: `m_v = I × m_r² / 156.444`. He claims agreement at the 0.1–0.2% level for stable isotopes throughout the periodic table; specific table CVII gives the full set but was not WebFetch-extractable in this pass.

## Hadron mass predictions

### NBM Ch 15 — Cosmic Ray Decay

| Particle | Larson calc | Observed (1979 era) | Agreement |
|---|---|---|---|
| Pion (π) | 137.95 MeV | 139.57 MeV | 1.2% |
| Muon (μ) | 106.42 MeV | 105.66 MeV | 0.7% |
| Cosmic H² (c-H²) base | 1848.61 MeV | n/a (theoretical) | — |
| J/ψ as c-H² + 2 charges | 3710.91 MeV | 3695 MeV | 0.43% |

### Satz, "Identification of Cosmic Particles"

| Particle | RS identification | RS calc | Observed | Agreement |
|---|---|---|---|---|
| ψ(3105) | c-He³ + 2 grav charges | 3104 MeV | 3105 MeV | 0.03% |
| J/ψ(3695) | c-H² + 2 isotope charges | 3710 MeV | 3695 MeV | 0.43% |

Plus older identifications mentioned (less detail in WebFetch pass): muon, pion, lambda, sigma, xi, omega — masses 106 to 1673 MeV/c², all identified as cosmic chemical-element isotopes.

### NBM Ch 16 — Cosmic Atom Building (kaon)

> "half of the c-Kr mass (52 MeV), and half of the 931 MeV mass" totaling 494 MeV

Kaon: predicted 494 MeV; observed K-meson mass is 493.7 MeV. Agreement: <0.1%.

> "successive members of which differ by 52 MeV"

A ladder of cosmic-element rest masses with 52 MeV spacing.

> "advancement from singly-charged to doubly-charged states requires **238 MeV** jumps"
> "motion into second dimensions demands **~1500 MeV** increments"

These are RS's predictions for resonance ladders.

## Atomic structure predictions

### Inter-atomic distances (BPM Ch 1, 2)

Eq 1-7: `s₀ = 2.914 × ln(t) Å`. Tables 2-6 in BPM Ch 2 give specific computed vs observed inter-atomic distances. Larson admits limits:

> "we are not yet in a position where we can determine specifically just what the inter-atomic distance will be for any given element under a given set of conditions"

So these are not strict predictions but post-hoc fits. Quality varies element to element.

### Periodic-table block sizes

NBM Ch 10's 2n² formula yields: 2, 8, 18, 32 magnetic + 4n² coupling: **inert gases at 2, 10, 18, 36, 54, 86, 118**. This matches the empirical periodic table exactly (with Z=118 being oganesson, only synthesized 2002, well after 1979 NBM publication). Strong prior-prediction.

### Stability limit

NBM Ch 11 + BPM Ch 24: "Maximum magnetic rotational displacement is four units. Elements up to atomic number 117 remain stable at zero ionization; element 118 and above are inherently unstable." Modern data: oganesson (Z=118) has half-life 0.7 ms. Tennessine (Z=117) is also short-lived. The boundary is approximately right.

## Cosmological predictions (UoM)

UoM is largely qualitative. Specific quantitative predictions are limited:

- **Background radiation isotropy** — predicted to come from cosmic-sector matter, not from Big Bang. Observed: yes, isotropic. Source mechanism: contested.
- **No primordial element abundance prediction** — RS rejects nucleosynthesis cosmology, replaces it with continuous cosmic-atom-building. No specific H/He abundance prediction.
- **Hubble redshift = recession beyond unit speed** — RS predicts redshift increases unbounded with distance, and at z > 1 corresponds to speeds > c (unobservable in space, observable as redshift only). Modern observations: high-z quasars exist with z > 7, consistent.
- **Galaxy speeds 8-10× c** — RS prediction for the most distant quasars. Observed: large redshifts consistent with this if interpretation is "scalar speed" not "Doppler velocity".

## Lifetime predictions

NBM Ch 15:

> "Each dimension of motion modifies the unit of time applicable to the particle life by approximately 10⁻⁸, while each gravitational charge modifies the unit by about 10⁻²."

So a particle with d dimensions of motion and g gravitational charges has lifetime ~ `t_unit × 10⁻⁸ᵈ × 10⁻²ᵍ`. With t_unit = 1.52 × 10⁻¹⁶ s:
- d=1, g=0: ~10⁻²⁴ s (typical hadron resonance)
- d=2, g=0: ~10⁻³² s (very short)
- d=1, g=1: ~10⁻²⁶ s
- d=2, g=2: ~10⁻³⁴ s

These align in order-of-magnitude with observed particle lifetimes across the resonance / strange / charm hierarchy. Not a precise prediction, but a structural one.

## Hydrogen atom

NFoS Ch 6:

> Net speed `1 - 1/n²` "explicitly recognized as matching hydrogen spectral frequencies"

This is the Balmer/Lyman/Rydberg formula. Larson's account: the n in `1/n²` is the rotational displacement quantum number, the `1 -` is the gravitational-progression baseline. Reproduces the entire hydrogen spectrum exactly (Rydberg formula).

## Compressibility and bulk properties (BPM Ch 4)

`P₀ = az / 938.67` and downstream equations give bulk modulus, compressibility, etc., in terms of atomic-number z and a tabulated structural parameter `a`. Tables in BPM Ch 4 show fits to observed data; quality ~1–10% typical. (Tables not fully extracted.)

## Gas constant

BPM Ch 5: `R = (2/3)` natural unit of specific heat. Modern: R = 8.314 J/mol/K. The (2/3) factor traces to dimensional reduction of three-dim gas thermal energy to one-dim measurement.

## What's MISSING — predictions Larson doesn't give

Critical for understanding what RS doesn't (yet) compute, and which Joe's cold rederivations would need to provide for Hard Problems #2–#6:

1. **Yang-Mills mass gap** — not in canonical Larson. The closest analog is "mass exists iff three-dim rotational displacement exists" (NBM Ch 11), but no quantitative bound on the gap.
2. **Hierarchy problem** — Larson does not address why electroweak/Planck scale ratio is 10⁻¹⁷. RS uses unit speed as the only fundamental scale; hierarchy emerges from displacement quantum numbers.
3. **Fine-structure constant** — Larson does not derive α = 1/137. RS has no equivalent of the EM coupling constant as a calculated number; α is implicit in the magnetic mass coefficient 0.006392 (which is approximately `α/(2π × something)` numerically, but not in canonical RS).
4. **Gauge couplings** — RS has no gauge group structure. The three forces (gravity, electric, magnetic) are not unified by a gauge group; they emerge from dimensional reduction of scalar motion.
5. **Master formula for particle physics** — Peret's RS2 attempts this; canonical Larson does not have a single master equation.

## Standard for Hard Problems Paper #2

Based on canonical Larson's track record, an RS-style YM mass-gap prediction should:

- **Predict the gap as a multiple of natural units**, ideally with the natural unit of mass (931 MeV) appearing as the base scale.
- **Achieve <1% agreement** with the lattice-QCD estimate (~1.5 GeV for the lightest glueball; experimental glueball mass not yet established).
- **Use only the canonical constants**: 156.444, 0.006392, 938.67, 931.146 MeV, integer rotational displacement quantum numbers.
- **Not introduce new free parameters**.
- **Optionally**: introduce an `ln` factor via Eq 1-1 if a confined-structure integration is invoked.
- **Optionally**: introduce a `2π` factor only if RS2/Peret birotation period is invoked — and flag that as an extension beyond canonical Larson.

Joe's `Δ = ln(2π) × 931.2 MeV ≈ 1.711 GeV` headline:
- The 931.2 is canonical (BPM Ch 20 natural unit).
- The ln(2π) requires either:
  - A Peret-RS2 birotation period 2π, or
  - A novel application of BPM Eq 1-1's ∫dt/t = ln(t) with some rotational period 2π set as the integration limit.

If the cold derivation can produce 2π without invoking RS2 — i.e., from canonical Larson alone — that is a stronger result than the prior Opus 4.5 derivation. If not, the result is RS2-dependent and should be cited as a Peret-Nehru extension (which is honest given Joe's RS2 framing).
