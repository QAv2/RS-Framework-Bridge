---
title: "Subatomic Mass, Recalculated — Peret 1995 (with 2018 Updated Units appendix)"
author_original: Bruce Peret
source_publication: Reciprocity Vol. 24 No. 2, 1995 (ISUS); update Reciprocity XXV/2 Autumn 1996; 2018 unit-calibration appendix
recovered_from: Google Drive (file IDs 1YbIVNsi… and 1FGlPvgD…)
distillation_by: Claude Opus 4.7 (cold read 2026-04-28 late evening)
purpose: Mass-derivation primary source for Hard Problems #2 Yang-Mills mass gap cold derivation
---

# Subatomic Mass, Recalculated — Peret 1995

Peret's quantitative validation of Larson's subatomic-mass framework. Calibrates RS natural mass units to atomic mass units (u) via the charged electron, computes proton, hydrogen, neutron, and neutrino masses to ≤ 0.011% error against 1994 *Physical Review D* values.

This is the **load-bearing quantitative result** for the RS framework: Larson's 1959 calculations for hadronic mass survived 35 years of experimental refinement at the 0.01%–0.04% level.

---

## §1. Mass component table (natural units)

Peret reproduces Larson's mass-component breakdown from *Nothing But Motion* p. 164:

| Symbol | Component | Natural units (n) | Atomic mass (u) |
|---|---|---|---|
| **p** | primary mass | 1.000 000 000 000 | 0.999 706 441 403 |
| **m** | magnetic mass | 0.006 392 045 455 | 0.006 390 169 015 |
| **p+m** | gravitational mass | 1.006 392 045 455 | 1.006 096 610 417 |
| **E** | electric mass (3D) | 0.000 868 055 556 | 0.000 867 800 730 |
| **e** | electric mass (2D) | 0.000 578 703 704 | 0.000 578 533 820 |
| **C** | mass of normal charge | 0.000 044 944 070 | 0.000 044 930 876 |
| **c** | mass of electron charge | −0.000 029 962 713 | −0.000 029 953 917 |

Conversion factor (u/n): **0.99970644** — derived from the ratio of the measured charged-electron mass to the calculated charged-electron mass:

$$\frac{0.000\,548\,579\,90\,u}{0.000\,548\,740\,99\,n} = 0.999\,706\,44\,u/n$$

`p` is the **t³/s³ inductance quantum** of one nucleon — RS calls this "primary mass." Numerically, p ≈ 1 amu by Larson's calibration choice. The conversion factor 1 u = **931.494 32 MeV/c²** is taken from Physical Review D 1994 footnote 4 (page 1396).

---

## §2. Particle composition formulas

Peret (after Larson) builds composite particles by addition:

| Particle | Composition | RS-calc (u) | Observed (u) | Error |
|---|---|---|---|---|
| Charged electron | e + c | 0.000 548 579 90 | 0.000 548 579 90 | 0.0000% |
| Charged positron | e + c | 0.000 548 579 90 | 0.000 548 579 90 | 0.0000% |
| (uncharged) Proton | p + m + 2e | 1.007 253 678 06 | 1.007 276 47 | −0.00228% |
| Charged proton | p + m + 2e + C | 1.007 298 608 93 | 1.007 276 47 | +0.00221% |
| **50/50 mix proton** | (avg of charged + uncharged) | 1.007 276 14 | 1.007 276 47 | **−0.000033%** |
| Hydrogen 1H | p + m + 3e | 1.007 832 211 88 | 1.007 940 00 | −0.0107% |
| Compound neutron | p + m + 3e + E | 1.008 700 012 61 | 1.008 664 904 | +0.0035% |

The 50/50 mixed-proton interpretation says the laboratory proton sample is 50.7% charged + 49.3% uncharged (mass-weighted), and reproduces the observed value to **0.000033%** — *below* the stated experimental error.

Hydrogen and compound-neutron predictions are within 0.011% of experiment.

**This level of agreement is the strongest quantitative case for the RS framework as physics**: an axiomatic system from 1959, with one calibration number (electron mass), reproducing nucleon masses to four decimal places.

---

## §3. Updated unit space and unit time (Peret 2018)

From the 2018 appendix at `peret_units_mass_integration_report.txt`. Calibration via the Rydberg constant rather than the Rydberg frequency of hydrogen:

| Quantity | Value (cgs) | Source |
|---|---|---|
| Speed of light c | 2.997 924 58 × 10¹⁰ cm/s (exact) | SI definition |
| Rydberg constant R∞ | 1.097 373 153 4 × 10⁵ cm⁻¹ | 1986 CODATA |
| **Unit space** | **4.556 335 267 × 10⁻⁶ cm** | (1/R∞) / 2 (half-cycle of Rydberg wavelength) |
| **Unit time** | **1.519 829 851 × 10⁻¹⁶ s** | Unit space / c |

Unit space corresponds to **half-cycle of the Rydberg wavelength** — the fundamental atomic spectroscopic length scale.

(Larson's 1959 values, derived from the Rydberg *frequency* of hydrogen, were 4.558 816 × 10⁻⁶ cm and 1.520 655 × 10⁻¹⁶ s — within 0.05% of the 2018 update.)

**Honest gap (from Joe's `peret_units_mass_integration_report.txt`)**: direct dimensional reduction from unit space + unit time to mass units does NOT reproduce the 1 amu = 1.66 × 10⁻²⁷ kg observed value. The naive `mass_unit = (unit space)³/(unit time)³ × (something)` calculation gives 1.69 × 10⁻³³ kg, off by ~10⁶. The mass-unit calibration is **empirical via the charged electron**, not via dimensional reduction. This is a load-bearing concession for any cold derivation that uses Larson's mass scale: the dimensionless mass *ratios* are first-principles, but the **absolute mass scale** requires one calibration input.

---

## §4. Implications for the Yang-Mills mass gap (Hard Problems #2)

The headline target for the RS2 cold derivation is `Δ_YM = ln(2π) × 931.2 MeV ≈ 1.711 GeV`. From this paper:

1. **Primary mass p = 1 amu (calibrated)** is the natural mass scale of the t³/s³ inductance quantum. This is the scale that goes into the YM formula, not the magnetic mass m (= 5.95 MeV per atomic constituent) and not the electric masses e, E, C, c (sub-MeV).

2. **Calibration is empirical (charged-electron-anchored), not dimensional.** The "1 amu" in `Δ_YM = ln(2π) × 1 amu` carries the empirical input. The dimensionless prediction `Δ_YM / m_amu = ln(2π)` is what's actually first-principles.

3. **Pure YM = no quark content.** In Peret's table, the proton has composition p + m + 2e + C. A pure-magnetic-only excitation (no electric, no charge) has composition p + m only — gravitational mass = 1.006 amu = 937.4 MeV. But the Yang-Mills mass gap is *additional excitation above vacuum*, not the rest mass of a constituent — it's the lowest stable energy quantum of pure magnetic-2D rotation.

4. **The lowest YM excitation is one full birotation of the primary mass quantum.** One full birotation = 2π in step measure (gluons are spin-1 = 1D rotation per Nehru's "Some Thoughts on Spin" §1, even when bound into 2D-rotation glueballs). Convert to growth measure (energy) via Δs = ln(Δt): the lowest growth-measure mass excitation is `ln(2π) × p`.

5. **Numerical prediction**: ln(2π) × 1 amu = 1.838 × 931.494 MeV = **1.711 GeV**. Compare lightest scalar glueball candidates: f₀(1500), f₀(1710), f₀(2020). The ln(2π) × m_amu value matches f₀(1710) most closely. Lattice-QCD 0++ glueball predictions cluster at 1.65–1.75 GeV (Morningstar–Peardon 1999, Chen et al 2006) — consistent.

---

## §5. Open questions raised by this paper

1. The neutrino mass corrections (e/2 muon-neutrino, (C+c)² electron-neutrino) are post-hoc fits to 1993 data. They flag conceptual issues in Larson's original neutrino interpretation. Are they re-derivable from first principles, or are they curve-fits?
2. The "50/50 charged/uncharged proton mixture" claim is a specific testable prediction (eliminating one charge-state should give exactly one of the two RS-calc values). Has anyone tested it experimentally?
3. The hydrogen mass error of 0.011% is larger than the proton error of 0.000033%. Why? Possible electron-binding-energy correction not accounted for in p+m+3e?
4. Direct dimensional reduction from unit space/time to mass unit fails by ~10⁶ — strong indicator that RS's unit-conversion chain is incomplete, and any physical prediction (including ln(2π) × 931 MeV) carries an empirical calibration element.
