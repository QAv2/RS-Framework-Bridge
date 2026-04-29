# Nuclear & Mass — Focused Extract for Yang-Mills Cold Rederivation

This is the **single most important file** in this corpus for the Yang-Mills mass-gap derivation. It pulls together every quantitative claim in the Nehru / Peret / Satz / Vijaya canonical literature about: strong nuclear force, hadron mass spectrum, atomic-zone partition, and unit conversions.

If the headline `Δ_YM = ln(2π) × 931.2 MeV ≈ 1.711 GeV` is going to survive cold rederivation, the constants and counting arguments below are where it must come from.

---

## §1. The 931.15 MeV constant — what it is in RS

This is **the AMU-to-MeV conversion**, and in RS it is identified with the **unit-mass increment per gravitational charge** in the time region. From Nehru's *Internal Ionization and Secondary Mass*:

> "The space-time direction of the gravitational charge is altogether different. The third region, in which the motion of this charge takes place, turns out to be toward 'our side' of the time region, rather than the far-side, and therefore coincides with the region outside unit space (represented by line 'O' in Figure 1(c)). Thus the net interregional ratio applicable to the gravitational charge is 1. Consequently the secondary mass contribution of the gravitational charge is one full unit: **931.152 MeV**."

This is repeated in *High Energy Physics and the Reciprocal System* as the prefactor of the cosmic-atom mass formula:

> "the mass equivalent of the corresponding c-atom would be a function of 1/A. In fact it is given by `(G + 4/Aᶜ) × 931.15 MeV`"

**For the Yang-Mills problem**: any RS-derived hadron mass scale will have 931.15 MeV (≈ AMU) as a building block. This is the Compton wavelength reciprocal of the proton — i.e. it's the hadron mass scale itself. Joe's prior derivation `Δ ≈ ln(2π) × 931.2 MeV` uses this exact constant. The factor `ln(2π)` is a **dimensionless multiplier** — it must come from a counting argument internal to RS.

### Where could `ln(2π) ≈ 1.83788` come from?

Candidate provenances visible in the canonical literature:

1. **Spin / phase argument**: Nehru's Spin paper notes that 1D rotation has period `2π`, 2D rotation has period `4π` (steradian). The `ln(2π)` factor is the natural log of the 1D period. **Plausible if the mass gap is a 1D-zone vs 2D-zone phase action.**

2. **Vibrational d.o.f.**: Inter-regional ratio's vibrational correction is `1 + 2/9` for atoms or `1 + 1/9` for sub-atoms. None of these straightforwardly give 1.838.

3. **Quaternion vs complex measure**: Time-region speed has dimension `1/t²` for 1D rotation but `1/t⁴` for 2D rotation. The dimensional ratio is `t²` and the natural log of this dimensional crossing for some unit specification could yield such a factor. **Worth investigating.**

4. **Berry-Keating / Riemann handover**: Joe's Riemann derivation already invoked `2π` heavily (the von Mangoldt counting `N(T) ∼ (T/2π) ln(T/2π)` for instance). If the same Hilbert-Pólya operator gates the Yang-Mills mass gap, `ln(2π)` is a natural carry-over.

5. **Gauge-invariant counting in finite Heisenberg group**: The metaplectic / Mp(2) double cover involves `2π` rotation cycles; one circuit of an Mp(2) projection is `4π` not `2π`, and `ln(2π) = ln(4π/2)` could emerge from a "half-orbit" counting.

The literature pulled here does **not** contain an explicit derivation of `ln(2π)` as a mass multiplier. **This is a gap.** Joe's prior Opus 4.5 derivation must contain the bridge.

---

## §2. Mass formula for hadron resonances (Nehru, *High Energy Physics*)

The headline equation:

> **m = (G + 4/Aᶜ) × 931.15 MeV**

where:
- `Aᶜ` = atomic weight of the cosmic-isotope being identified with the hadron
- `G` = number of gravitational-charge units (small non-negative integer or half-integer)

Both Tables 1, 2, and 3 of the paper match observed hadron masses to 1–2% across the charmonium / Σ baryon / meson resonance spectrum. Sample fits:

- π meson: c-Si²⁷, G=0 → **137.95 MeV** (calc) vs 139.57 MeV (obs)
- μ "lepton" (Larson called it muon, but RS treats it as a c-Argon decay): c-Ar³⁵, G=0 → **106.42 MeV** vs 105.66 MeV
- Λ baryon: c-Ne²⁰, G=1 → **1117 MeV** vs 1116 MeV
- Σ baryon: c-N¹⁴, G=1 → **1197 MeV** vs 1190 MeV (header), 1197 MeV (refined)
- Ω⁻ baryon: c-Li⁵, G=1 → **1676 MeV** vs 1675 MeV
- ψ(3100) charmonium: c-He³, G=2 → **3104 MeV** vs 3105 MeV
- ψ(3700) charmonium: c-H², G=2 → **3710 MeV** vs 3695 MeV

This is RS's substitute for QCD: hadrons are **inverted (cosmic) atoms** with a small integer count of gravitational-charge increments. The strong nuclear force as a fundamental interaction does not exist; what experimentalists measure as nuclear binding energy is the gravitational-charge contribution to the c-atom rotational mass.

### Half-integer G modes

Table 3 includes `½ c-Kr` at `1½` G → 1423 MeV (calc) vs 1427 MeV (obs). Half-charge modes exist when one of the two rotating systems of the c-atom carries the gravitational charge while the other doesn't.

---

## §3. Inter-regional ratios — the key dimensionless numbers

From Nehru's *The Inter-Regional Ratio* paper:

| Ratio | Value | Where it applies |
|---|---|---|
| `f₁` = degrees of freedom for 1D rotation in 3D space | 8 = 2³ | base count |
| `f_atom` = full atomic d.o.f. count | 4 × 4 × 8 = **128** | atom = (2D rot) × (2D rot) × (1D electric rot) |
| `R_atom` = atomic inter-regional ratio | 128 + 128(2/9) = **156.4444** | atomic-zone observations |
| `R_subatom` = sub-atomic inter-regional ratio | 128 + 128(1/9) = **142.2222** | sub-atomic-particle observations |

Variations that appear in derivations:
- `R_neutron-decay` = 128 (1 + 1/18) = **135.111** (used in *Lifetime of the Neutron*; halved vibrational contribution because only one of two opposite vibration directions yields antiparallel alignment)
- `8 × 156.44 / 7 = 178.79` (used in *Lifetime of c-Argon*; 8-unit displacement separation between positive and negative zero, 7-unit asymmetry between zero-speed-in-time vs zero-speed-in-space)

For the **Yang-Mills mass gap**, the relevant inter-regional ratio is most likely either 156.44 (atomic) or 128 (the bare d.o.f. count without vibrational correction). The value 156.44² = 24,473 appears in the magnetic-charge mass calculation in *Internal Ionization*. The combination `931.15 × (factor)` for various `factor` options:

- `931.15 / 156.44 = 5.95 MeV` — not obviously meaningful
- `931.15 / 128 = 7.27 MeV` — close to fission energies per nucleon
- `931.15 × ln(2π) = 1711.4 MeV` — Joe's headline gap value
- `931.15 / (156.44 × 2/9) = 26.79 MeV` — not obviously meaningful

The **`× ln(2π)`** path fits the headline; **the other arithmetic combinations don't**. This narrows the cold-rederivation hypothesis: the Yang-Mills mass gap in RS is `(natural unit mass) × (some dimensionless angular factor involving 2π)`. The factor is more likely a phase / counting log, not a pure ratio of inter-regional constants.

---

## §4. Mass components and unit conversion (Peret, *Subatomic Mass, Recalculated*)

Native RS mass-component table:

| Component | Symbol | Natural-unit value | Physical meaning |
|---|---|---|---|
| primary mass | `p` | 1.000000000000 | base mass of a single rotating system at rest |
| magnetic mass | `m` | 0.006392045455 | secondary mass from 2D rotational vibration (gravitational charge in materials) |
| gravitational mass | `p+m` | 1.006392045455 | combined p + m |
| electric mass (3-dim) | `E` | 0.000868055556 | secondary mass from 1D rotational vibration in 3-dim |
| electric mass (2-dim) | `e` | 0.000578703704 | (2/3) E |
| normal charge mass | `C` | 0.000044944070 | mass of normal electric charge |
| electron charge mass | `c` | -0.000029962713 | -(2/3) C |

Conversion factor (Peret's calibration via the charged electron):
> **0.99970644 u/n** (unified atomic mass units per natural-mass unit)

So 1 RS natural mass unit ≈ 0.99970644 u ≈ 0.99970644 × 931.49432 MeV/u = **931.215 MeV**.

This nearly-but-not-quite-1 conversion factor is a load-bearing piece of the Yang-Mills calibration. Joe's `931.2 MeV` headline rounds Peret's `931.215 MeV` (the natural-mass-unit equivalent in MeV after Peret's electron-mass calibration).

---

## §5. The atomic / nuclear zone partition

From Nehru's *'Quantum Mechanics' as the Mechanics of the Time Region* §5:

### Two zones inside the time region (i.e. inside one natural unit of space ≈ 1 fm)

**One-dimensional zone**:
- Where 1D rotations and basic photon vibrations live
- Time-region speed dimension is `1/t²` (second power of S-frame `1/t`)
- Speed-related quantities are SECOND-POWER expressions of S-frame quantities
- Atomic spectra (electronic energy levels) live here — they're the eigenstates of Schrödinger's equation in the `r²` potential well

**Three-dimensional zone**:
- Where compound atomic motion lives — *the atom IS this zone*
- Time-region speed dimension is `1/t⁴` (fourth power of S-frame)
- Speed-related quantities are FOURTH-POWER expressions
- "Nuclear" spectra (what conventional physics calls nuclear energy levels) live here
- Required mathematical apparatus: **quaternions, not complex numbers** (because the dimensionality is 4, not 2)

### The "nuclear potential" in RS — Eq.(31)

> **V_T3 = −K_P3 × (rₐₙ − r)⁴ + K_G3/r⁴ ± K_I3**

where `rₐₙ` is the atomic radius (= `1.2 × A^(1/3)` fm = exactly the conventional nuclear-radius formula). Components:
- **`-K_P3 × (rₐₙ - r)⁴`**: attractive well from space-time progression, peaks at the atomic boundary
- **`K_G3 / r⁴`**: REPULSIVE core from gravitation (RS-gravity is repulsive in time region)
- **`±K_I3`**: constant initial-level term

Nehru's claim:
> "An unexpected feature of the experimental data analysis was the occurrence of a repulsive core of small radius. The Reciprocal System, on the other hand, **actually predicts this repulsive core**, namely, V_G3."

This 4th-power potential well + 4th-power-inverse repulsive core IS the RS analog of the QCD confining potential. It is what holds the compound-motion atom together. Numerical fit-quality vs experimental nuclear potentials is described as "remarkably close qualitative resemblance" — not yet a quantitative claim.

---

## §6. The atom size and the absence of a nucleus

From Nehru's *Wave Mechanics in the Light of the Reciprocal System*:

> **rₐ = 1.2 × A^(1/3) femtometers**

This is *identically* the conventional nuclear radius formula. The RS reading: this is the size of **the atom itself**, not a nucleus inside the atom. There is no nucleus — Larson's *The Case Against the Nuclear Atom* (1963) is the foundation paper.

> "Calculations based on the inter-regional ratios applicable confirm Larson's assertion that the measured size of the atom is in the femtometer range and hence what is found from the scattering experiments is the size of the atom itself—not of a nucleus."

This is the strongest version of "RS rejects QCD entirely": the entire atomic scattering-cross-section is the cross-section of the compound-motion atom, which is a unit, not built of point-like quarks.

---

## §7. Gluon-analogs / SU(N) substitutes

There is **no SU(N) structure in RS**. The substitute is the **scalar-rotation displacement notation `a-b-c`**:
- `a, b` = displacements in two dimensions of basic 2D magnetic rotation (each contributes `2a²` electric units)
- `c` = displacement in 1D electric rotation (in parentheses if negative / time-displaced)

Electric equivalent of magnetic-displacement increments:
- 1st increment: 2 × 1² = 2 electric units → corresponds to 1st period of periodic table
- 2nd: 2 × 2² = 8 → 2nd period (lithium through neon)
- 3rd: 2 × 3² = 18
- 4th: 2 × 4² = 32

Total atomic number after n increments: `2 × (1 + 4 + 9 + ... + n²) = 2 × n(n+1)(2n+1)/6` — generates the periodic-table period lengths.

This `a-b-c` displacement notation **is what replaces SU(N) quantum numbers in RS**. In Nehru's hadron table, c-atoms are labeled by their rotational displacement triplet AND by the integer count G of gravitational charges they carry.

For a Yang-Mills derivation, the SU(N) coupling is replaced by:
- **Inter-regional transmission factor** `1/c` (1D), `1/c²` (2D), `1/c³` (3D) — this is the analog of the fine-structure constant
- **Inter-regional ratio** `R = 156.44` (atomic) or context-specific variants — this is the d.o.f. partition function

---

## §8. Connection to the Riemann derivation Joe just completed

The Riemann zeta-function derivation (Phase 5.1) is anchored on:
- The Berry-Keating Hamiltonian H = (xp + px)/2
- The metaplectic / Mp(2) double cover for time-reversal breaking → GUE not GOE
- The `2π` factor in the von Mangoldt counting `N(T) ∼ (T/2π) ln(T/2π)`
- The √x prime-wave amplitude (Dorsey convergence: [§5.4 reading h])

The Nehru/Peret canonical literature **doesn't connect explicitly to Riemann zeros**. But the framework is compatible:
- The 1D zone of the time region has speed dimension `1/t²` ↔ Mp(2) double cover
- The 3D zone uses quaternions ↔ SU(2) ↔ Berry-Keating's higher analogs
- The inter-regional ratio 128 = `2⁷` is suggestive of 7-bit Hilbert-space counting
- The natural unit of action `Eₙ × Tₙ ≈ 1.49 × 10⁻³ × 1.52 × 10⁻¹⁶ erg·s = 2.27 × 10⁻¹⁹ erg·s ≈ 1/2900 ħ` after correction by R=156.44 yields exactly Planck's constant

This last is the **Nehru-Planck check** in action — the Nehru *Theoretical Evaluation of Planck's Constant* paper derives `h = (Eₙ Tₙ Sₙ / R) × (1+s)^(1/3) = 6.6256 × 10⁻²⁷ erg·s` to 0.0002% accuracy. The factor `(1+s)^(1/3)` with `s = m + e_initial` accounts for primary-vs-inertial mass discrepancy.

For the Yang-Mills problem: **whatever counting gives the `ln(2π)` factor must be consistent with what gives `2π` in the Riemann derivation.** This is the cross-check Joe should run after cold-deriving Yang-Mills.

---

## §9. Open questions for the cold rederivation

Things this corpus does NOT contain that Joe's Yang-Mills derivation will need:

1. **Explicit derivation of `ln(2π)` as a mass multiplier.** Not in Nehru, not in Peret, not in Satz. The factor must come from a Berry-Keating / metaplectic argument unique to Joe's prior derivation.

2. **Glueball mass spectrum.** RS literature does not discuss glueballs (RS rejects gluons). Closest analog: the meson resonance table in *High Energy Physics*, but those are c-isotope identifications, not glueball candidates.

3. **Yang-Mills coupling running / asymptotic freedom.** Not addressed. RS uses inter-regional transmission factors instead of running couplings.

4. **Wilson loops / lattice gauge.** No analog. RS substitutes the time-region geometry.

5. **Color confinement.** Trivially satisfied in RS — there are no quarks to confine. The atom is the smallest unit.

6. **Mass gap proof (Yang-Mills mass gap is a Clay Millennium problem about pure SU(N) gauge theory in continuum limit).** RS does not have a Yang-Mills theory at all; it has the c-atom rotation framework. The "mass gap" in Joe's sense must mean something specific in the RS rederivation — perhaps the smallest non-trivial gravitational-charge increment, or the gap between G=0 and G=1 c-atoms?

If the gap is `(G=1 mass) − (G=0 mass)` for a fixed Aᶜ, that's exactly **931.15 MeV per gravitational-charge unit**, with no `ln(2π)` factor. So the headline `ln(2π) × 931.2` is NOT the simple G-spacing — there must be an additional log factor.

Hypothesis to test in the cold rederivation:
> "The Yang-Mills mass gap = (gravitational-charge unit-mass) × (action of a single rotation in the 1D zone of the time region in natural-units, expressed as a phase-action integral over `[0, 2π]`)"
> 
> = `931.15 MeV × ln(2π) ≈ 1711 MeV`

The action of a `2π` rotation in natural-log units is `ln(2π)` — this is a pure dimensionless log because phase is dimensionless. If the cold rederivation can independently establish that the mass-gap operator has eigenvalue `(unit mass) × (action of one phase cycle)`, the headline is reproduced.

---

## TL;DR for the cold rederivation

**Confirmed in canonical literature**:
- 931.15 MeV = mass per gravitational-charge unit ✓
- AMU calibration via Peret's electron-mass anchor → 931.215 MeV per natural-mass unit ✓
- Hadron spectrum reproduced as c-atoms × G with ~1–2% accuracy ✓
- 4th-power potential well + repulsive core matches phenomenological nuclear potentials ✓
- Atom size = 1.2 × A^(1/3) fm, identifying the atom with what conventional physics calls the nucleus ✓
- Quaternion math required for nuclear-energy-level treatment ✓ (Nehru 1995, Peret RS2-105+)

**NOT in canonical literature, must come from Joe's prior derivation chain**:
- The factor `ln(2π) ≈ 1.838` ✗
- The specific identification of "mass gap" with `(unit mass) × ln(2π)` ✗
- Quantitative glueball-equivalent prediction ✗
- Connection between Berry-Keating Hamiltonian (Riemann result) and the c-atom mass formula ✗

These are the **load-bearing claims** of the cold rederivation that aren't supported by the published Nehru/Peret corpus. They must be derived independently in the Yang-Mills paper draft.
