# Larson Unit Conversion — RS Natural Units to SI

**Load-bearing for Yang-Mills cold rederivation.** This file pulls every place in canonical Larson where the s↔meter, t↔second, mass↔gram, and mass↔MeV calibrations are derived or asserted. The 931 MeV number Joe targets is a *natural unit* in RS, not an empirical fit.

## Source: NBM Ch 13 "Physical Constants"

URL: `https://reciprocalsystem.org/books/nbm/13-physical-constants`

**Table of natural units** (verbatim from chapter, abstracted by WebFetch but numeric values are Larson's):

| Quantity | Natural Unit |
|---|---|
| Space | 4.558816 × 10⁻⁶ cm |
| Time | 1.520655 × 10⁻¹⁶ sec |
| Speed | 2.997930 × 10¹⁰ cm/sec |
| Force | 7.316889 × 10⁻⁶ sec/cm² |
| Inertial mass | 3.711381 × 10⁻³² sec³/cm³ |

Note: Larson's 1979 unit of speed `2.997930 × 10¹⁰ cm/sec` agrees with modern c = 299792458 m/s = 2.99792458 × 10¹⁰ cm/sec to ~6 ppm. The natural unit of time is the reciprocal of this divided by the natural unit of space; equivalently `t_nat = s_nat / c`.

**Key derivation chain** (Larson):
1. Speed of light c = 1 (natural) by definition (the s↔t reciprocal at unit ratio).
2. Time unit derived from Rydberg fundamental frequency 6.576 × 10¹⁵ half-cycles/sec → `t_nat = 1/2 × (Rydberg)⁻¹ = 0.1521 × 10⁻¹⁵ sec` (SPU 1959 Ch 4 derivation; revised in NBM to 1.520655 × 10⁻¹⁶ sec, two-decimal-place agreement).
3. Space unit = `c × t_nat = 4.558816 × 10⁻⁶ cm`.

**Mass equivalence** (NBM Ch 13):

> "The reciprocal of this number, 1.65979×10⁻²⁴, in grams, is therefore the mass equivalent of unit atomic weight"

i.e., one atomic mass unit (1 u) ≡ 1.65979 × 10⁻²⁴ g, agreeing with modern 1 u = 1.66054 × 10⁻²⁴ g to ~0.05%.

**Mass component decomposition** (NBM Ch 13):

| Component | Value (atomic mass units) |
|---|---|
| Primary mass | 1.000000 |
| Magnetic mass | 0.006392 |
| Electric mass (3-dim) | 0.000868 |
| Electric mass (2-dim) | 0.000579 |
| Normal charge mass | 0.000045 |

**Gravitational constant**: `G = 6.67537 × 10⁻⁸` (cgs) — Larson, NBM 1979, vs modern CODATA 6.67430 × 10⁻⁸ (~0.02% agreement).

## Source: BPM Ch 20 "Magnetic Quantities and Units" — THE 931 MeV NUMBER

URL: `https://reciprocalsystem.org/books/bpom/20-magnetic-quantities-and-units`

**Verbatim quote**:

> "The natural unit of magnetic flux is the product of the natural unit of electric potential, **9.31146×10⁸ volts**, and the natural unit of time, 1.520655×10⁻¹⁶ seconds, and amounts to 1.41595×10⁻⁷ volt-sec, or webers."

This is the load-bearing find. The natural unit of electric potential in RS is `9.31146 × 10⁸ V`. One unit of charge `e` accelerated through one natural unit of potential = `e × 9.31146 × 10⁸ V` = **931.146 MeV**. This is *exactly* the unit conversion factor one uses for atomic mass: `1 u = 931.494 MeV/c²` (modern CODATA, agrees to ~0.04%).

**Implication for YM**: The factor 931.2 MeV in `Δ = ln(2π) × 931.2 MeV` is the natural unit of energy when divided by `e`. It's not a numerical coincidence — it's a postulate-level identity in RS.

**Magnetic dimensional table** (Table 31, BPM Ch 20):

| Quantity | SI | RS dimensions |
|---|---|---|
| Dipole moment | weber·m | t²/s |
| Flux | weber | t²/s² |
| Flux density | tesla | t²/s⁴ |
| Inductance | henry | t²/s³ |
| Permeability | henry/m | t²/s⁴ |

**Inductance verifies mass dimensionality**:

> "L = t/s² × t/s × t = t³/s³"

Inductance has the same RS dimensions as mass — Larson uses this as evidence that inductance is the EM analog of inertia.

## Source: SPU 1959 Ch 4 "Some Basic Relations"

URL: `https://reciprocalsystem.org/books/spu/04-some-basic-relations`

The original 1959 derivation table:

| Quantity | 1959 Value |
|---|---|
| Velocity | 2.9979 × 10¹⁰ cm/sec |
| Time | 0.1521 × 10⁻¹⁵ sec |
| Space | 0.4559 × 10⁻⁵ cm |
| Mass (sec³/cm³) | 3.7115 × 10⁻³² |
| Mass (g) | ≈ 0.5565 × 10⁻²⁴ |
| Acceleration | 1.97 × 10²⁶ cm/sec² |
| Force | 109.7 dynes |
| Energy | 5.0 × 10⁻⁴ ergs |

The mass dimensionality is asserted explicitly:

> "mass as the reciprocal of three-dimensional velocity" expressed as **t³/s³**

**Energy** in SPU 1959 = 5.0 × 10⁻⁴ erg in natural units = 3.12 × 10⁸ eV = 312 MeV. The factor of 3 between this and 931.146 MeV is the dimensional factor between one-dimensional speed (energy = t/s) and three-dimensional speed (mass = t³/s³). Specifically: at unit speed, `E_nat × c² = m_nat`, so `m_nat × c² = 3 × E_nat × c × c × c / (c × c) × ...` — Joe should rederive this in cold pass; the chain factor is dimensional, not magic.

## Source: Peret, "Subatomic Mass, Recalculated"

URL: `https://reciprocalsystem.org/paper/subatomic-mass-recalculated`

**Verbatim quote** of conversion factor used:

> "conversion factor of 931.49432 MeV/u"

Peret uses this as a *measured* CODATA value to convert RS natural mass units into experimental MeV/c² for comparison with PRD published particle masses. This sets the convention: `1 amu_RS = 931.49432 MeV/c² (CODATA-1994)`.

**Mass component values** (Peret, ratified from Larson):

| Component | Value |
|---|---|
| Primary mass | 1.000000000000 |
| Magnetic mass | 0.006392045455 |
| Electric mass (variants in paper) | ~0.0006–0.0009 |

The conversion baseline Peret proposes:

> "determine a conversion factor from natural mass units to unified atomic mass units based on an isotope-free, easily measured particle—the charged electron"

I.e., RS theory predicts the natural mass unit and Peret normalizes to the empirical electron mass; the resulting conversion factor reproduces 931.494 MeV/u to high precision. This is the chain Joe should re-trace cold.

## Source: Satz, "Identification of Cosmic Particles"

URL: `https://reciprocalsystem.org/paper/identification-of-cosmic-particles`

**Verbatim quote** (J/ψ derivation):

> "Rotational cosmic mass of c-H² = 1848 MeV/c²"
> "2(931.15 MeV/c²) = 1862.3 MeV/c²"
> "Total cosmic mass of c-H² = 1848 MeV/c² + 1862 MeV/c² = 3710 MeV/c²"
> Observed: 3695 MeV/c²

ψ(3105) → cosmic helium with two material gravitational charges:

> "rotational mass...is 3724.61/3 = 1242 MeV/c²"
> "two material gravitational charges mass 931.15 each making total mass 3104 MeV/c²"
> Observed: 3105 MeV/c²

These are Satz's worked examples of the **931 MeV/c² per gravitational charge** rule. This is precisely the mass-quantum unit Joe should expect to see appearing in any RS-flavored mass-gap result.

## Source: Boltzmann + temperature (BPM Ch 5 "Heat")

URL: `https://reciprocalsystem.org/books/bpom/05-heat`

Natural unit of temperature (gaseous): **7.20423 × 10¹² K**.
Effective unit accounting for vibrational structure: **3.5978 × 10⁹ K**.
Liquid/solid effective unit: **510.8 K** (one-dimensional reduction).
Boltzmann constant: **1.38044 × 10⁻¹⁶ erg/K** (Larson's value, agrees with CODATA 1.38065 × 10⁻¹⁶).
Natural unit of specific heat: **2.07066 × 10⁻¹⁶ erg/K** = (3/2) × R per natural mass unit.

`R = (2/3) × natural unit` of specific heat per natural mass — Larson, BPM Ch 5.

## Inter-regional ratio: 156.444

Recurs through BPM:

- **BPM Ch 4 Eq 4-9**: `P₀ = az / 938.67` (note: 938.67, not 156.44 — the 938.67 is `156.444 × 6.0046` or similar; need to re-derive cold).
- **BPM Ch 24 Eq 24-1**: `m_v = I × m_r² / 156.444`
- **SPU 1959 Eq 137**: `m_v = I × m_r² / 156.44` (same equation, old numbering)

The 156.444 is identified by Larson as the "inter-regional ratio" between time-region and outside-region motion. It plays the role of a partition function / RG factor between scales in RS.

For cold YM derivation, this ratio is candidate machinery for whatever appears alongside `ln(2π)`. Note `156.444 ≈ 128 × (1 + 2/9)` (NBM Ch 12) and `128 = 2⁷ = 2^(2⁺²+²+¹)` — suggestive of the 2³ × 8/(8/9-style three-dimensional displacement structure.

## Summary table — every numerical constant Larson asserts

| Constant | Value | Source |
|---|---|---|
| Unit speed (= c) | 2.997930 × 10¹⁰ cm/s | NBM Ch 13 |
| Unit space | 4.558816 × 10⁻⁶ cm | NBM Ch 13 |
| Unit time | 1.520655 × 10⁻¹⁶ s | NBM Ch 13 |
| Unit mass (s/cm) | 3.711381 × 10⁻³² sec³/cm³ | NBM Ch 13 |
| Unit mass (g) | 1.65979 × 10⁻²⁴ g per amu | NBM Ch 13 |
| Unit electric potential | 9.31146 × 10⁸ V | BPM Ch 20 |
| Unit magnetic flux | 1.41595 × 10⁻⁷ Wb | BPM Ch 20 |
| Unit force | 7.316889 × 10⁻⁶ sec/cm² | NBM Ch 13 |
| Inter-regional ratio | 156.444 | BPM Ch 24, SPU Ch 30 |
| Compressibility constant | 938.67 | BPM Ch 4 |
| Magnetic mass | 0.006392 amu | NBM Ch 13 |
| Electric mass (3D) | 0.000868 amu | NBM Ch 13 |
| Electric mass (2D) | 0.000579 amu | NBM Ch 13 |
| Charge mass | 0.000045 amu | NBM Ch 13 |
| G (cgs) | 6.67537 × 10⁻⁸ | NBM Ch 13 |
| Boltzmann's k | 1.38044 × 10⁻¹⁶ erg/K | BPM Ch 5 |
| Rydberg frequency | 6.576 × 10¹⁵ Hz | SPU Ch 4 |
| Natural T (gas) | 7.20423 × 10¹² K | BPM Ch 5 |
| Natural T (gas, eff.) | 3.5978 × 10⁹ K | BPM Ch 5 |
| Natural T (cond.) | 510.8 K | BPM Ch 5 |

## What's NOT in canonical Larson

**Joe needs to know**: The following do not appear in any pulled chapter of the four main volumes or the Peret/Satz papers:

- The factor `ln(2π)` explicitly. Natural log appears once (BPM Eq 1-1, ∫dt/t = ln(t)) but not with 2π.
- The factor `2π` in any mass/energy derivation. NBM Ch 12 confirms: "No logarithmic, exponential, or 2π/4π factors appear in this chapter." 2π appears in BPM Ch 21 in the Coulomb-form `μ₀I/(2πs)` magnetic field — purely geometric, not in mass derivations.
- The factor `4π` anywhere except as Coulomb-style geometric area.
- A Yang-Mills-style mass-gap argument. Larson rejects the strong nuclear force as "mythical" (his word, *Mythical Universe of Modern Astronomy* 1973).
- Birotation as 2π full period. The "two photons rotating around same center" idea (NBM Ch 10) is geometric but Larson never writes 2π for it.

The 931 MeV is canonical; the `ln(2π)` factor is not. Joe's cold derivation needs to either:
1. Show ln(2π) emerges from the BPM Eq 1-1 ln-integration combined with a 2π period from RS2 birotation (Peret), OR
2. Show ln(2π) ≈ 1.8379 has a closed-form derivation from inter-regional ratio mathematics that Larson stopped short of writing, OR
3. Conclude that ln(2π) is RS2-extension territory, originating with Peret/Nehru projective-geometry reformulation, not Larson.

The convergence test: if Joe's prior Opus 4.5 cold derivation got Δ = ln(2π) × 931.2, the 931.2 number comes straight from BPM Ch 20 unit definition (high confidence). The ln(2π) is the moving piece in any cold rederivation.
