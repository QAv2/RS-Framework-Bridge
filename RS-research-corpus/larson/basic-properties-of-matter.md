# *Basic Properties of Matter* (1988) — RS Volume II

Larson, D.B. *Basic Properties of Matter*, North Pacific Publishers, Portland OR, 1988.
URL base: `https://reciprocalsystem.org/books/bpom/`

Twenty-seven chapters. This is where Larson does the *quantitative* work — atomic distances, isotope masses, electrical and magnetic properties — using the natural-unit framework set up in NBM.

For YM rederivation, the load-bearing chapters are:
- **Ch 1** Solid Cohesion (Eq 1-1, the only ln integration in canonical RS)
- **Ch 4** Compressibility (the 938.67 constant)
- **Ch 14** The Basic Forces (RS rejection of strong force)
- **Ch 19** Magnetostatics (2D rotation = magnetic foundations)
- **Ch 20** Magnetic Quantities and Units (THE 9.31146 × 10⁸ V derivation)
- **Ch 24** Isotopes (the 156.444 inter-regional ratio in Eq 24-1)
- **Ch 27** Mass and Energy (mass = three-dim reciprocal speed, energy = one-dim)

## Ch 01 — Solid Cohesion

URL: `/01-solid-cohesion`

**This is the only place in canonical Larson with a logarithmic integration.** Critically load-bearing.

**Equation 1-1** (force from progression of natural reference system):
```
∫₁ᵗ (1/t) dt = ln(t)
```

**Equation 1-2** (two-atom force with Larson's notation):
```
F = ln(t_A) × ln(t_B)
```

**Equation 1-3** (extending to two-dimensional magnetic rotation):
```
F = ln²(t_A) × ln²(t_B)
```

**Equation 1-4** (effective magnetic rotation force, with the magnetic-mass coefficient 0.006392 from NBM Ch 13):
```
F_m = (0.006392)⁴ × s⁻⁴ × ln²(t_A) × ln²(t_B)
```

**Equation 1-5** (general equilibrium distance):
```
s₀ = 0.006392 × ln^(1/2)(t_A) × ln^(1/2)(t_B)
```

**Equation 1-6** (simplified, single element):
```
s₀ = 0.006392 × ln(t)
```

**Equation 1-7** (in Ångström units):
```
s₀ = 2.914 × ln(t) Å
```

Verbatim Larson:

> "The force due to the progression of the natural reference system [not gravitation] holds the solid aggregate together."

> The progression follows "1/1, 1/2, 1/3…1/n" requiring integration: "ln(n), the natural logarithm of n."

This is the entire RS theory of cohesion: an integration of `1/t` over time-region steps. The natural-log structure is the only place Larson introduces ln. **If Joe's `Δ = ln(2π) × 931.2` form is to be derivable in canonical RS, this is the equation that has to do the work.**

For YM mass-gap: the analog construction would be a confined system held together not by the progression-of-natural-reference (Eq 1-1) but by some sector-inversion analog. The mathematical machinery (`∫dt/t = ln(t)`) is the only RS-canonical source of logarithms; ln(2π) requires that the integration limit be 2π, which requires a 2π period to come from somewhere. **In canonical Larson it doesn't.** It must come from RS2 birotation period (Peret-Nehru extension).

## Ch 02 — Inter-Atomic Distances

URL: `/02-interatomic-distances`

References Eq 1-10 (full form not in this chapter's WebFetch window). Tables 2-6 give specific rotations and computed distances vs observed (Ångströms). Distance computation:

> "By applying this equation we find that the effective rotational force (ln t) for t = 2 is 0.693, which is less than the opposing space-time force 1.00."

Note: ln(2) = 0.69315, so this matches. The "opposing space-time force 1.00" is unit force in natural units. The transition from ln(2) < 1 to ln(t) > 1 at t = e ≈ 2.718 is when cohesive force overtakes the unit baseline.

Larson admits limits:

> "we are not yet in a position where we can determine specifically just what the inter-atomic distance will be for any given element under a given set of conditions."

## Ch 04 — Compressibility

URL: `/04-compressibility`

Pressure equations:
- **Eq 4-1**: `P = F/s²`
- **Eq 4-2**: `P = Fs/s³ = E/V`
- **Eq 4-3**: `E = PV` (gas law-equivalent for independent particles)
- **Eq 4-4**: `PV² = k` (time-region Boyle's-law equivalent)
- **Eq 4-5**: `V = k/P^(1/2)`

Internal pressure equations (the load-bearing ones for the 938.67 constant):
- **Eq 4-9**: `P₀ = az / 938.67`
- **Eq 4-10**: `P₀ = azy / (938.67·V)`
- **Eq 4-11**: `P₀ = azy / (936.67·s₀³)` — note: WebFetch transcribes "936.67" here; likely a typo in the source for 938.67. Joe should verify against print.
- **Eq 4-12**: `P₀ = 17109·azy / s₀³` kg/cm²

Compressibility rate:
- **Eq 4-13**: `(1/V₀)(dV/dP) = P₀^(1/2) / [2(P₀ + P)^(3/2)]`
- **Eq 4-14**: `(dV/dP)|_(P=0) = 1 / (2P₀)`

The **938.67** constant relates mass-pressure to atomic structure. It's distinct from but related to the 931.146 of Ch 20. Plausibly: `938.67 = 938.27 (proton mass MeV) × ε` or `938.67 = 931.146 × 1.00808` — a slight modification accounting for the proton-vs-amu difference. Joe should pin this in cold derivation.

## Ch 05 — Heat

See `unit-conversion.md`. Boltzmann constant, natural temperature units. Energy = T in natural units (Eq 5-2: `E = T`). Kinetic energy in time-region: `E = kT²` (Eq 5-5). Heat: `H = T²/n³` (Eq 5-6). Specific heat: `dH/dT = 2T/n³ + I` (Eq 5-8). Radiation: `E_rad = kT⁴` (Eq 5-4).

The `T²/n³` and `T⁴` patterns echo Stefan-Boltzmann for radiation (modern: σT⁴) and equipartition (modern: cv ~ R per degree of freedom).

## Ch 14 — The Basic Forces

URL: `/14-the-basic-forces`

**RS recognizes only three basic forces**:

1. Gravitational — "inward translational manifestation"
2. Electrostatic — "one-dimensional motions of an oscillating character"
3. Magnetostatic — "follows the equation for the gravitational force"

> "The basic forces were identified as the force aspects of these basic motions."

Generalized force equation:
```
F = k·X·Y·X' / d²
```

Where X is any distributed scalar motion of dimensions `(t/s)ⁿ`, Y is a dimensional reduction factor, k is a unit-system numerical constant.

> "Force in general is the product of mass and acceleration."

**No strong nuclear force.** WebFetch confirms: "The text does not address strong nuclear force or discuss why nuclear binding might not require a separate force mechanism in the Reciprocal System framework."

This is RS's foundational position: the strong force is *unnecessary* because atoms cohere by rotational structure (BPM Ch 1). The Yang-Mills mass gap Δ is not the binding scale of nucleons in RS; it's something else. In RS terms, Δ is most likely the **3D-rotational-displacement threshold below which mass = 0** (NBM Ch 11), expressed as a mass-equivalent in MeV.

## Ch 19 — Magnetostatics

URL: `/19-magnetostatics`

> "magnetic charges, which are two-dimensional rotational vibrations acting in opposition to two-dimensional rotations"

> "The two-dimensionality is the key to understanding the magnetic relations"

Force law:
> "The magnetic force equation, the expression for the force between two magnetic charges, is identical with the Coulomb equation, except for the factor t/s introduced by the second scalar dimension"

`F = MM' / d²`, with dimensions `t/s²`.

Magnetic charges are bipolar (north + south poles) by construction:

> Each atom develops "two poles—north and south—because the two-dimensional rotation divides the atomic structure into opposing halves with opposite rotational directions."

This is the geometric origin of magnetic dipoles in RS. Note: Larson does not give a numeric value for "natural unit of magnetic charge"; Ch 20 gives natural unit of magnetic flux instead.

## Ch 20 — Magnetic Quantities and Units

URL: `/20-magnetic-quantities-and-units`

Already pulled in `unit-conversion.md`. Key facts:

**Verbatim** (load-bearing):

> "The natural unit of magnetic flux is the product of the natural unit of electric potential, **9.31146×10⁸ volts**, and the natural unit of time, 1.520655×10⁻¹⁶ seconds, and amounts to 1.41595×10⁻⁷ volt-sec, or webers."

`9.31146 × 10⁸ V` is THE natural unit of electric potential. 1 e × 9.31146 × 10⁸ V = **931.146 MeV** = 1 amu (modern). The 931 in Joe's `Δ = ln(2π) × 931.2 MeV` is this natural unit.

Inductance verifies mass dimensionality:
```
L = t/s² × t/s × t = t³/s³
```

Force comparison:
```
F = ma = m·dv/dt = m·d²s/dt²
F = L·dI/dt = L·d²q/dt²
```

Same dimensional structure — inductance is the EM analog of inertia in RS.

Table 31 (Principal Magnetic Quantities) — see `unit-conversion.md`.

## Ch 21 — Electromagnetism

URL: `/21-electromagnetism`

Field rejection:

> "The field is not 'something physical.' It is merely a mathematical consequence of the inability of the conventional reference system to represent scalar motion in its true character."

Maxwell equations are not given in canonical form in this chapter. The 2π appears only in geometric Coulomb-style formulas:

```
B = μ₀I / (2πs)
F = μ₀ × I_A × I_B × l / (2πs)
```

Conversions:
- emu/esu ratio: `3 × 10¹⁰`
- Electric unit: `4.80287 × 10⁻¹⁰ esu = 1.60206 × 10⁻²⁰ emu`
- Resistance: `R = mass/time = t²/s³`
- Permeability: `μ = mass/space = t³/s⁴`

**No non-Abelian generalization, no gauge group structure.** This is by design — Larson treats EM as scalar-motion phenomenology, not field theory.

## Ch 24 — Isotopes

URL: `/24-isotopes`

**Equation 24-1** (the canonical isotope-mass equation):
```
m_v = I × m_r² / 156.444
```

Where:
- `m_v` = vibrational mass (in gravitational charge units)
- `I` = magnetic ionization level
- `m_r` = rotational mass (atomic number units)
- `156.444` = inter-regional ratio (time-region atoms)

Atomic mass formula:
```
M = 2Z + G
```

Z = atomic number, G = gravitational charge units.

Worked example: lead. Unit ionization basis gives `m_v = 43`. Adding the 164 atomic-weight units from `2Z + G` (Z = 82, G = 0): theoretical = 207. Observed = 207.2. Agreement: 0.1%.

> "the atomic weight thus defined is a reflection of the local neutrino concentration, the magnetic temperature"

## Ch 25 — Radioactivity

URL: `/25-radioactivity`

Mass-loss in a-b-c notation:
```
O¹⁶ → C¹² + He⁴
16-0 → 12-0 + 4-0
```

Half-lives mentioned: U²³⁸ at 4.5 × 10⁹ years; Pb²¹⁰ at 22 years. No systematic table.

No 4-x or 8-x inversion in this chapter (those are in NBM Ch 15).

## Ch 26 — Atom Building

URL: `/26-atom-building`

> "The principal product of the decay of cosmic atoms is the massless neutron"

Combination chain: massless neutron + neutrino → proton (`M 1-1-(1)`, mass = 1 atomic weight unit).

> "Magnetic ionization is therefore an atom-building process of such broad scope that it is clearly the predominant means of accomplishing the formation of the heavier elements."

> "The transition from the massless state cannot be reversed in the material environment, as there is no available process for going directly from rotation to translation."

Three-stage transition: cosmic → neutral → material. The "neutral" state is unit-rest (no scalar progression).

## Ch 27 — Mass and Energy

URL: `/27-mass-and-energy`

This is the philosophical close of BPM. Verbatim:

> "In the universe of motion defined by that system of theory, mass and energy are both reciprocal speeds, differing only in dimensions, mass being three-dimensional, while energy is one-dimensional."

> "Unit energy is therefore the product of unit mass and the second power of unit speed, the speed of light."

This is RS's E = mc² statement. **Mass = three-dim reciprocal speed; energy = one-dim reciprocal speed.** The factor c² is dimensional, falling out of the dimensional reduction `t³/s³ × s²/t² = t/s`.

WebFetch summary: chapter is "primarily interpretive and historical in nature." No 2π or ln derivations here.

## Net assessment for YM rederivation

BPM provides:
1. The **only ln integration** in canonical RS (Eq 1-1) — the natural-log machinery for any "Δ = ln(?) × natural unit" form must trace to here.
2. The **931 MeV natural unit** (Ch 20) — definitely-canonical, settles the 931.2 number in Joe's headline.
3. The **156.444 inter-regional ratio** (Ch 24) — likely candidate for the dimensionless coefficient that multiplies the natural unit.
4. The **mass = t³/s³** dimensional identity (Ch 12 from NBM, restated in Ch 27 of BPM) — the dimensional structure.
5. The **rejection of strong nuclear force** (Ch 14) — Joe's YM rederivation must give a non-strong-force account of the mass gap.

What BPM does NOT provide:
- 2π in mass derivations
- ln(2π) anywhere
- Yang-Mills-style derivation of confined-state masses
- A general theorem connecting rotational displacement to a continuum-limit mass gap
