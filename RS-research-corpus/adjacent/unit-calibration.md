# RS Unit Calibration — From Postulates to MeV

**Sub-area 3**. The single most important practical question for the YM cold derivation: where does the 931 MeV in `Δ = ln(2π) × 931.2 MeV` actually come from in RS?

## Headline answer

**931 MeV is not a derived quantity in RS — it is the standard textbook conversion factor `1 amu = 931.494 MeV/c²`** (CODATA), inherited by RS via Larson's choice to make his natural unit of inertial mass equal to one atomic mass unit (1 amu). The choice is made by anchoring on Avogadro's constant.

The non-trivial derivations in RS are:

1. The **natural unit of time** `t₀ = 1.520655 × 10⁻¹⁶ s`, calibrated to the Rydberg fundamental frequency.
2. The **natural unit of space** `s₀ = 4.558816 × 10⁻⁶ cm`, derived from `t₀ × c`.
3. The **natural unit of inertial mass** `m₀ = 1.65979 × 10⁻²⁴ g`, calibrated to Avogadro's constant.

Once you have `m₀ = 1 amu`, the conversion `m₀ × c² = 931.494 MeV` follows from Einstein, not from RS-internal logic.

## Sources

| URL | Status | What's there |
|-----|--------|--------------|
| https://reciprocalsystem.org/books/nbm/13-physical-constants | FREE | Larson, *Nothing But Motion*, ch. 13 — the canonical Larson unit-calibration chapter |
| https://reciprocalsystem.org/welcome | FREE | The two RS postulates verbatim |
| https://reciprocalsystem.org/paper/identification-of-cosmic-particles | FREE | Cosmic-particle paper with explicit `1 amu = 931.15 MeV/c²` calibration |
| https://reciprocalsystem.org/papers (Satz, "Further Mathematics of the Reciprocal System") | FREE | Mathematical formalism extension |
| https://reciprocalsystem.org/papers (Peret, "Subatomic Mass, Recalculated") | FREE | Updated 1995 mass values for Larson particles |
| https://archive.org/details/nothing-but-motion | FREE | Full book PDF |
| http://www.lrcphysics.com/rst/ | INACCESSIBLE: TLS certificate name mismatch (lrcphysics.com cert is invalid) — content recovered via search-engine summaries |

## The Reciprocal System postulates (Larson, 1959)

Verbatim from reciprocalsystem.org/welcome:

> **Postulate 1 (Composition):** "The physical universe is composed of one component, motion, existing in three dimensions, in discrete units and with two, reciprocal aspects, space and time."
>
> **Postulate 2 (Mathematical structure):** "The physical universe conforms to the relations of ordinary, commutative mathematics, its primary magnitudes are absolute and its geometry is Euclidean."

So before any calibration: motion is primary, space and time are reciprocal aspects, both are quantized in "natural units," geometry is Euclidean, and arithmetic is commutative.

## Larson's calibration chain (Nothing But Motion, ch. 13)

### Step 1 — Identify a known physical constant as the calibration anchor

Larson chose the **Rydberg fundamental frequency**:
> R_∞ × c = 6.576115 × 10¹⁵ half-cycles per second.

This is the frequency associated with the Rydberg constant, the most precisely measured physical constant available to him in 1959.

### Step 2 — The natural unit of time is the reciprocal of this frequency

> **t₀ = 1 / (6.576115 × 10¹⁵) ≈ 1.520655 × 10⁻¹⁶ s**

(verbatim from Larson, *Nothing But Motion*, ch. 13, as recovered from reciprocalsystem.org online edition)

### Step 3 — The natural unit of space follows from c (Postulate 1: motion = space/time)

Because the speed of light c = 2.997930 × 10¹⁰ cm/s is "1 unit space / 1 unit time" in natural units (Larson's identification — the fundamental motion is at unity-velocity in natural units),

> **s₀ = c × t₀ = 4.558816 × 10⁻⁶ cm**

### Step 4 — The natural unit of inertial mass is the reciprocal of Avogadro's number

This is the step where Larson "smuggles in" the conventional 1 amu unit. Verbatim from ch. 13:

> "The reciprocal of [Avogadro's] number, 1.65979 × 10⁻²⁴, in grams, is therefore the mass equivalent of unit atomic weight."

That is:

> **m₀ = 1 / N_A = 1 / (6.02486 × 10²³) ≈ 1.65979 × 10⁻²⁴ g = 1 atomic mass unit**

(The Avogadro value Larson used, 6.02486 × 10²³, is slightly different from the modern CODATA value 6.02214 × 10²³ — this is a 1959 calibration that has not been refreshed in the published RS literature.)

### Step 5 — Convert to MeV via E = mc²

This step is **not** in *Nothing But Motion* — Larson uses ergs throughout. But in the cosmic-particle paper (Satz, on reciprocalsystem.org), the conversion is explicit:

> "1 atomic mass unit = 931.15 MeV/c²"

> "Mass of each gravitational charge is one atomic weight unit = 931.15 MeV"

The 931.15 MeV is the older value; modern CODATA is 931.494 MeV. The difference is ~0.04% and not load-bearing for the YM derivation at the precision Joe is targeting.

## What Larson actually predicts (and why)

Once the unit chain `t₀ ↔ s₀ ↔ m₀` is locked, Larson's claim is that all physical constants reduce to combinations of these natural units (and pure numbers), with the constants disappearing when expressed in natural units. The sub-atomic and atomic mass calculations in *Basic Properties of Matter* (1988) and the recalculations by Peret (1995, on reciprocalsystem.org) work in natural units, then convert back via `m₀ = 1 amu = 931.15 MeV/c²` for comparison to particle physics data.

The cosmic-particle paper gives the formula:

> Cosmic mass = 3724.61 / m, in MeV/c² (where m is atomic weight)

and uses 1 amu = 931.15 MeV/c² as the conversion factor throughout Table 1.

## Comparison to mainstream natural unit systems

For perspective, mainstream particle physics natural units:

- **Planck units**: c = ℏ = G = 1, derived from G, ℏ, c.
- **Atomic units (Hartree)**: ℏ = m_e = e = 1, derived from electron mass, Planck's constant, electron charge.
- **Atomic units (Rydberg)**: similar but with energy unit half of Hartree.
- **High-energy physics units**: ℏ = c = 1, energies in eV/MeV/GeV.

Larson's natural units are **Rydberg-based** (anchored on the Rydberg frequency, not on G or m_e). This puts them in the same family as atomic units.

## Implications for the YM cold derivation

### The 931.2 MeV in `Δ = ln(2π) × 931.2 MeV`

This is just the conversion factor 1 amu → MeV. Larson's natural unit of mass is *defined* to equal 1 amu, so any RS prediction expressed in "natural mass units" trivially converts to MeV by multiplying by 931.494.

**Therefore, the load-bearing physics of the cold derivation is the dimensionless factor `ln(2π) ≈ 1.83788`.** That number must come out of the RS axioms (or RS2 quaternion mechanics) without invoking experimental input. The 931.2 is just the textbook conversion factor.

### What `ln(2π)` could mean in RS terms

Some structurally suggestive identities to look at when cold-deriving:

- `ln(2π) = ln(2) + ln(π)` — a sum of two scale factors
- `2π` is the natural-unit period for a full rotation (a full quaternion cycle in the (i,j) plane)
- `ln(2π)` appears in the **Stirling expansion** `ln Γ(n+1) ≈ n ln n − n + (1/2) ln(2π n)` — the half-log-2π is the leading sub-asymptotic correction to factorials
- `(1/2) ln(2π)` appears in the **functional equation of the Riemann zeta function** — the same ξ(s) symmetry that drove the Riemann derivation
- `ln(2π)` appears as `arg Γ(1/2 + iT)` integrated over T in the explicit-formula treatment of zeta zeros (Riemann-von Mangoldt term)

The Riemann derivation that Joe completed in Phase 5.1 already produced `ln(2π)` as the leading-order term in the zero-counting function `N(T) = (T/2π) ln(T/2πe) + ...` — the same `2π` appears. So one promising attack on Paper #2: derive `ln(2π)` as the same Hilbert-Pólya–Berry-Keating gap structure, but applied to the YM Hamiltonian rather than the Riemann zeta operator.

### Independent check: the value 1711 MeV

> ln(2π) × 931.494 MeV = 1.837877 × 931.494 = 1711.36 MeV

versus the mainstream lattice 0++ glueball band 1.6 – 1.75 GeV. 1711 lands almost exactly on the f₀(1710) candidate and within 1σ of the Chen-2006 lattice central value 1730 MeV.

This is what makes Paper #2 worth doing: the prediction is sharp, falsifiable, and currently within the experimental band.

## What's missing from the public RS literature

1. **No published RS derivation of `ln(2π)` as a fundamental factor.** The cold pass would be original work.
2. **No published RS treatment of the Yang-Mills mass gap as such.** RS atomic/sub-atomic masses are tabulated, but the strong-force binding scale Δ is not isolated as a quantity in the RS literature.
3. **No update of Avogadro's constant to modern CODATA values.** The published RS literature still uses 6.02486 × 10²³. This is a minor issue (~0.04% drift).
4. **No published RS-to-QCD bridge.** Larson treats nucleons as compound rotations but does not engage with QCD or color confinement.

These four gaps are precisely what Paper #2 needs to fill.

## Recommended derivation chain for Paper #2

1. State the two RS postulates.
2. Derive the natural unit chain: `t₀` from Rydberg, `s₀` from `c × t₀`, `m₀` from Avogadro.
3. Argue (in the spirit of Phase 5.1 Riemann) that the *gap* in any RS Hamiltonian is set by a Hilbert-Pólya-style operator with a symmetry-fixed leading coefficient.
4. Show that for the Yang-Mills sector — i.e., the strong/color sector in RS2's three-quaternion-axes language — the leading coefficient is `ln(2π)`.
5. Convert to MeV via the standard `m₀ × c² = 931.494 MeV` conversion.
6. Compare to mainstream lattice glueball values (Morningstar-Peardon, Chen, Lucini).

This mirrors the Riemann §7 numerical-confirmation structure: cold derivation first, mainstream comparison second.
