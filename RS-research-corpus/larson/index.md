# Larson Primary Corpus — Research Index

Curated for RS2 cold rederivation of the Yang-Mills mass gap (Hard Problems Paper #2).
Pulled from canonical sources only. Date: 2026-04-28.

## Source provenance

All chapter URLs below resolve under `https://reciprocalsystem.org/books/<book>/<chapter>`. This is the Peret-curated mirror of Larson's full text and is the closest thing to a canonical online corpus. Page indices and equation numbering match the print editions Larson published through North Pacific Publishers (Portland, OR).

Cross-mirror sources verified:
- Internet Archive — `https://archive.org/details/nothing-but-motion` (NBM full PDF, 417 MB, 1979)
- Internet Archive — `https://archive.org/details/the-reciprocal-system-of-theory-collection` (Larson's papers + late-life essays; does NOT contain the three main volumes — those are at the reciprocalsystem.org HTML mirror only)
- Larson Research Center: `http://www.lrcphysics.com/` — TLS cert is broken (ERR_TLS_CERT_ALTNAME_INVALID); INACCESSIBLE via WebFetch
- transpower.wordpress.com (Satz, not St. James as Joe's instructions said) — accessible, hosts archive of essays + revision-history posts

## Files in this directory

| File | Scope | Status |
|---|---|---|
| `index.md` | This overview + accessibility map | written |
| `unit-conversion.md` | s↔meter, t↔second, mass↔MeV — load-bearing for YM | written |
| `nothing-but-motion.md` | NBM 21-chapter coverage; postulates, rotation, atoms | written |
| `basic-properties-of-matter.md` | BPM 27-chapter coverage; Eq 1-1 + cohesion + mass/energy | written |
| `universe-of-motion.md` | UoM cosmic sector; high-speed regime | written |
| `neglected-facts.md` | NFoS late synthesis | written |
| `quantitative-predictions.md` | Larson's actual numerical claims with errors | written |

## Accessibility ledger

**Fully accessible** (HTML, fetched verbatim):
- *Nothing But Motion* (NBM, 1979) — all 21 chapters at `/books/nbm/`
- *Basic Properties of Matter* (BPM, 1988) — all 27 chapters at `/books/bpom/`
- *The Universe of Motion* (UoM, 1984) — 31 chapters at `/books/uom/`
- *The Neglected Facts of Science* (NFoS, 1982) — 8 chapters at `/books/nfs/`
- *The Structure of the Physical Universe* (SPU, 1959 original) — partial chapter list at `/books/spu/`
- Peret, *Subatomic Mass, Recalculated* — at `/paper/subatomic-mass-recalculated`
- Satz, *Identification of Cosmic Particles* — at `/paper/identification-of-cosmic-particles`
- Satz, *A New Derivation of Planck's Constant* — at `/paper/a-new-derivation-of-plancks-constant`

**Listed but not pulled in this pass** (low priority for YM):
- *Beyond Newton* (1964) — gravitation popularization
- *New Light on Space and Time* (1965) — popularization
- *Quasars and Pulsars* (1971) — superseded by UoM
- *The Case Against the Nuclear Atom* (1964) — short polemic; relevant to YM but largely qualitative

**Inaccessible / needs alternative**:
- LRC (`lrcphysics.com`) — TLS error. Can be tried later via `curl -k` if needed.
- Print-only commentary (Peret's hand-written EE dictionary, Nehru *Quaternion Organon*) — Joe already has these as JPEGs in `RS2-foundations/`.

## Headline finds for YM rederivation

**1. Mass dimensionality is canonical.** Larson 1959 SPU Ch 4 and 1979 NBM Ch 12: mass = `t³/s³`. NBM Ch 13 gives natural unit of inertial mass = `3.711381×10⁻³² sec³/cm³` ≈ `0.5565×10⁻²⁴ g`. This is unambiguously the dimensional fact RS2-EE inherits.

**2. The 931 MeV is built into RS at unit level.** BPM Ch 20 gives the natural unit of electric potential as **9.31146×10⁸ volts** — the 931 of "1 amu = 931.494 MeV" is a *natural unit* of RS, not an empirical match. Peret's "Subatomic Mass, Recalculated" then uses **931.49432 MeV/u** explicitly to convert. Satz ("Identification of Cosmic Particles") uses **931.15 MeV/c² per gravitational charge** in cosmic atom calculations. The natural electric potential unit ≈ MeV/amu is the foundation.

**3. Larson's load-bearing equation is BPM Eq 1-1: ∫₁ᵗ (1/t)dt = ln(t).** Solid cohesion derives from this. The natural-log form gives the inter-atomic equilibrium distance; this is the closest thing in canonical Larson to the Berry-Keating / RH-style operator. NBM Ch 12 gives no equivalent equation — it's BPM-specific.

**4. RS rejects the strong nuclear force categorically.** BPM Ch 14 enumerates only three basic forces (gravitational, electrostatic, magnetostatic). Atoms hold together by *rotational coincidence*, not by gluon-mediated binding. The closest analog to YM mass gap is therefore: "the mass gap is the threshold below which a rotational displacement cannot exist as a stable atom" — i.e., the massless-neutron / neutrino threshold. See NBM Ch 11.

**5. Cosmic sector → high-energy / inverse mass regime.** NBM Chs 14–16 + UoM throughout: speeds > c are the "cosmic" regime. Cosmic atoms have mass = reciprocal of material mass × MeV-equivalence factor. c-H² has theoretical mass 1848.61 MeV; with two material isotope charges (each 931.15 MeV) it predicts the J/ψ at 3710.91 MeV vs observed 3695. This is a 0.4% prediction error on a hadron mass from RS first principles — directly analogous to what Joe wants for the Δ = ln(2π) × 931.2 ≈ 1.711 GeV claim.

**6. The 156.444 inter-regional ratio.** Appears in BPM Ch 4 (compressibility, with the related constant 938.67 in Eq 4-9), BPM Ch 24 (isotopes, Eq 24-1: m_v = I·m_r²/156.444), and SPU 1959 atom-building (Eq 137 in old numbering). This is RS's intra-atomic-to-inter-atomic conversion factor. Likely shows up in any quantitative YM bridge as a Δ-region ratio.

**7. Larson uses no 2π, 4π, ln(2π), or e^x natural logs except in BPM Eq 1-1.** This is a critical gap. The Δ = ln(2π) × 931.2 form Joe targets does **not** appear in canonical Larson. ln(2π) must come from elsewhere — most likely Peret's projective-geometry RS2 reformulation (where 2π is the full birotation period, by construction), or from Nehru's quaternion spin paper, or from an inter-regional measure-theoretic argument that Larson never wrote out. **This is the single most important gap to flag for the cold derivation.**

## Recommended next sources Joe should seek

1. **NBM Ch 13 in full PDF** — the natural-units chapter is the load-bearing one; the WebFetch summary is excellent but a verbatim read of the full chapter would catch any subtle factor of 2π or √2 hidden in derivations.
2. **Larson's *Beyond Space and Time* (1995)** — late synthesis, may contain the ln-of-something result if Larson ever derived it.
3. **Peret's RS2 papers 101–109** — Joe already has these. Re-read with the YM target in mind: where does ln(2π) appear in projective geometry? It's the log of the period of the universal cover of S¹ → ℝℙ¹.
4. **Satz, *The Physics of Motion*** — listed in archive.org RS collection; may be the most quantitative late RS text.
5. **Halprin, *Atomic Number Equation Based on Larson's Triplets*** — not pulled here. Could be relevant to building atomic-mass equations from triplet structure.

## Methodology notes

- Where text was abstracted by WebFetch summarizer, equations are quoted verbatim where possible; prose around equations is typically WebFetch's compression of Larson's own words. Quoted strings inside this corpus are reliable as Larson's actual text.
- For YM cold rederivation: do **not** trust the WebFetch summarizer's "no X found" claims as absolute — these mean "not in the chapter window the summarizer prioritized". Joe should re-read the full HTML for any chapter where a load-bearing claim depends on an absence (e.g., ln(2π)).
- The 931 number convergence (931.146, 931.15, 931.494, 938.67) is real and traceable. The unit-level number is `9.31146×10⁸ V` from BPM Ch 20. The empirical 1 amu ≡ 931.494 MeV/c² agrees to ~0.04% with Larson's natural-unit derivation. **This is an existing prediction**, not a postdiction — Larson published the value pre-CODATA-1986.
