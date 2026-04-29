# RS-Adjacent Research Corpus — Index

**Purpose**: Curated notes on frameworks, sources, and physics that border the Larson / Nehru / Peret canon but are not strictly inside it. Background reading for the Yang-Mills mass gap cold re-derivation (Hard Problems Paper #2). Headline target: `Δ = ln(2π) × 931.2 MeV ≈ 1.711 GeV`.

**Compiled**: 2026-04-28 by research agent for Joe Van Horn.
**Methodology rule observed**: No fabricated quotes, page numbers, or URLs. Inaccessible sources marked.

## Files

| File | Sub-area | Contents |
|------|----------|----------|
| `lebon-evolution-of-matter.md` | 1 — Le Bon (1907) | Chapters, mass-as-inertia equation `M = P/g`, etheric vortex ontology, dematerialization. Note: `p/v = M` does NOT appear verbatim in Le Bon's text — discussion below. |
| `bridge-frameworks.md` | 2 — Sister frameworks | Hestenes (geometric algebra / zitterbewegung), Rowlands (nilpotent quaternion), Hamilton & Tait, Larmor (rotational ether), Doran-Lasenby Clifford-QFT, Lerner (plasma cosmology), Carezani (autodynamics) |
| `unit-calibration.md` | 3 — RS unit calibration | Larson's natural units: t = 1.520655 × 10⁻¹⁶ s; s = 4.558816 × 10⁻⁶ cm; m = 1.65979 × 10⁻²⁴ g (= 1/Avogadro). Rydberg-frequency calibration. 931 MeV = 1 amu via E = mc². |
| `critics-and-objections.md` | 4 — Skeptical engagement | RationalWiki, Quora, Cosmoquest, Wright's plasma-cosmology critique. The "no math" charge, status of testable predictions. |
| `glueball-and-yang-mills.md` | 5 — Mainstream YM bridge | Clay Millennium statement, lattice QCD glueball spectrum (Morningstar-Peardon 1999, Chen 2006, Lucini 2014, Ochs 2013), 2024-2025 reviews. Best mainstream value for lightest scalar 0++ glueball: ~1.6-1.75 GeV. |

## Key findings at a glance

### Unit calibration (sub-area 3) — the most important practical finding

Larson's unit-of-mass calibration is **explicit and external**: the natural unit of inertial mass is the reciprocal of Avogadro's constant, `1.65979 × 10⁻²⁴ g`. This is exactly `1 / (6.02486 × 10²³)`. By construction, **one Larson natural unit of mass = 1 atomic mass unit (1 amu)**, and the conversion to MeV via `E = mc²` is therefore the *standard* textbook factor:

> 1 amu = 931.494 MeV/c² (CODATA value)
>
> Or, as quoted in RS literature: 1 amu = 931.15 MeV/c² (older value, still used in some Larsonian papers).

**Implication for the cold re-derivation**: the "931.2 MeV" in `Δ = ln(2π) × 931.2 MeV` is not a free parameter or RS-internal quantity — it is the same conversion factor every nuclear physics textbook uses. RS doesn't *derive* 931 MeV from postulates; it *inherits* it by choosing Avogadro's number as the calibrating constant. The interesting thing to derive in the YM cold pass is therefore the **dimensionless factor `ln(2π)`**, not the 931.2 MeV.

The natural unit of time `t₀ = 1.520655 × 10⁻¹⁶ s` is calibrated to the Rydberg fundamental frequency `R = 6.576115 × 10¹⁵ half-cycles/s` (Larson, *Nothing But Motion*, ch. 13).

### Mainstream lightest scalar glueball (sub-area 5)

The mainstream value for the lightest scalar 0++ glueball mass from pure-gauge SU(3) lattice QCD has converged in the band ~1.6–1.75 GeV across multiple independent calculations:

- **Morningstar & Peardon (1999)** [hep-lat/9901004]: m(0++) ≈ 1730 ± 50 ± 80 MeV (anisotropic lattice)
- **Chen et al. (2006)** [PRD 73, 014516]: m(0++) ≈ 1730(50)(80) MeV (anisotropic, larger lattice)
- **Lucini & collaborators (2014)** [arXiv:1401.1494]: m(0++) ≈ 1.6 GeV in pure Yang-Mills (with scale-setting caveats)
- **Wikipedia/PDG-style summary value**: 1730 ± 80 MeV/c² for 0++; 2400 ± 120 MeV/c² for 2++; 2590 ± 130 MeV/c² for 0-+

**Larson/Peret prediction** `Δ = ln(2π) × 931.2 MeV = 1.83788 × 931.2 = 1711 MeV` lands squarely **inside** the mainstream lattice 0++ band — closer to f₀(1710) than f₀(1500) and within 1σ of the Chen-2006 central value.

**Recent caveat (2024-2025 reviews)**: Some unquenched lattice studies (e.g., dynamical-quark calculations) find that *no* single resonance below 2 GeV is predominantly a glueball — the f₀(1500)/f₀(1710) candidates are mixtures. The "scalar glueball as small object" lattice result of 2024 (gravitational form factor, mass radius 0.263 fm) is the newest twist.

### Le Bon "p/v = M" status (sub-area 1)

The verbatim equation **`p/v = M`** does NOT appear in either *Evolution of Matter* (1907) or *Evolution of Forces* (1908). What Le Bon actually writes:

- "The mass which serves to characterize matter is only the measure of its inertia."
- The fundamental equation he gives is `M = P/g` (mass = weight / gravitational acceleration), in *Evolution of Forces* — i.e., dimensionless re-arrangement of Newton's second law.
- "All mass is mass of the ether, all momentum is momentum of the ether, and all kinetic energy is kinetic energy of the ether" — paraphrase from the Borderland Sciences review of *Evolution of Forces*.

So Peret's hand-written EE-dictionary entry "p/v = M" is **Peret's own re-derivation in the spirit of Le Bon**, not a Le Bon quote. Le Bon's actual claim is that mass is *purely* a measure of inertia, and that all such inertial properties are properties *of the ether*, not of the matter itself. This is consistent with — and ancestor to — the Larson/Peret position that mass is a derived quantity (in RS, motion is primary; mass = scalar component of compound-rotational motion).

### Critics & objections (sub-area 4)

The unanimous mainstream charge is "lacks mathematics." Quora and RationalWiki both make this case. Specific testable predictions:

- **Larson's quasar prediction** (1959) — exploding galaxy cores at superluminal speeds — predates mainstream identification of quasars but is *not* validated in its mechanism (mainstream attributes apparent superluminal motion to relativistic beaming on accretion disks).
- **The "no peer review" charge** is correct: Larson never published in a refereed physics journal.
- **No specific Larson prediction has been formally falsified** in the same sense that, e.g., autodynamics was disproved by particle accelerator data — because the predictions are typically too qualitative.

The cold-re-derivation methodology (Joe's: re-derive from postulates, then compare to mainstream) is therefore especially valuable here: it forces RS claims into a quantitative, falsifiable form that the original literature avoids.

## What's accessible vs paywalled

| Source | Status |
|--------|--------|
| Le Bon, *Evolution of Matter* (1907) | FREE — archive.org and rexresearch.com, public domain |
| Le Bon, *Evolution of Forces* (1908) | FREE — rexresearch.com, public domain |
| Larson, *Nothing But Motion* | FREE — archive.org, reciprocalsystem.org online edition |
| Larson, *Basic Properties of Matter* | FREE — nanopdf.com mirror; also reciprocalsystem.org |
| Bruce Peret RS2 papers (RS2-101..109) | FREE — reciprocalsystem.org PDFs |
| Morningstar & Peardon 1999 | FREE — arxiv.org/abs/hep-lat/9901004 |
| Chen et al. 2006 (PRD 73, 014516) | PARTIALLY PAYWALLED — APS abstract free; full PDF subscription |
| Lucini 2014 | FREE — arxiv.org/abs/1401.1494 |
| Crede & Meyer 2009 (Status of Glueballs) | FREE — arxiv.org/abs/0812.0600 |
| Ochs 2013 (Status of Glueballs) | FREE — arxiv.org/abs/1301.5183 |
| 2024-2025 reviews (arXiv 2502.02547, 2503.03821, 2504.09120, 2508.21821) | FREE — arXiv |
| PDG 2024 Review of Particle Physics | FREE — pdg.lbl.gov |
| Hestenes papers on zitterbewegung | FREE — davidhestenes.net, arxiv |
| Rowlands, *Zero to Infinity* (2007) | INFORMAL FULL PDF circulating; official Amazon/WorldScientific paywalled |
| Doran & Lasenby, *Geometric Algebra for Physicists* (2003) | INFORMAL PDF mirror; official Cambridge paywalled |
| RationalWiki, "Reciprocal Theory" / "Dewey Larson" | INACCESSIBLE: 503 errors during this research session — content recovered from earlier search-engine snapshots |

## See also (companion files in this corpus)

- `~/RS-Framework-Bridge/RS-research-corpus/larson/` — primary Larson sources
- `~/RS-Framework-Bridge/RS-research-corpus/nehru-peret/` — Nehru and Peret RS2 distillation
- `~/RS-Framework-Bridge/cold-derivation/` — Joe's cold re-derivation work (Riemann complete, YM next)
- `~/RS-Framework-Bridge/Riemann-Hypothesis/` — completed Phase 5.1 reference
- `~/RS-Framework-Bridge/RS2-foundations/` — Peret EE-dictionary (local-only, photo originals)
