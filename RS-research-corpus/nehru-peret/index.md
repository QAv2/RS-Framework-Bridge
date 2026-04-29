# Nehru / Peret Research Corpus — Index

**Compiled**: 2026-04-28 for Yang-Mills Hard Problem #2 (RS2 cold rederivation).
**Working directory**: `~/RS-Framework-Bridge/RS-research-corpus/nehru-peret/`
**Source-of-truth site**: [reciprocalsystem.org](https://reciprocalsystem.org/) (the merged successor to old `reciprocalsystem.com`).
**Owning society**: International Society of Unified Science / Reciprocal System Research Society (ISUS, Inc.), Salt Lake City UT — non-profit. The journal is **Reciprocity** (quarterly).

This corpus extends Joe's existing `RS2-foundations/` (Peret 101–109 + EE-dictionary + Nehru spin photo) by pulling everything else from the canonical archive that bears on **strong nuclear force, mass derivation, hadron spectrum, atomic-zone partition, and unit-conversion to MeV**.

---

## Files in this corpus

| File | Purpose |
|---|---|
| `index.md` | This file — map of the corpus, accessibility log, methodology notes |
| `nehru-bibliography.md` | Complete list of 49 Nehru papers with PDF URLs and accessibility status |
| `nehru-key-papers.md` | Distilled content of the most load-bearing Nehru papers (high-energy physics, time region QM, wave mechanics, non-locality, inter-regional ratio, internal ionization, neutron lifetime, Planck's constant) |
| `peret-extended.md` | Peret papers beyond 101–109 (Subatomic Mass Recalculated, Tao of Larson, Periodicals Master Index, etc.) |
| `quaternion-organon.md` | What "Quaternion Organon" actually is — important correction to prior memory |
| `reciprocity-journal-index.md` | Master Index of the Reciprocity journal — issue-by-issue inventory of every Nehru/Peret/Satz paper that bears on mass / nuclear physics |
| `nuclear-and-mass.md` | Focused extract: every quantitative claim about strong nuclear force, mass spectrum, glueball-equivalents, MeV calibration, found in the Nehru/Peret/Satz/Vijaya corpus |
| `_pdfs/` | 18 source PDFs downloaded from reciprocalsystem.org (canonical) |
| `_text/` | Layout-preserved text extractions of the same PDFs |

---

## Author roster (canonical, from the 184-paper alphabetical list)

The paper authors active in *Reciprocity* (1971–~2010s):

- **Dewey B. Larson** — founder of the Reciprocal System (~30+ papers in the journal beyond his books)
- **Prof. KVK Nehru, Ph.D.** — JNTU India, ISUS member; the most prolific *extension* author. **49 papers** on the canonical site
- **Dr. Bruce Peret** — current Director of Research, ISUS; **22 papers** on the site, including the entire **RS2 Tutorial Series 101–109**
- **Dr. Ronald Satz** — author of mass-derivation extensions (cohesive energy, dissociation energy, equation of state, magnetic charge, time region particle dynamics, cosmic-particle identification). Resigned from ISUS at some point (per master index)
- **Gopi Krishna Vijaya** — derives lepton magnetic moments in RS2 (e/μ/τ g-factors to 12 decimal places without QED loops)
- **Frank H. Meyer** — early Reciprocity editor; absolute magnitudes
- **David Halprin** — atomic-number equation, ionization levels, four magnetisms
- **Thomas Kirk** — photon model, sub-atomic particle array, motion fundamentals
- **George Hamner** — author of *Quaternion Organon* novel (2001) — see correction note in `quaternion-organon.md`
- **Daniel Phoenix III** — the "--daniel papers"; metaphysics-leaning

**No "Antoine St. James"** appears in the canonical inventory. If Joe's prior notes have that name, it may be a hallucination from earlier sessions.

---

## Accessibility status (summary)

- **Fully extracted to text in this corpus** (18 papers): high-energy-physics, quantum-mechanics-time-region, wave-mechanics, non-locality, inter-regional-ratio, internal-ionization, lifetime-neutron, planck-constant, on-the-nature-of-rotation-and-birotation, scalar-rotation, spin, photon-birotation, superconductivity, electric-ionization, lifetime-c-argon-muon, lifetimes-c-atom, periodicals-master-index, vijaya-magnetic-moments-leptons, subatomic-mass-recalc, tao-of-larson.
- **Listed but not pulled** (the remaining ~28 Nehru papers + most Satz papers): URLs known, can be fetched any time. Not pulled this session because they are second-priority to the mass / nuclear axis.
- **Partial — HTML-rendered on site, no PDF**: a handful of older Satz papers (e.g. *Cosmic Rays and Elementary Particles*, *Cohesive Energy*) live as page text rather than PDF. Summaries captured by WebFetch are in `nuclear-and-mass.md`.
- **Inaccessible**: nothing critical. The KVK Nehru biography page (`/kvk-nehru`) returns 404; biographical data is partial. The old `reciprocalsystem.com` is fully merged into `.org` per the site's own merge note.

---

## Methodology used

1. Pulled the two author-index pages: [`/papers-for-author/4`](https://reciprocalsystem.org/papers-for-author/4) for Nehru, [`/papers-for-author/3`](https://reciprocalsystem.org/papers-for-author/3) for Peret. These are the canonical lists.
2. Cross-checked against the [alphabetical paper list](https://reciprocalsystem.org/papers) (184 entries) — gives author attribution for cases where the page-3/4 listings are short.
3. Pulled the **Reciprocity Master Index** PDF (Peret) — gives volume/issue/page coordinates for every paper in the journal back to 1971, which lets us cite them properly (`Reciprocity XXVI, № 2, page 7` etc.).
4. Downloaded source PDFs with `curl` (WebFetch can't parse PDF binary; it does cache them, but direct curl is cleaner).
5. Converted with `pdftotext -layout` to preserve tabular data (Larson's mass tables, Nehru's hadron resonance tables).
6. Where the reciprocalsystem.org page renders the paper *as HTML* rather than serving PDF, used WebFetch to extract the text (Satz papers).

## Source primacy

Per Joe's standing instruction, only canonical/upstream sources used:
- ✅ `reciprocalsystem.org` — the canonical ISUS site
- ✅ NewSouth Books / publisher metadata for the Hamner *Quaternion Organon* novel
- ❌ No `isunivers.org` or `rstheory.org` content was actually pulled — both came up as referenced names but the canonical material on `reciprocalsystem.org` covers their content (per ISUS's site-merge statement). The forum discussions there are secondary sources by Joe's standing rule.
- ❌ No blog/Medium/Substack chatter pulled.
