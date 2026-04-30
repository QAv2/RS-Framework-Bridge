# Navier-Stokes Regularity (Hard Problem #1) — RS2 Bridge

This folder contains three artifacts on the Navier-Stokes existence-and-smoothness problem (Clay Mathematics Institute Millennium Prize Problem) under the Reciprocal System (RS2) framework, ordered chronologically.

## Reading order

For someone coming to this folder fresh:

1. **`prior-art/rs2_navier_stokes_paper/`** — the original January 12 2026 publication package (paper + 8 computational experiments). This is the prior pass, generated under Anthropic Claude Opus 4.5. Self-contained, with its own README, SUMMARY, and runnable Python experiments (`experiment_01..08`).

2. **`cold-derivation/`** — a sealed cold re-derivation written three and a half months later under Anthropic Claude Opus 4.7, with the prior-art conversation files quarantined until the cold pass + honest assessment + prior-art comparison were committed in writing.
   - `01-navier-stokes-cold.md` — the cold derivation
   - `02-honest-assessment.md` — Levels 1–4 epistemic assessment
   - `03-prior-art-comparison.md` — comparison after seal-break
   - `04-experiment-validation.md` — re-runs `experiment_08` from the prior-art package, identifies γ as the central open question

3. **`supplement/gamma-derivation-supplement.md`** — a first-principles derivation attempt for the γ coefficient that is hard-coded across the prior-art experiments. Closes the cold pass's open question (l) at Level 2 under stated assumptions, and reports three alternative geometric normalizations within ~20× spread.

## Two distinct directory roles in `prior-art/`

The `prior-art/` directory plays two roles, and the distinction matters:

- **Sealed conversation files at the root level** (`2026-01-13_*.md`, `Hard_Problems_01_NavierStokes.docx.txt`) — gitignored, private; the source material the cold pass was sealed against. `SEALED.md` documents the seal protocol that applied to those files during the cold pass.
- **The `rs2_navier_stokes_paper/` subdirectory** — the publishable Jan 12 publication package; tracked in git; intended for distribution under CC BY 4.0 per the package's own README.

The cold-rederivation seal applied to the conversation files. The Jan 12 publication package was generated from those conversations but is itself a finished, redistributable artifact.

## Methodology context

HP#1 is one of seven Hard Problems instances completing a cross-version cold re-derivation cycle under Anthropic Claude Opus 4.5 → 4.7. The methodology paper documenting the full cycle is in `../methodology-paper/` (currently local-only, pre-arXiv). HP#1 is the qualitative-class instance — a structural reframing argument rather than a closed-form numerical prediction.

For the other six instances, see sibling folders at the repository root: `Riemann-Hypothesis/` (HP#7), `Yang-Mills/` (HP#2), `Hierarchy-Problem/` (HP#3), `Fine-Structure-Constant/` (HP#4), `Gauge-Coupling-Constants/` (HP#5), `Master-Formula/` (HP#6).
