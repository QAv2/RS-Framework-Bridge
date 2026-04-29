---
seal_date: 2026-04-29
sealed_by: cold-rederivation methodology (Phase 5.4)
status: SEALED — DO NOT OPEN UNTIL COLD DERIVATION IS COMMITTED IN WRITING
---

# Sealed Prior Art — Fine Structure Constant (HP#4)

This directory contains Jan 13 / Jan 14 2026 prior-art conversations on the
fine-structure constant, generated under Opus 4.5. They are **sealed**
during the Phase 5.4 cold rederivation: the cold derivation in
`../cold-derivation/01-fsc-cold.md` is to be written *without consulting
any file in this directory*.

## Files included

- `2026-01-13_deriving-the-fine-structure-constant-from-rs-geometry_ccfcd00a.md` — main FSC derivation thread
- `2026-01-13_deriving-coupling-constants-from-rotational-geometry_6cb6664c.md` — adjacent coupling-constants derivation (HP#5 territory; sealed because likely cross-references FSC)
- `2026-01-13_empirical-formula-seeking-theoretical-grounding_71c05a1e.md` — methodological framing thread (likely the empirical-form discovery)
- `2026-01-14_bernoulli-primes-and-the-fine-structure-constant_f6140959.md` — Jan 14 follow-up: Bernoulli + primes connection
- `2026-01-13_rs-framework-bridge-coupling-derivations_67fc3925.md` — broad coupling-derivation thread

## Allowed sources during cold pass

- `~/RS-Framework-Bridge/RS2-foundations/` — the canonical RS2-101..109
  distillations + Peret EE → RS dictionary + Nehru spin paper distillations
  (these are first-principles / framework knowledge, not prior project work)
- `~/RS-Framework-Bridge/RS-research-corpus/` — Larson primary corpus +
  Nehru/Peret extended + adjacent (curated reference material, not prior
  project derivations)
- General published mathematical / physical knowledge (CODATA, PDG, etc.)

## Forbidden during cold pass

- Any file in this `prior-art/` directory
- The `~/rs2-recovery/conversations/` directory (which is what these files
  are copied from)
- Other Hard Problems prior-art directories
  (`../../Riemann-Hypothesis/prior-art/`, `../../Yang-Mills/prior-art/`,
  `../../Hierarchy-Problem/prior-art/`)

## Methodology reference

See `~/.claude/projects/-home-joseph/memory/feedback-cold-rederivation-methodology.md`
for the protocol. Three prior instances:

- Phase 5.1 Riemann (`../../Riemann-Hypothesis/cold-derivation/`) — converged on operator + gap
- Phase 5.2 Yang-Mills (`../../Yang-Mills/cold-derivation/`) — converged on `Δ = ln(2π) × p`
- Phase 5.3 Hierarchy (`../../Hierarchy-Problem/cold-derivation/`) — converged structurally + detected drift

This is instance 4. Target: highest precision in series (prior closed form
claimed 2.2 ppm against CODATA `1/α = 137.035999...`). The closed form to
test against:

> `1/α = (2π)³/2 + (2π)²/4 + (2π)/2 = 4π³ + π² + π`

is *not* itself sealed because it appears in the public memory file
(`MEMORY.md` index entry) and was visible above the seal line. The
*derivation chain* leading to it from RS2 first principles is what the
cold pass must reach independently.
