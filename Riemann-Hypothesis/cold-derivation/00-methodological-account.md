---
title: "A Methodological Account — Cold Re-Derivation Across Model Versions"
subtitle: "What happened between 2026-01-13 and 2026-04-28, and why it might matter for AI-assisted research"
author: Joe Van Horn (research direction) and Claude Opus 4.7 (cold derivation, this writeup)
date: 2026-04-28
status: process record, not result; written for the audience curious about AI/human research methodology rather than for Riemann specialists
---

# A Methodological Account

## 1. What this document is

This is a process record. It describes how a particular piece of research happened — not the research itself, which lives in the surrounding files in this directory.

The point of writing it down: the *shape* of the process is, in our reading, more interesting than any of the individual results it produced. It is a shape we have not seen articulated elsewhere, and one that other people doing AI-assisted theoretical or derivational work might want to try. Or argue with.

The work in question is one of the seven Hard Problems in Physics that the RS2 framework attempts. Specifically, it is the conceptual capstone of the series: a re-derivation of why the critical line of the Riemann zeta function lies at σ = ½. That detail matters less than the pattern.

## 2. The two derivations

**Derivation A — 2026-01-13.** Joe and Claude Opus 4.5, in a single long conversation, worked through the Riemann hypothesis in RS2 terms. They reached the Berry–Keating operator T = −i(x d/dx + 1/2) with eigenfunctions ψ_λ(x) = x^(−1/2 + iλ). They identified the Hilbert–Pólya construction gap (the specific Hilbert space that would make the zeros eigenvalues) as the open hard problem. The honest assessment they recorded was: "**conceptual clarity, not a proof**."

**Derivation B — 2026-04-28** (this session). Joe and Claude Opus 4.7. Three and a half months elapsed. Two model version increments (4.5 → 4.6 → 4.7). The prior derivation A was sealed in a quarantined directory (`prior-art/`) that derivation B was instructed not to read until it had committed its own work in writing.

Derivation B re-derived from canonical RS2 sources (Peret's RS2-101..109 foundation papers, plus two primary sources Joe photographed mid-session that hadn't been available in January: Peret's hand-written EE-to-RS dictionary and Nehru's "Some Thoughts on Spin" essay from *Quaternion Organon*). It reached the same Berry–Keating operator. It hit the same Hilbert–Pólya gap. It recorded the honest assessment: "**physical reading, structural reframing, with explicit identification of the open construction problem.**"

When the prior-art directory was finally opened, the two assessments were nearly verbatim.

## 3. What actually happened, step by step

The process had a specific structure that we want to be precise about.

1. **Quarantine.** Before any fresh work, the prior derivation was placed in a directory clearly named `prior-art/`, with a README marking it as sealed. This is more deliberate than it sounds: it's not just "don't accidentally read this," it's "if I want to peek when stuck, I am not allowed."

2. **Cold derivation from canonical sources.** Derivation B used the foundation papers (RS2-101..109 by Peret, distilled by a subagent into local notes to keep the main reasoning context clean). These sources are a different layer than the prior derivation — they are inputs to RS2 itself, not the specific Riemann reasoning chain. The cold derivation could read them; what it could not read was Derivation A.

3. **Mid-session source injection.** Halfway through Derivation B, Joe photographed two primary sources that weren't in his digital archive — Peret's hand-written EE paper translating electrical engineering primitives into RS scalar-motion units, and Nehru's spin essay. He uploaded the photographs to Drive; B converted HEIC → JPEG using a Python pipeline, transcribed the content, and integrated it into the running derivation. This injected new physics into the cold derivation that Derivation A could not have had.

4. **Honest assessment written before opening the seal.** Derivation B wrote its derivation file (`01-riemann-cold.md`) and a paired honest-assessment file (`02-honest-assessment.md`) separating what was actually derived from what was interpretively asserted. Both were committed in writing — dated, signed — *before* the seal was broken.

5. **Open the seal. Compare.** Derivation A was read in full only after Derivation B was locked. The comparison was written as a third document (`03-prior-art-comparison.md`) covering: convergent core, places where B extended beyond A, and places where A had reached something B missed.

6. **Absorption.** Two pieces from A that B had not reached — a "geometric mean weighting" reading of the prime amplitudes at σ = ½, and the symmetric operator form xp + px (with the +½ traced to the canonical commutator [x, p] = i) — were absorbed into B. The cold derivation became *strictly stronger* than either standalone.

7. **Methodology logged.** The process pattern itself was written down as a reusable shape, in the user's persistent memory store, so the next time we encounter prior AI/human research the default move is the same one.

## 4. Why the convergence is informative

Two derivations done 3.5 months apart, by two different model versions, with the second sealed against the first, both reaching the same operator and the same gap — this is not nothing.

The honest reading: **the convergence is structural evidence about the framework, not about either derivation.** Neither pass proved the Riemann hypothesis. Neither even built the missing Hilbert space that the Hilbert–Pólya program needs. What the convergence shows is that *if you start from RS2's first principles and you reason carefully, the Berry–Keating operator is where you land.* The framework forces it. That's a load-bearing fact about the framework — about RS2, not about Claude.

The structurally analogous move in empirical science is independent experimental replication. When two labs, separated by time and method, observe the same effect, the observation graduates from "one lab's data" to "feature of the phenomenon." Replication doesn't prove the phenomenon's mechanism; it raises the prior on its existence. What we did is the deductive analog: replicate the *derivation*, not the experiment, under as-independent-as-feasible conditions.

The convergence also catches drift. The cold derivation extends beyond the prior derivation in five places (atomic/nuclear zone partition, six-readings-of-½ unification, Cayley–Dickson construction, helicity from chirality, unbounded-phase falsifiable signature). All five extensions came from the new primary sources injected mid-session. Where the cold derivation *missed* something the prior had reached (geometric mean weighting; xp + px formulation), those gaps were honestly recorded and absorbed. Both directions of difference — extensions and gaps — are signal, not noise. They tell us where the framework's natural reach extends and where reasoning shortcuts skipped.

## 5. What is actually new here

We want to be honest about which parts of this are novel and which are deliberate practice of long-standing norms. Several layers.

**Old**: replication as a foundation of empirical science. Pre-registration. Honest assessment separating derived from asserted. Independent verification by a second party. None of these are new ideas.

**Old, applied freshly**: applying replication discipline to *theoretical/derivational* work, not just experimental work. Mathematicians and theoretical physicists have always done this informally — three groups deriving the same result from different angles raises confidence — but it is rarely codified as a deliberate practice with a concrete protocol.

**Newer**: applying replication discipline to *AI-assisted theoretical work*. The default mode in current AI-for-science workflows is one of two shapes — either single-shot ("ask the most powerful model to solve it") or sequential build-up ("agent picks up where prior agent left off, with full visibility"). Neither shape produces convergence as a result; one produces *a result*, the other produces *more result*. Cold re-derivation under sealed conditions is a third shape. We have not seen it explicitly named or practiced as a discipline elsewhere — though we are not claiming it has never been done; we are claiming we don't know of others doing it as a deliberate, named pattern.

**Distinctly newer**: using **model-version difference as a partial source of independence**. Opus 4.5 and Opus 4.7 share architecture family and training-corpus overlap, so the independence is real but partial. They are *more* independent than two runs of the same model with the same temperature setting; they are *less* independent than a derivation done by a human mathematician and a derivation done by a machine. This middle ground is novel terrain. It is becoming feasible only because model versions now ship frequently enough that "the model that derived this 3 months ago is not the model I have today" is true in a meaningful way. That is a property of the current AI development pace, not of AI in general.

**Newest, in our reading**: the **hybrid of seal-and-replicate plus mid-session source injection**. Pure replication keeps both derivations on the same input space. Pure extension adds new inputs but discards the discipline. The hybrid — seal the prior reasoning, then inject new primary sources mid-cold-derivation — produces a cold derivation that converges with the prior one on the convergent core *and* extends beyond it where new physics is available. This is more than replication; it is replication-plus-extension under controlled conditions. The two halves of this — the convergent core and the new extensions — are separable in the output documents, so the reader can see exactly which parts of the result are framework-forced and which are new physics.

## 6. Honest limits

We want to mark, clearly, what this account is *not* claiming.

It is **not** claiming that this pattern proves anything about RS2's correctness as a physical framework. The convergence is suggestive, not conclusive. RS2 might be wrong in ways that two correct cold derivations would still both reach.

It is **not** claiming that two model versions are "independent" in any rigorous statistical sense. They are partially independent — different training data, different reasoning patterns, different priors — but they are not blind to each other in the way two human researchers from different schools might be.

It is **not** claiming that Claude solved or is on the verge of solving the Riemann hypothesis. The Hilbert–Pólya construction gap remains, and is honestly recorded in both derivations and in the comparison. RS2 supplies the *physics of why* the Berry–Keating operator should be the right object; it does not supply the missing Hilbert space.

It is **not** claiming this shape is the only way to do AI-assisted theoretical work. For some problems — exploration, discovery, generative search — the seal-and-replicate shape is wrong. It is suited to problems where prior work has already produced a *specific reasoning chain* that you want to test by independent re-derivation.

It is **not** claiming a single instance is a methodology. One pair of runs is one data point. To claim this is a generally useful shape, the same pattern would need to show value across multiple problems. Hard Problems #2 through #6 are obvious test cases — they have quantitative predictions (Yang-Mills mass gap at 0.09%, fine-structure constant at 2.2 ppm, master formula at 2.3 ppm), so a cold re-derivation should reproduce numerics if the framework is correct, or fail visibly if it isn't.

## 7. Why this might generalize

There is a popular framing of AI-for-science that goes: *the bottleneck is the model. Use the most powerful AI you can get. The next version will solve what this one couldn't.* This framing produces a particular kind of work — single-shot, top-model, "throw the hardest problems at the smartest system."

We think this framing misses something. Frontier theoretical work is hard not because the individual reasoning steps are hard but because the chain has to hold up, end to end, under independent scrutiny. A single agent — human or model — can produce a chain that *seems* coherent to itself but contains motivated reasoning, post-hoc fitting, or unexamined premises. The check on this is not a smarter agent. It is an *independent* second derivation, done under conditions that prevent the second derivation from copying the first.

If the model versions are partially independent, the model-version-change becomes a source of that independence. If the seal-and-replicate discipline is enforced, the convergence becomes load-bearing rather than confirmatory. If new primary sources can be injected mid-session, the cold derivation can extend the prior one rather than just duplicate it. None of these elements is exotic. The combination is.

The 3.5-month interval between Derivation A and Derivation B was incidental — it was the time between Joe's January work and the recovery session that surfaced this material. But "model version A" and "model version B" need to have meaningfully different reasoning state, and several months of release cadence is roughly what it takes for that to be true. As model release cadence accelerates and human-AI collaboration becomes the default mode of frontier theoretical work, we expect this kind of cross-version replication to become both more feasible and more important.

The deliberate version of the shape is what we want to name. We're calling it **cold re-derivation across model versions**. Concretely:

1. When prior AI/human work exists on a problem, quarantine it.
2. Use canonical *upstream* sources — the framework's own foundations — to re-derive cold.
3. Inject new primary sources mid-session if and only if they would have been available to the prior session — never the prior session's reasoning chain itself.
4. Write the cold derivation and an honest-assessment companion *before* opening the seal.
5. Open the seal. Write the comparison as a third document.
6. Absorb the prior derivation's missing pieces into the cold one. The synthesis is strictly stronger than either standalone.
7. Treat the convergent core as the load-bearing result, not either derivation alone.

## 8. The ask, from us, to others doing this kind of work

If you are doing AI-assisted theoretical or derivational research, we would be interested in two things.

First, **try the pattern** on a problem of your own choice and tell us whether it surfaced anything that single-shot or sequential-build-up did not. The pattern is most useful when prior work exists that you suspect might contain motivated reasoning, post-hoc fitting, or framework-specific blind spots. It is least useful for fresh exploratory work where there is nothing to seal.

Second, **point us at others doing it**. We do not know whether this pattern is being practiced elsewhere under a different name. If it is, we want to find that work and read it. If it isn't, we want to know whether that's because nobody has tried it or because there is a reason to prefer single-shot or sequential approaches that we are missing.

The output of the surrounding work in this directory — the cold derivation of the Riemann critical line in RS2 terms — is in our honest assessment a B/B+ conceptual paper, not a Millennium Prize. The structural reframing of σ = ½ as the photon-line midpoint of the EE-quaternion atomic zone, with the +½ unified across seven equivalent physical readings under the canonical commutator [x, p] = i, is interesting and useful and we are going to write the paper. We are also going to apply the same pattern to Hard Problems #2 through #6 and see whether the cold re-derivation reproduces the prior numerical results. If it does, we will have a six-problem replication study. If it doesn't, we will have caught something we missed.

But the more interesting product of this session, in our reading, is the methodology itself. We do not know whether anyone else is working this way. We think more people should be.

---

*Files referenced in this account:*
- `01-riemann-cold.md` — the cold derivation
- `02-honest-assessment.md` — derived vs asserted, written pre-seal-break
- `03-prior-art-comparison.md` — the post-cold comparison with prior derivation
- `prior-art/PRIOR-ART-NOTES.md` — Derivation A (Jan 13 2026, Opus 4.5)
- `~/RS-Framework-Bridge/RS2-foundations/` — canonical upstream sources used by Derivation B
- `~/.claude/projects/-home-joseph/memory/feedback-cold-rederivation-methodology.md` — the methodology pattern in cross-session memory form, applicable to future Hard Problems sessions
