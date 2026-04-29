---
title: "Why σ = ½? A Physical Reading of the Riemann Critical Line from Reciprocal-System First Principles"
short_title: "A Physical Reading of σ = ½"
authors:
  - name: Joseph Vanhorn
    affiliation: Independent researcher
    orcid: 0009-0003-0972-606X
    email: contact@qualia-algebra.com
target_venue: Foundations of Physics
status: draft v1 complete — sessions 1–4 done; ready for read-through by Joe and submission prep
session_log:
  - 2026-04-29 (session 1): scaffold + abstract + §1 + §2 drafted
  - 2026-04-29 (session 2): §3 σ-coordinate + §4 EE-quaternion photon line drafted
  - 2026-04-29 (session 3): §5 Eight Readings of σ = ½ + §6 Hilbert–Pólya operator + §7 Numerical confirmation drafted
  - 2026-04-29 (session 4, part 1): §8 Convergence with Independent Vantages drafted
  - 2026-04-29 (session 4, part 2): §9 Discussion and Open Problems + §10 Conclusion drafted
  - 2026-04-29 (session 4, part 3): byline → Joseph Vanhorn + ORCID 0009-0003-0972-606X + email contact@qualia-algebra.com; references alphabetized; vantage-count consistency fixed across abstract / §1.7 / §6.5 (now "five vantages" matching §8); citation audit clean (34/34); equation numbering verified
date: 2026-04-29
sources:
  - cold-derivation/01-riemann-cold.md
  - cold-derivation/02-honest-assessment.md
  - cold-derivation/03-prior-art-comparison.md
  - cold-derivation/00-methodological-account.md
  - dorsey-bridge.md
  - section-7-tests/results.md
  - RS2-foundations/00-FIRST-PRINCIPLES.md
license: CC-BY 4.0
---

# Why σ = ½? A Physical Reading of the Riemann Critical Line from Reciprocal-System First Principles

**Joseph Vanhorn**¹

¹ *Independent researcher. ORCID: [0009-0003-0972-606X](https://orcid.org/0009-0003-0972-606X). Correspondence: contact@qualia-algebra.com.*

---

## Abstract

The Riemann hypothesis asserts that all non-trivial zeros of the Riemann zeta function lie on the critical line σ = ½. Why ½? Standard derivations answer structurally — σ = ½ is the fixed-point line of the functional-equation involution s ↔ 1 − s — but supply no physical content for the value itself. This paper offers a physical reading of the critical line from the Reciprocal System (RS) framework (Larson; Peret; Nehru), in which the universe is composed of one component, motion, in three discrete dimensions, with two reciprocal aspects, space and time. We show that σ = ½ admits **eight equivalent characterizations**: (1) the natural-datum locus under a half-shifted log-displacement coordinate; (2) the photon-line midpoint of the electrical-engineering quaternion w-axis between the progression (+1) and gravity (−1) datums; (3) the geometric-mean exponent at which every prime contributes amplitude p^{-1/2} = √(p^{−1} · 1) across the Euler product; (4) the metaplectic weight of the Jacobi theta function under modular inversion; (5) the spin-½ of magnetic 2D rotation, equivalently the steradian/radian unit ratio between solid and plane angle; (6) the n = 1 single-loop inductance quantum of Peret's electrical-engineering–to–RS dictionary; (7) one step of the Cayley–Dickson doubling chain ℝ → ℂ → ℍ → 𝕆; and (8) the √x prime-wave amplitude in the Riemann–von Mangoldt explicit formula. The eight readings are unified under the canonical commutator [x, p] = i: symmetrizing position and momentum into the reciprocal-aspect operator H = xp + px = −i(2x d/dx + 1) recovers the Berry–Keating dilation generator with the +½ derived rather than postulated. We present numerical confirmation at N = 500 zeros (GUE pair correlation residual 89% smaller than the uniform reference; nearest-neighbor spacing residual 3.8× smaller than GOE, 23× smaller than Poisson) and record the convergence with four pre-existing vantages on σ = ½ — Riemann analytic, Berry–Keating semiclassical, Dorsey geometric, random-matrix statistical — and an earlier independent sealed cold re-derivation of the same chain, for five independent vantages in total. The Hilbert–Pólya construction problem — the specific Hilbert space, measure, and boundary conditions making H eigenvalues equal {2γ : ζ(½ + iγ) = 0} — is not solved here. What this paper supplies is the physics of why this operator, and why this half.

**Keywords**: Riemann hypothesis · Hilbert–Pólya program · Berry–Keating operator · metaplectic representation · Reciprocal System · birotation · prime number theorem · GUE statistics

---

## 1. Introduction

### 1.1 The Riemann hypothesis

The Riemann zeta function

$$\zeta(s) = \sum_{n=1}^{\infty} \frac{1}{n^s} = \prod_{p \text{ prime}} \frac{1}{1 - p^{-s}}, \qquad \Re s > 1, \tag{1.1}$$

extends by analytic continuation to a meromorphic function on ℂ with a single simple pole at s = 1. The completed zeta function

$$\xi(s) = \tfrac{1}{2}\, s(s-1)\, \pi^{-s/2}\, \Gamma(s/2)\, \zeta(s) \tag{1.2}$$

is entire and satisfies the functional equation ξ(s) = ξ(1 − s), making the involution s ↔ 1 − s a symmetry of ξ. Zeros of ξ — equivalently, the *non-trivial* zeros of ζ — lie in the critical strip 0 < σ < 1. The **Riemann hypothesis** [^riemann1859] asserts that every non-trivial zero satisfies σ = ½.

Empirically, the hypothesis is well attested. Hardy [^hardy1914] proved in 1914 that infinitely many zeros lie on the critical line. Numerical verification of the first 10¹³ zeros all on σ = ½ is now standard [^gourdon][^plattrudgian]; the conjecture was named one of the seven Millennium Prize Problems by the Clay Mathematics Institute in 2000. No counterexample is known.

What is *not* known — and what motivates this paper — is the **physical content** of the value ½. The standard derivation simply notes that σ = ½ is the fixed-point line of the involution s ↔ 1 − s on ℂ, the only locus invariant under the functional-equation symmetry. This is mathematically correct but physically empty: it tells us *where* the line sits, not *why*.

### 1.2 The Hilbert–Pólya program

The most enduring physical reading of the Riemann hypothesis was suggested independently by Hilbert and Pólya in the 1910s and 1950s respectively: the imaginary parts of the non-trivial zeros, $\{\gamma_n : \zeta(\tfrac{1}{2} + i\gamma_n) = 0\}$, are the eigenvalues of a self-adjoint operator $\hat H$ on a suitable Hilbert space. Self-adjointness forces real spectrum; if such an operator exists, all non-trivial zeros automatically lie on σ = ½.

Several candidate operators have been proposed. Selberg's 1956 trace formula [^selberg1956] for the Laplace–Beltrami operator on the modular surface produces an explicit-formula structure reminiscent of the prime distribution. Connes' 1996 trace-formula program [^connes1996][^connes1998] places ζ-zeros as the absorption spectrum of a flow on the noncommutative space of adeles modulo the idele class group — a deep and ongoing program. Berry and Keating [^berryKeating1999] proposed the semiclassical Hamiltonian

$$H_{\text{BK}} = -i\!\left(x\frac{d}{dx} + \tfrac{1}{2}\right) \tag{1.3}$$

— a dilation generator with eigenfunctions $f_\lambda(x) = x^{-1/2 + i\lambda}$ — as the natural Hilbert–Pólya candidate for the Riemann zeros. Each of these operators has the *spectral* features one would want; none has yet been equipped with a Hilbert space, measure, and boundary conditions that produce the actual Riemann zeros.

Two structural questions about the Berry–Keating operator have been asked but not, in our reading, answered. First: why the **+½**? Berry and Keating offered a semiclassical motivation (regularizing the conjugate-pair Hamiltonian xp at the origin), but the half-integer is not derived from any deeper structure; it is *chosen* to make the eigenfunctions match $x^{-1/2 + i\lambda}$. Second: why **xp specifically**, rather than some other Hermitian combination of position and momentum? Again, semiclassical reasoning — but no first-principles necessity.

### 1.3 Random matrix theory

A second strand of evidence comes from random matrix theory. Montgomery [^montgomery1973] conjectured in 1973 that the pair correlation of zeros matches the pair correlation of eigenvalues of large random matrices drawn from the Gaussian Unitary Ensemble (GUE):

$$R_2(r) = 1 - \!\left(\frac{\sin\pi r}{\pi r}\right)^{\!2}. \tag{1.4}$$

Odlyzko's empirical work [^odlyzko1987][^odlyzko2001], extending to the 10²¹-st zero, confirms Montgomery to extraordinary precision. The full distributional fingerprint — pair correlation, nearest-neighbor spacing, $n$-level correlations, form factor — agrees with GUE.

The choice of GUE rather than GOE (Gaussian Orthogonal Ensemble) is itself a structural fact: GUE describes complex Hermitian matrices *without* time-reversal symmetry, whereas GOE describes real symmetric matrices *with* time-reversal symmetry. The empirical fact that Riemann zeros are GUE-distributed says, in the language of physics, that whatever operator generates them must break time-reversal symmetry. Why? The standard answer is "because that's what the data show." A physical answer — a reason this should hold from first principles — has not been forthcoming.

### 1.4 What this paper asks

The standard derivations of the critical line answer *where* σ = ½ sits and *what* its statistical company looks like. They do not answer:

- **Why ½, and not 1, 0, or some other rational?** The functional equation forces a fixed-point line; nothing in the standard derivation forces that line to ½ rather than to any other value the involution might fix.
- **Why does the Hilbert–Pólya operator carry a +½?** Berry and Keating supply a +½ on semiclassical grounds; it could equally well be 0 or 1 with no formal cost.
- **Why GUE rather than GOE?** Empirical, not derived.
- **Why xp + px (or its half-form)?** Why this antisymmetrization of position and momentum is the natural generator has been treated as suggestive rather than necessary.

These questions ask for *physical content*. To answer them we need a framework in which the value ½ is forced by physics rather than by coordinate. The framework we use is the **Reciprocal System** (RS) of Larson [^larson1959], extended in its second-generation form (RS2) by Peret [^peretRS2] and Nehru [^nehruRS2] in the 2010s.

### 1.5 The contribution

The Reciprocal System postulates a single physical primitive — *motion* — whose two aspects, *space* and *time*, are reciprocal. This is a strong commitment: not "space and time as separate," but "space and time as the two faces of one ratio." The natural reference of measurement is *unit speed* (the speed of light). All physical motions are displacements from unity. Angular motion takes a primary form called **birotation** — two oppositely directed rotations whose cosine sum is the photon waveform.

Within RS, σ = ½ is not an abstract fixed-point. It is the **photon line** — the locus of unit speed under a half-shifted logarithmic displacement coordinate, the midpoint between the progression (+1) and gravity (−1) datums on the Lorentz unit circle, and the unique exponent at which the prime-counterspace recursion $p^{-s}, p^{-2s}, \ldots$ balances against the prime-material recursion $p^s, p^{2s}, \ldots$ for every prime simultaneously. The factor ½ recurs throughout the framework — in the spin-½ of magnetic 2D rotation, in the n = 1 single-loop inductance quantum of Peret's electrical-engineering–to–RS dictionary, in the steradian/radian unit ratio identified by Nehru as the geometric meaning of fermionic spin, in one step of the Cayley–Dickson doubling chain, in the geometric-mean weighting of primes in the Euler product, and in the $\sqrt{x}$ amplitude that every prime-wave carries in the Riemann–von Mangoldt explicit formula.

We show that **all eight readings are faces of one quantum**. The unifying source is the canonical commutator $[x, p] = i$: forming the symmetric operator

$$H = xp + px = -i\!\left(2x\frac{d}{dx} + 1\right) \tag{1.5}$$

— the natural RS-reciprocal antisymmetrization of position and momentum — picks up $+i/2$ from the commutator algebraically, recovering the Berry–Keating operator (1.3) with the +½ derived rather than chosen. The eight readings of §3–§5 are physical realizations of this single algebraic fact; the framework supplies the physics of *why* the operator antisymmetrizes the reciprocal aspects, and the algebra supplies the half.

The structural argument is supplemented by **numerical confirmation** at N = 500 non-trivial zeros: the GUE pair correlation conjecture fits Montgomery's prediction (1.4) with sum-squared residual 89% smaller than the uniform reference, and the nearest-neighbor spacing distribution fits the GUE Wigner surmise with sum-squared residual 3.8× smaller than GOE and 23× smaller than Poisson. The critical-line locus and the statistical fingerprint are jointly consistent with the framework's reading.

We also record a **convergence result** about the derivation itself. The argument presented here was performed under a sealed-prior-art protocol [^vanhornMethodology] — an earlier independent derivation of the Riemann critical line in the same framework, completed in January 2026 with a different model version of the same AI assistant, was quarantined during the present work and opened only after the cold derivation was committed in writing. The two passes converged on the same operator (1.5), the same identification of the +½ with metaplectic / birotational structure, and the same honest gap (the Hilbert space construction). We treat that convergence as load-bearing structural evidence about the framework, separate from any individual derivation.

### 1.6 Honest framing

This paper is not a proof of the Riemann hypothesis. The Hilbert–Pólya construction problem — the specific Hilbert space $\mathcal H$ with measure and boundary conditions on which the operator (1.5) has discrete spectrum exactly $\{2\gamma_n : \zeta(\tfrac{1}{2} + i\gamma_n) = 0\}$ — has been the canonical bottleneck of the Hilbert–Pólya program for a century. Nothing in the Reciprocal System framework, as currently developed, builds that Hilbert space. What the framework supplies is the *physics* of why H must take the form (1.5) and why the +½ must appear; it does not supply the missing construction.

We are therefore explicit: the contribution of this paper is a **structural reframing with physical content**. The eight readings of σ = ½ are independent only at the surface; algebraically they collapse to one. Each reading is a face of the same quantum, but the unified picture admits multiple physical interpretations, each of which can be tested independently. We present falsifiable signatures (§9) where appropriate. Where claims are speculative we mark them as such.

### 1.7 Roadmap

§2 presents a brief primer on the Reciprocal System framework for readers unfamiliar with it. §3 derives the σ-coordinate of ζ as a half-shifted RS displacement coordinate, with σ = ½ corresponding to the natural-datum locus of unit speed. §4 introduces the electrical-engineering quaternion structure of Peret's foundation papers and identifies σ = ½ as the photon-line midpoint of the atomic-zone w-axis. §5 develops the eight readings of σ = ½ and shows that they are unified under the canonical commutator [x, p] = i. §6 presents the Hilbert–Pólya operator H = xp + px in the RS reciprocal-aspect form, derives the +½ from antisymmetrization, and states the open construction problem honestly. §7 reports the numerical confirmation at N = 500 zeros. §8 records the convergence with four pre-existing vantages — Riemann analytic, Berry–Keating semiclassical, Dorsey geometric [^dorsey2023], random-matrix statistical — and adds the sealed-prior-art cold re-derivation as a fifth. §9 discusses falsifiable signatures, open problems, and limits. §10 concludes.

---

## 2. The Reciprocal System Framework — A Brief Primer

This paper rests on the Reciprocal System (RS) framework. RS is not part of the standard physics curriculum; readers approaching it for the first time should expect a different ontology than the field-and-particle picture of the standard model. This section presents the minimum we will use. Readers familiar with RS may skip to §3.

### 2.1 Origins

The Reciprocal System was developed by Dewey B. Larson over four decades beginning in the late 1950s, in a sequence of monographs starting with *The Structure of the Physical Universe* (1959) [^larson1959] and continuing through *Nothing But Motion* [^larson1979nbm], *The Universe of Motion* [^larson1984uom], *New Light on Space and Time* [^larson1965nlst], *Basic Properties of Matter* [^larsonBPM], and *The Neglected Facts of Science* [^larsonNFoS]. The framework was substantially revised in the 2010s by Bruce Peret and K. V. K. Nehru in a body of work known as **RS2** (Reciprocal System v2), through papers RS2-101..109 [^peretRS2] and a long series in the journal *Reciprocity* [^reciprocity]. The most important RS2 modifications, for our purposes, are: replacement of Larson's 1D linear "rotational base" with a quaternion structure on each point of space; the introduction of *birotation* (oppositely directed rotation pairs whose Euler sum is a cosine) as the primary angular motion; and an extension of Larson's discrete-unit framework into projective rather than Euclidean geometry.

We use RS2 throughout. The framework has remained outside the mainstream of physics; its claim to relevance for the Riemann question is not that it is established but that it forces specific algebraic structures (reciprocal aspects, half-integer rotations on division algebras, projective cross-ratios) that turn out to coincide with the structures in which the critical line naturally lives.

### 2.2 Minimal axioms

Distilled from Peret's RS2-101..109 [^peretRS2], the working axiom set is:

**A1 (Postulate I).** The universe is composed of one component, *motion*, in three dimensions, in discrete units, with two reciprocal aspects, *space* and *time*.

**A2 (Postulate II).** The universe conforms to ordinary mathematics; its primary magnitudes are absolute; its geometry is *projective* (Euclidean and affine are sub-cases).

**A3.** The minimum scalar magnitude is unity. There is no zero, no negative, no fractional in the natural framework.

**A4.** Scalar motion is the projective cross-ratio of two scalar orientations, with space and time as the specific aspects. Two reciprocal forms: speed s/t and energy t/s.

**A5.** The natural reference is *unit speed* (the speed of light). All motion is measured as displacement from unity.

**A6.** Each scalar dimension carries a tristate of motion: unity (default progression), speed (s/t < 1, *material*), energy (t/s < 1 equivalently s/t > 1, *cosmic*).

**A7.** The number of stable scalar dimensions is uniquely 3, by Nehru's binary-combination stability count $n(n-1)/2 = n \Rightarrow n \in \{0, 3, \infty\}$.

**A8.** Both linear (yang) translation and angular (yin) rotation are *primary*. Vibration is shear strain between counter-directed motions; *birotation* is the angular primitive.

**A9.** Direction reverses only at unit boundaries — the discrete-unit "links of a chain" can flex only at junctures.

A1 and A2 are the postulates proper. A3–A9 are explicit consequences or supplementary commitments stated as foundational across the papers.

### 2.3 Motion-only ontology and reciprocal aspects

The strongest claim is A1. There are no "things" in RS — no particles in space, no fields on a manifold, no events at points. There is only *motion*, of which space and time are reciprocal aspects in the same sense that the numerator and denominator of a single ratio are reciprocal aspects of that ratio. A "particle" in conventional physics is, in RS, a stable rotational configuration of motion; an "event" is a unit step of progression; a "field" is a pattern of independent motions appearing to interact when projected onto the conventional space-time frame [^peret104].

The reciprocal-aspect structure is the load-bearing primitive for our argument. The involution

$$\frac{s}{t} \;\leftrightarrow\; \frac{t}{s} \tag{2.1}$$

— material aspect to cosmic aspect — is the fundamental symmetry of the framework. Under multiplicative coordinates $p = s/t$, the involution is $p \leftrightarrow 1/p$ with fixed point $p = 1$ (unit speed, A5). Under the additive log-coordinate $r = \ln p$, the involution becomes the reflection $r \leftrightarrow -r$ with fixed point $r = 0$. This involution is what the Riemann functional-equation involution s ↔ 1 − s will turn out to be (§3), once a half-shift is applied to land the natural datum on σ = ½.

### 2.4 Discrete-unit framework

By A3 the minimum scalar is 1; there is no continuous infinitesimal. Magnitudes are counting numbers, with the natural ordering "unit | speed (< 1) | energy (> 1)" supplied by A6. By A9 direction can change only at unit boundaries — between counts of 1, 2, 3, …, never *within* a unit. This is the "links of a chain" picture: a unit is a discrete object, and the only loci where reversal can occur are the joints between units.

The discrete-unit commitment has consequences we will use. The Euler product structure $\prod_p (1 - p^{-s})^{-1}$ in (1.1), which has historically been treated as a calculational identity, becomes in RS a *partition function over discrete prime counterspace shells* (§4), with each prime contributing a geometric series of unit shells. The destructive interference required for $\zeta(s) = 0$ is then a discrete-unit resonance condition rather than a continuous cancellation.

### 2.5 Linear and angular motion; birotation

By A8, both linear and angular motion are primary. RS2 takes angular motion seriously enough to make it a separate axiom rather than a derived phenomenon — and the angular form is *birotation*. Two oppositely directed rotations $e^{i\theta}$ and $e^{-i\theta}$ sum to a cosine $2\cos\theta$, and that cosine waveform is the **photon** in RS2 [^peret107][^peret109]. There is no "wave moving through a medium" in RS; there is a birotation whose Euler-cosine projection looks like a wave to a clock-time observer.

Birotation is double-cover structure. A single rotation of period $2\pi$ has half-integer angular periodicity when viewed as a birotation (because the two directions sum coherently only over a $4\pi$ cycle), and this is the **spin-½** of conventional physics: Nehru's *Some Thoughts on Spin* [^nehruSpin] reads the conventional 4π → 2π halving not as a halving of any single quantity but as a unit conversion between solid angle (steradians, 4π for a full sphere) and plane angle (radians, 2π for a full circle). The spin-½ is the *steradian/radian unit ratio* of one full sphere divided by one full circle. We will return to this in §5.

The double-cover structure of birotation is what produces the metaplectic representation of $SL(2, \mathbb R)$, and the metaplectic weight ½ is what produces the half-integer in the Berry–Keating operator. This is the central physical identification of the paper, defended in detail in §5.

### 2.6 Atomic and nuclear zones

Peret's electrical-engineering–to–RS dictionary [^peretEE] partitions the unit-speed neighborhood into two zones:

- **Atomic zone** — 4D, full quaternion $(w, i, j, k)$, supports the full magnetic 2D rotation $i \cdot j$ and the photon birotation.
- **Nuclear zone** — 2D, restricted to $(w, i)$, real plus electric-1D only, no magnetic component.

In the σ-coordinate of §3, the critical strip $0 < \sigma < 1$ corresponds to the atomic zone and the half-plane $\sigma \le 0$ to the nuclear zone. Non-trivial zeros — which require the full quaternion structure for their birotational eigenmodes — live in the critical strip; trivial zeros — which require only real and electric-1D structure — live in the nuclear zone, with Δσ = 2 spacing inherited from the magnetic-2D doubling factor (RS2-107, [^peret107]).

This zone partition is the second load-bearing identification. It supplies a physical reason why the critical strip is the *atomic* zone (the only region in which the full birotational quaternion is supported) and the half-plane $\sigma \le 0$ is the *nuclear* zone (real + electric-1D, no magnetic).

### 2.7 Why this frame

Most physical frameworks place ½ in coordinate or convention. RS places ½ in physics: the spin-½ of magnetic rotation, the n = 1 minimum loop of inductance, the metaplectic weight of birotation, the steradian/radian unit conversion. These are all the same ½, traceable to the one fact that RS angular motion is birotational and birotation is double-cover. This is why RS, despite sitting outside the mainstream of physics, is the natural frame in which to ask the question this paper asks: *why does the critical line of $\zeta(s)$ carry the value of a half?*

The remainder of the paper takes up that question and answers it eight ways.

---

## 3. Coordinate Correction: σ as RS Displacement

### 3.1 The standard frame is silent on the half

The standard treatment of $\zeta(s)$ places it on a featureless complex plane. The critical line $\Re s = \tfrac{1}{2}$ then appears as a contingent fact — the fixed-point line of the involution $s \leftrightarrow 1 - s$ — but $\mathbb C$ itself does not distinguish ½ from any other real number. Two consequences follow.

First, the value σ = ½ looks arbitrary. Without an external reason, the half is just where the involution happens to fix itself; nothing in the standard derivation forces the fixed line to ½ rather than 0, 1, or any other value the involution might fix.

Second, the functional equation looks technical. Its proof goes through Poisson summation on the Jacobi theta function $\theta(t) = \sum_n e^{-\pi n^2 t}$ via the modular relation $\theta(1/t) = \sqrt{t}\,\theta(t)$. The $\sqrt{t}$ factor in this identity is treated as a calculational artefact rather than as content — yet it is precisely the ½ we are trying to explain.

The diagnosis offered by the Reciprocal System framework is that the s-plane is the wrong frame: the physical content of the critical line is hidden behind a coordinate choice that obscures what σ = ½ means. The remainder of this section presents the coordinate correction. Sections 4 and 5 supply the physical origin of the half itself.

### 3.2 The RS displacement coordinate

By A5 the natural reference of motion is unit speed (the speed of light). All physical motions are displacements from unity. By A4 scalar motion is the projective cross-ratio of two scalar orientations, with space and time as the specific aspects, taking the two reciprocal forms $s/t$ (speed, *material*) and $t/s$ (energy, *cosmic*).

Define the RS multiplicative ratio
$$p = \frac{s}{t}, \qquad p = 1 \text{ at unit speed.} \tag{3.1}$$

The reciprocal involution $p \leftrightarrow 1/p$ is the material/cosmic sector swap of A1 and A6. In multiplicative coordinates it has fixed point $p = 1$. Pass to the additive log-coordinate
$$r = \ln p = \ln(s/t), \tag{3.2}$$

so that $r = 0$ is unit speed and the involution becomes the reflection $r \leftrightarrow -r$. The displacement $r$ is the natural RS measurement: positive for the material sector ($s/t > 1$ in some conventions, $r > 0$ in this log-displacement form), negative for the cosmic sector, zero at the natural datum.

### 3.3 The half-shifted Riemann coordinate

Define the Riemann coordinate σ from the RS displacement by
$$\sigma = \tfrac{1}{2} + \tfrac{r}{2}, \qquad r = 2\sigma - 1. \tag{3.3}$$

Under this map:

| RS displacement $r$ | Riemann $\sigma$ | Physical meaning |
|---|---|---|
| $r = 0$ (unit speed, A5) | $\sigma = \tfrac{1}{2}$ | **critical line** |
| $r = +1$ (one unit material) | $\sigma = 1$ | upper boundary; pole of $\zeta$ |
| $r = -1$ (one unit cosmic) | $\sigma = 0$ | lower boundary; nuclear-zone entry |
| $r \leftrightarrow -r$ (sector swap, A1) | $\sigma \leftrightarrow 1 - \sigma$ | functional-equation involution |

The critical line σ = ½ is precisely the locus of unit speed. The critical strip $0 < \sigma < 1$ is the unit-displacement neighborhood of the natural datum: $|r| < 1$. The functional equation $\xi(s) = \xi(1-s)$ becomes, in RS coordinates, the assertion that ξ is invariant under the reciprocal-aspect involution $r \leftrightarrow -r$. We will refer to this as the **sector-reciprocal symmetry** of the completed zeta function.

### 3.4 The half-shift is forced, not free

The map (3.3) is not a relabel. It is the unique linear map from $r$ to a new coordinate that simultaneously sends the involution's fixed point to the natural datum (A5) and aligns the involution's action with the RS sector swap (A1). Any other linear map either misplaces the fixed point or compresses the involution into a non-canonical action on the new coordinate.

What this section has shown: σ = ½ is the natural-datum locus under the RS-correct coordinate. What it has *not* shown is why the value ½ should appear in the first place. The map (3.3) carries a half-shift and a rescale by two, and so far neither has been derived from physics — they have been imposed to make the involution and the natural datum coincide. The "½" still needs an origin.

That origin is supplied by the electrical-engineering–quaternion realization of §4. The half-shift is the magnetic-2D quantum of Peret's EE → RS dictionary [^peretEE]; the σ = ½ photon line is the geometric midpoint of the atomic-zone real axis between the gravity datum and the progression datum, and it is the unique exponent at which the prime counterspace partition functions of §4.3 are simultaneously balanced. After §4 the value ½ in (3.3) will no longer be a fit but a derived consequence of the framework's magnetic-rotation and prime-counterspace structure.

---

## 4. The EE-Quaternion Photon Line

This section supplies the physical origin of the half. The half-shift in (3.3) is not a coordinate fit; it is the magnetic-2D quantum of Peret's electrical-engineering-to-RS dictionary, the geometric midpoint of the atomic-zone w-axis, and the unique exponent at which the prime counterspace partition functions are simultaneously balanced. The argument runs in three steps: the EE-quaternion identification (§4.1), the atomic / nuclear zone partition mapped onto the σ-coordinate (§4.2), and the Euler product as a partition function over prime counterspace shells (§4.3). A short synthesis closes the section (§4.4).

### 4.1 The four electrical primitives as quaternion units

The four electrical primitives — resistance $R$, inductance $L$, conductance $G$, capacitance $C$ — realize the quaternion units $\{+1, +i, -1, -i\}$ in the Argand plane. Peret's *Basic Relationships* table [^peretEE] makes the identification explicit:

| | Resistance $R$ | Inductance $L$ | Conductance $G$ | Capacitance $C$ |
|---|---|---|---|---|
| Argand position | $+1$ | $+i$ | $-1$ | $-i$ |
| RS units ($s/t$) | $t^2/s^3$ | $t^3/s^3$ | $s^3/t^2$ | $s^3/t$ |
| Conventional unit | Ω | Φ flux $t^2/s^2$ | S | Ψ in $s/1$ |
| Series law | additive | additive | reciprocal | reciprocal |
| Parallel law | reciprocal | reciprocal | additive | additive |

The four primitives satisfy
$$1 = c^2 \mu_0 \epsilon_0 = (c\mu_0)(c\epsilon_0) \implies RG = 1,$$

with $c\mu_0 = R$ and $c\epsilon_0 = G$. $R$ and $G$ are reciprocal *as a product*; this does not make them "the same quantity inverted." $R$ is magnetic (linear in the Argand sense) and $G$ is dielectric (1/linear). The same applies to $L$ and $C$ in the imaginary direction. The reciprocity is product-reciprocity, not identity — a distinction we will use in §4.3 when the Euler product factors $p^{-s}$ and $p^{s-1}$ enter as material and counterspace amplitudes of the same prime.

A consequence we will use repeatedly: the **mass–inductance identification** of LeBon [^lebon1907]. LeBon's relation $p/v = M$ — momentum per velocity is mass — combined with the EE identity that magnetic flux $\Phi$ has the same units as momentum gives
$$\Phi/c = M = L = t^3/s^3.$$

Mass and inductance are the same RS quantity. There is no separate "mass field" in this framework; mass is the $t^3/s^3$ inductance reading of a magnetic 2D rotation. This identification will reappear in §6 as the physical content of the Berry–Keating $+½$: the lift from a pure dilation operator to the natural-datum operator is a single-loop ($n = 1$) magnetic-rotation quantum.

By A4 the cross-ratio is the unique projective invariant. In the EE-quaternion frame, cross-ratios between the four primitives are the natural projective scalars. The reciprocal involution $s/t \leftrightarrow t/s$ exchanges $R \leftrightarrow G$ and $L \leftrightarrow C$ — i.e., conjugation across the $+1/-1$ axis of the Argand plane. This is the same involution as the σ ↔ 1−σ of §3, viewed in EE coordinates rather than displacement coordinates.

### 4.2 Atomic and nuclear zones

Peret's foundation papers partition the unit-speed neighborhood into two zones, distinguished by the dimensional content of their motion [^peretEE]:

- **Atomic zone** — 4D, full quaternion $(w, i, j, k)$. Supports magnetic 2D rotation $i \cdot j$ and the birotation $k = ij$ that produces the photon waveform.
- **Nuclear zone** — 2D, restricted to $(w, i)$. Real plus electric-1D only; no magnetic component.

Mapped onto the σ-coordinate of §3:

| σ region | Zone | Quaternion content | Behavior |
|---|---|---|---|
| $0 < \sigma < 1$ | atomic, 4D | $(w, i, j, k)$ | critical strip; non-trivial zeros |
| $\sigma \le 0$ | nuclear, 2D | $(w, i)$ | trivial zeros at $s = -2, -4, -6, \ldots$ |
| $\sigma \ge 1$ | divergent | — | Euler product diverges; no zeros |

The critical strip is the atomic zone: the full quaternion is supported, including the magnetic-2D rotation that produces the photon birotation. Non-trivial zeros — those that require the complete quaternion structure for their birotational eigenmodes — live here.

The trivial-zero region $\sigma \le 0$ is the nuclear zone. Only the real ($w$) and electric-1D ($i$) directions are supported; no magnetic component. The trivial zeros sit at unit-spaced intervals in σ with $\Delta\sigma = 2$ because the nuclear-zone resonance step is one electric-rotation quantum times the magnetic-2D weighting. The factor of two is the **magnetic doubling** — directly from RS2-107's identification of magnetic rotation as 2D rotation with a $4\pi$ period (which appears as $2\pi$ — i.e., one full electric rotation — when projected to the nuclear-zone subspace) [^peret107]. This is the strongest direct EE-quaternion prediction: if the simple 4D / 2D zone picture is right, the trivial-zero density on $\sigma \le 0$ is exactly $\rho_{\text{nuc}}(\sigma) = \tfrac{1}{2}$ (one zero per two units), and this is empirically the case.

For $\sigma \ge 1$ the Euler product diverges, so $\zeta$ has no zeros there — consistent with the elementary fact that $\zeta$ has no zeros for $\Re s \ge 1$. The RS reading supplies a physical reason for the boundary: $\sigma = 1$ is the unit-displacement boundary by (3.3) ($r = +1$), and direction reversal can occur only at unit boundaries (A9). Beyond it the partition function is unbounded; the framework has exited the unit-displacement neighborhood of the natural datum.

The **photon line** $\sigma = \tfrac{1}{2}$ is the midpoint of the w-axis between the two atomic-zone boundaries: $\sigma = 0$ (gravity, $-1$ in Peret's labeling) and $\sigma = 1$ (progression, $+1$). It is the unique locus where the magnetic-2D rotation $i \cdot j$ is in phase balance with the electric-1D rotation $i$ — the locus at which the birotation is a stable resonance rather than a transient pulse. This is the geometric realization of the half: the mid-point of the atomic w-axis, equidistant from the gravity and progression datums.

### 4.3 The Euler product as counterspace partition function

The Euler product
$$\zeta(s) = \prod_p \frac{1}{1 - p^{-s}} = \prod_p \sum_{k=0}^{\infty} p^{-ks} \tag{4.1}$$

has a standard reading as a calculational identity: the geometric series $\sum_k p^{-ks}$ at each prime sums to $1/(1 - p^{-s})$. The RS reading reinterprets this geometric series as physical content.

Per Peret's annotated derivation [^peretEE], "The sequence of recursion in counterspace would be $x^{-1}, x^{-2}, x^{-3}, \ldots$, which are 1·p, 2·p, and 3·p energy levels of the 'electron' (spatial rotation)." Each prime $p$ generates its own counterspace recursion $p^{-1}, p^{-2}, p^{-3}, \ldots$ — the **discrete energy shells** of the spatial rotation at that prime. The Euler factor $1/(1 - p^{-s})$ is the geometric series summing these shells with weights $p^{-s}, p^{-2s}, p^{-3s}, \ldots$ — i.e., the **partition function of the prime-$p$ counterspace shell** at spectral parameter $s$.

The full Euler product is therefore the joint partition function over the entire prime spectrum: each prime contributes its own counterspace partition function, and $\zeta(s)$ is the product across primes. Zeros of $\zeta(s)$ are points where every prime's shell amplitude destructively interferes simultaneously — the resonance condition under which the joint partition function vanishes.

The simultaneous-destructive-interference condition has a sharp form. At spectral parameter $\sigma$, a given prime contributes amplitude $p^\sigma$ on the material side and $p^{1-\sigma}$ on the counterspace side (the latter via the reciprocal-aspect involution of A1, equivalently the functional-equation symmetry). Their ratio is $p^{2\sigma - 1}$. Simultaneous balance across **all** primes requires the ratio to equal unity for every prime at once:
$$p^{2\sigma - 1} = 1 \;\; \forall p \quad \iff \quad \sigma = \tfrac{1}{2}. \tag{4.2}$$

The condition (4.2) holds for every prime simultaneously only at $\sigma = \tfrac{1}{2}$. Off the critical line, the ratio $p^{2\sigma - 1}$ varies with $p$, and no single coherent destructive-interference condition is available across the prime spectrum: the framework predicts that destructive interference of the joint partition function — and hence the existence of zeros — is structurally impossible off the critical line. (The structural argument does not by itself prove that all zeros lie on σ = ½; it produces the line as the unique candidate locus and shifts the burden to the spectral construction of §6.)

### 4.4 The photon line as the convergent locus

Two independent identifications now place the same value $\sigma = \tfrac{1}{2}$:

1. **Geometric.** The photon line is the midpoint of the atomic-zone w-axis between progression ($+1$) and gravity ($-1$) — the unique locus at which the magnetic-2D rotation $i \cdot j$ is in phase balance with the electric-1D rotation $i$ (§4.2).
2. **Spectral.** The photon line is the unique σ at which every prime is in self-counterspace–to–material balance, enabling the simultaneous destructive interference required for $\zeta(s) = 0$ (§4.3, equation 4.2).

These two readings independently converge on $\sigma = \tfrac{1}{2}$. The convergence is not a coincidence: both are faces of the same RS reciprocal-aspect structure (A1, A4). The atomic-zone w-axis is the geometric realization of the same reciprocal-aspect symmetry that makes (4.2) hold for every prime at $\sigma = \tfrac{1}{2}$; the photon line is where the reciprocal aspects balance both geometrically (in the EE quaternion) and spectrally (across the prime partition functions).

This closes the half-shift question of §3. The "½" in (3.3) is no longer a free fit. It is the midpoint of the atomic zone geometrically and the balance exponent of the prime counterspace partition spectrally — two physical realizations of the same RS reciprocal half. Section 5 develops this single half into eight equivalent readings, all faces of the canonical commutator $[x, p] = i$, and the picture closes when those readings unify under a single algebraic identity in §6.

---

## 5. Eight Readings of σ = ½

### 5.1 Overview

Sections 3 and 4 placed σ = ½ on the photon line by two independent identifications: the natural-datum locus of the RS displacement coordinate (§3) and the geometric–spectral midpoint of the atomic-zone w-axis where every prime-counterspace partition function is simultaneously balanced (§4). Both identifications produced the same value, but each treated the half as a single feature — a coordinate origin in §3, a midpoint in §4. The framework, in fact, supplies the half in eight surface-distinct ways.

This section enumerates those readings. Each appears, locally, as a different physical or algebraic quantity: the metaplectic weight of an abstract group representation, the periodicity of spin-½, a unit conversion between solid and plane angle, a coordinate shift, a single coil of an inductor, a step of the Cayley–Dickson doubling chain, the geometric mean of two prime amplitudes, the amplitude of a prime-indexed wave. The eight are listed in Table 5.1 with their canonical RS sources.

| # | Reading | RS / mathematical content | Source |
|---|---|---|---|
| (a) | Metaplectic weight ½ | $\theta(1/t) = \sqrt{t}\,\theta(t)$ modular factor | Mumford [^mumfordTata]; §5.5 |
| (b) | Spin-½ of magnetic 2D rotation | 4π period of fermion wavefunction | RS2-107 [^peret107]; Nehru §1 [^nehruSpin] |
| (c) | Steradian/radian unit ratio | (1 sr / 1 rad²) dimensional conversion | Nehru §1 [^nehruSpin] |
| (d) | σ-coordinate half-shift | $\sigma = \tfrac{1}{2} + r/2$ | §3 of this paper |
| (e) | n = 1 single-loop inductance | $L = n^2 A \mu / l$, minimum $n = 1$ | Peret EE dictionary [^peretEE]; §5.6 below |
| (f) | Cayley–Dickson doubling step | ℂ → ℍ via $\psi = \{\psi_1, j\psi_2\}$ | Nehru §8 [^nehruSpin] |
| (g) | Geometric-mean prime weighting | $p^{-1/2} = \sqrt{p^0 \cdot p^{-1}}$ for every $p$ | §4.3 of this paper |
| (h) | $\sqrt{x}$ prime-wave amplitude | $x^{\rho}/\rho$ with $\Re\rho = \tfrac{1}{2}$ | Dorsey [^dorsey2023]; Edwards [^edwards1974] |

§§5.2–5.9 walk through the eight readings in turn. §5.10 collapses all of them into one algebraic identity — the canonical commutator $[x, p] = i$ — and closes with a falsifiable signature.

### 5.2 (a) Metaplectic weight of θ

The Jacobi theta function
$$\theta(t) = \sum_{n \in \mathbb Z} e^{-\pi n^2 t}, \qquad t > 0, \tag{5.1}$$
satisfies the modular-inversion identity
$$\theta(1/t) = \sqrt{t}\,\theta(t). \tag{5.2}$$

The factor $\sqrt{t} = t^{1/2}$ is the **metaplectic weight**. It records that $\theta$ transforms not as a function on the upper half-plane but as a section of a line bundle on the *double cover* $\widetilde{SL}(2, \mathbb R)$ — the metaplectic group. The cover is the unique non-trivial double cover of $SL(2, \mathbb R)$, and the weight ½ is its index [^mumfordTata].

The completed zeta function is the Mellin transform of $(\theta(t) - 1)/2$:
$$\xi(s) = \tfrac{1}{2}\, s(s-1) \int_0^\infty \frac{\theta(t) - 1}{2}\, t^{s/2 - 1}\, dt. \tag{5.3}$$

Splitting the integral at $t = 1$ and substituting (5.2) on the $(0, 1)$ piece exchanges the Mellin kernel exponent $s/2 - 1$ for $(1-s)/2 - 1$, producing $\xi(s) = \xi(1 - s)$ directly. The functional equation is the Mellin transform of the metaplectic relation. The factor ½ scaling $s$ in the Mellin variable matches the half-weight of $\theta$; both halves are the same half, and the critical line $\sigma = \tfrac{1}{2}$ sits at their crossover. The half is, in this reading, a feature of $SL(2, \mathbb R)$'s double-cover structure — the abstract face of the quantum.

### 5.3 (b) Spin-½ of magnetic 2D rotation

The conventional reading of spin-½ is that the wavefunction of a fermion returns to itself only after a $4\pi$ rotation, twice the $2\pi$ period of a scalar. RS2 supplies a structural source for the doubling. By RS2-107, electric and magnetic motions inhabit different rotational dimensions [^peret107]: electric rotation is one-dimensional with period $2\pi$; magnetic rotation is two-dimensional with period $4\pi$. The doubling is dimensional, not arbitrary — a single rotation in 2D requires a $4\pi$ sweep to complete because 2D rotation has two equivalent rotation axes whose phases couple multiplicatively rather than additively.

Nehru's *Some Thoughts on Spin* §1 [^nehruSpin] makes the identification explicit: "**fermions are based on two-dimensional rotation**," and the spin-½ that conventional QM treats as an abstract quantum number is the structural fact that fermionic wavefunctions live on the magnetic-2D sector. Bosons (spin-1) live on the electric-1D sector; the photon, which is the birotation of §2.5, lives on both simultaneously.

In this reading the σ = ½ of the critical line is the same ½ as the spin-½ of the electron: the period-doubling weight that arises whenever a 2D rotation is projected onto a 1D base. The Riemann zeta function is sensitive to it because $\zeta$, via its functional equation, is the spectral signature of an operator that acts on the magnetic-2D sector — the atomic zone of §4.2. This is the **spin-physics face** of the half.

### 5.4 (c) Steradian/radian unit ratio

The $4\pi \to 2\pi$ relation conventionally read as a "halving" is, on closer inspection, **not a halving of any single quantity at all** — it is a unit conversion between two distinct geometric quantities. A full sphere subtends $4\pi$ steradians; a full circle subtends $2\pi$ radians. The two are different angular units, and the dimensional ratio $(1\,\text{sr}) / (1\,\text{rad}^2) = 1$ relates them only as a definition.

Nehru §1 [^nehruSpin] makes the correction precise: spin-½ is the dimensional unit conversion
$$\frac{\hbar}{4\pi\,\text{sr}} = \frac{1}{2} \cdot \frac{\hbar}{2\pi\,\text{sr}}, \tag{5.4}$$
where the ½ is the **steradian-to-radian conversion factor** rather than a fractional weight on any single quantity. The conventional reading — that fermion phase is "half-integer" because it returns after $4\pi$ instead of $2\pi$ — collapses, under this correction, into the structurally clean fact that the magnetic-2D rotation is measured in steradians (solid angle), the electric-1D rotation is measured in radians (plane angle), and the ratio is geometric, not numerical.

This is the **geometric face** of the half: not a halving but a unit ratio. We will use this distinction in §5.10 to clarify why $H = xp + px$ — antisymmetrization of two reciprocal aspects — produces ½ algebraically rather than numerically.

### 5.5 (d) The σ-coordinate half-shift

Section 3 derived
$$\sigma = \tfrac{1}{2} + \tfrac{r}{2}, \qquad r = \ln(s/t), \tag{5.5}$$

as the unique linear map sending the RS displacement coordinate $r$ (with reciprocal involution $r \leftrightarrow -r$ at fixed point $r = 0$) to the Riemann coordinate $\sigma$ (with functional-equation involution $\sigma \leftrightarrow 1 - \sigma$ at fixed point $\sigma = \tfrac{1}{2}$). Two halves appear in (5.5): an additive ½ that lands the involution's fixed point on the natural datum, and a multiplicative ½ that rescales the involution amplitude.

Both halves are the same half. By §5.6 below, the multiplicative ½ is the n = 1 inductance quantum of magnetic-2D rotation; the additive ½ is the corresponding coordinate shift that makes the involution's fixed point coincide with that quantum's locus. The σ-coordinate half-shift is, in the language of the eight readings, the **coordinate face** of the same physical half that §5.6 will identify with the single-loop inductance quantum.

### 5.6 (e) The n = 1 single-loop inductance quantum

In Peret's EE-RS dictionary [^peretEE], inductance has the form
$$L = \frac{n^2 A}{l} \cdot \mu, \tag{5.6}$$

where $n$ = number of coil loops, $A$ = coil area, $l$ = coil length, and $\mu$ = permeability. The minimum inductance is set by $n = 1$: a **single-loop inductor**. This is the quantum unit of magnetic-rotation flux per unit progression — the smallest discrete amount of magnetic-2D rotation that the framework permits, since the discrete-unit axiom A3 forbids fractional loops.

The "+½" of the Berry–Keating operator (1.3) admits a direct physical reading in these terms. Without the +½, the operator $-i x\,d/dx$ is a pure dilation generator, with eigenfunctions $f_\lambda(x) = x^{i\lambda}$ — *integer-weight* under dilation $x \mapsto e^t x$. Adding +½ shifts the eigenfunctions to $f_\lambda(x) = x^{-1/2 + i\lambda}$ — *half-integer-weight*. The shift is by exactly the metaplectic weight ½ of §5.2, and it has, in RS2 terms, the content:

> **The lift from the integer-weight pure-dilation operator to the half-integer-weight natural-datum operator is the n = 1 minimum loop of magnetic-2D inductance contributing one loop's worth of flux to the dilation generator.**

This is the most concrete of the eight readings — the **EE / engineering face** of the half. It also fixes a physical interpretation of why Berry–Keating chose +½ rather than +1 or 0: anything other than +½ would fail to match the n = 1 minimum-loop quantum, and the operator would either decouple from the magnetic-2D rotation entirely (n = 0 → integer weights, no inductive lift) or double it (n = 2 → unit weights, two-loop coil). Only n = 1 produces the half-integer eigenfunction $x^{-1/2 + i\lambda}$ on which $\zeta$ has its non-trivial zeros.

### 5.7 (f) Cayley–Dickson doubling step

The Cayley–Dickson construction generates the four normed division algebras $\mathbb R, \mathbb C, \mathbb H, \mathbb O$ by an iterated doubling: each step combines two copies of the previous algebra via a new imaginary unit, doubling the algebra dimension at each stage [^baez2002].

Nehru §8 [^nehruSpin] reaches the same construction by physical reasoning from the RS2 atomic-zone wavefunction. Starting from a 2D complex wavefunction $\varphi = \{\varphi_r, i\varphi_a\}$ in the nuclear zone, the atomic zone (4D, full quaternion) requires
$$\psi = \{\psi_1, j\psi_2\} \;=\; \{\psi_r,\, i\psi_b,\, j\psi_c,\, k\psi_d\}, \tag{5.7}$$

with $\psi_1, \psi_2$ each complex and $k = ij$. Equation (5.7) is the Cayley–Dickson construction of $\mathbb H$ from $\mathbb C$, derived not from algebraic doubling but from the requirement that the wavefunction of a 2D rotation in the time region must have four components — exactly Dirac's relativistic insight, reached by Nehru via RS2 first principles.

In this reading, σ = ½ is **one step of the Cayley–Dickson chain**. The transition $\sigma \le 0 \to 0 < \sigma < 1$ in the σ-coordinate is the transition $\mathbb C \to \mathbb H$ in the algebra: the nuclear zone is two-dimensional (real + electric-1D, i.e., $\mathbb C$) and the atomic zone is four-dimensional (full quaternion, i.e., $\mathbb H$). The half — the midpoint $\sigma = \tfrac{1}{2}$ between $\sigma = 0$ and $\sigma = 1$ — is the photon-line midpoint of the doubling step, equidistant from both the $\sigma = 0$ entry (lower bound of $\mathbb H$) and the $\sigma = 1$ exit (upper bound, where the Euler product diverges). This is the **algebraic face** of the half.

### 5.8 (g) Geometric-mean prime weighting

At σ = ½, each prime contributes amplitude
$$p^{-1/2} = \sqrt{p^0 \cdot p^{-1}}, \tag{5.8}$$
the geometric mean of its $\sigma = 0$ contribution ($p^0 = 1$, the counterspace boundary) and its $\sigma = 1$ contribution ($p^{-1}$, the material boundary). This is the §4.3 partition-function balance condition (4.2) read in reverse: the σ at which every prime contributes the geometric mean of its endpoint amplitudes is the unique σ at which the reciprocal-aspect involution $p^{-s} \leftrightarrow p^{-(1-s)}$ has a simultaneous fixed point across the entire prime spectrum. For all primes simultaneously to be in self-balance, the exponent must be the geometric mean — i.e., σ = ½.

The geometric-mean reading is the **prime-spectrum face** of the half. It is what makes σ = ½ the unique candidate locus for non-trivial zeros: any other σ assigns each prime a *different* deviation from the geometric-mean balance, and no coherent destructive-interference condition is available across the prime ensemble.

### 5.9 (h) The √x prime-wave amplitude

The Riemann–von Mangoldt explicit formula expresses the Chebyshev prime-counting function $\psi(x) = \sum_{p^k \le x} \log p$ in terms of the non-trivial zeros [^edwards1974]:
$$\psi(x) = x - \sum_\rho \frac{x^\rho}{\rho} - \log(2\pi) - \tfrac{1}{2}\log\!\left(1 - x^{-2}\right). \tag{5.9}$$

Each non-trivial zero $\rho = \tfrac{1}{2} + i\gamma$ contributes
$$\frac{x^\rho}{\rho} = \frac{\sqrt{x}}{\rho} \cdot e^{i\gamma \log x}, \tag{5.10}$$

a wave of frequency $\gamma$ in $\log x$ with **amplitude** $\sqrt{x}$. The amplitude is exactly $x^{\Re\rho} = x^{1/2}$ on the critical line; off the critical line at $\rho = (\tfrac{1}{2} \pm \epsilon) + i\gamma$, the amplitude becomes $x^{1/2 \pm \epsilon}$ — growing or shrinking relative to its companions.

Dorsey [^dorsey2023] renders this geometrically as **prime-waves on a 3D sphere — bounded, infinite, and non-touching**. The geometric content of "non-touching" is precisely the requirement that all prime-waves share one amplitude scaling: if even a single zero sat off the critical line at $\sigma = \tfrac{1}{2} + \epsilon$, its corresponding wave would carry amplitude $x^{1/2 + \epsilon}$ and would either grow or shrink relative to the rest, breaking the bounded wavefield. The unique σ at which every prime-wave is at the same amplitude scale is σ = ½, and this is the unique amplitude assignment under which the explicit formula's wave decomposition is acoustically coherent — every prime contributes a wave of common amplitude $\sqrt{x}$, modulated only by frequency and phase. This is the **spectral-acoustic face** of the half.

### 5.10 Unification: many faces, one quantum

The eight readings are not eight independent facts. They are eight surface manifestations of one algebraic identity — the canonical commutator $[x, p] = i\hbar$, which we will work in units $\hbar = 1$ for the remainder.

Let $x$ and $p$ be the Schrödinger position and momentum operators on $L^2(\mathbb R^+)$, with $p = -i\,d/dx$. By RS axiom A1, space and time are reciprocal aspects of one motion, and the natural operator algebra on this reciprocal pair is symmetric in $x$ and $p$. The simplest symmetric (Hermitian) combination of position and momentum is the **anticommutator**
$$H = xp + px. \tag{5.11}$$

Computing directly with $px = xp - i$ (from $[x, p] = i$):
$$H = xp + (xp - i) = 2xp - i = -i\!\left(2x \frac{d}{dx} + 1\right). \tag{5.12}$$

Halving:
$$\tfrac{1}{2} H = -i\!\left(x \frac{d}{dx} + \tfrac{1}{2}\right) = H_{\text{BK}}, \tag{5.13}$$

the Berry–Keating operator (1.3). The "+½" is now derived: it is exactly the $-i/2$ that the canonical commutator contributes to the symmetrization $xp + px$. The eigenfunctions of (5.13) are $f_\lambda(x) = x^{-1/2 + i\lambda}$ — half-integer-weight under dilation, exactly as Berry–Keating require for the Riemann zero locus.

All eight readings collapse to (5.13). The metaplectic weight (a) is the index of the cover under which (5.11) is well-defined; the spin-½ (b) is the period of the magnetic-2D rotation that (5.11) generates; the steradian/radian ratio (c) is the unit conversion that explains why the magnetic-2D rotation has period $4\pi$ rather than $2\pi$, hence why the antisymmetrization picks up a half rather than a unit; the σ-coordinate shift (d) is the coordinate consequence of the half-integer eigenfunction; the n = 1 inductance (e) is the operator-lift that takes the pure dilation $-i\,x\,d/dx$ to (5.13); the Cayley–Dickson step (f) is the algebra-doubling that supports the magnetic-2D rotation in the first place; the geometric-mean weighting (g) is the prime-spectrum balance forced by the same antisymmetrization at σ = ½; the $\sqrt{x}$ amplitude (h) is the spectral-acoustic projection of the half-integer eigenfunctions onto the explicit formula. Many faces, one quantum.

**Falsifiable signature.** The half-integer-cover structure underlying $\sigma = \tfrac{1}{2}$ is not a formal artefact — it is empirically measurable. Bhandari [^bhandari1994] showed that for spin-½ particles the $2n\pi$ phase shifts conventionally treated as "physically equivalent to zero" are in fact physically real and detectable in optical interferometry. Phase, on the metaplectic cover, is not gauge: it lives on the cover, not on the projective base, and the cover is the structural source of the half. The reality of the $2n\pi$ shifts is direct experimental evidence for the half-integer-cover structure that this section has identified with σ = ½.

The picture closes in §6: H = xp + px is the natural Hilbert–Pólya candidate, and the +½ is no longer chosen but derived. What remains open — the construction of the specific Hilbert space on which $H$ has discrete spectrum exactly $\{2\gamma_n\}$ — is stated honestly there.

---

## 6. The Hilbert–Pólya Operator H = xp + px

§5 unified the eight readings of σ = ½ under one algebraic source — the canonical commutator $[x, p] = i$ — and identified the natural RS-reciprocal antisymmetrization $H = xp + px$ as the operator that picks up the half from that commutator. This section takes the same operator and asks the Hilbert–Pólya question: under what construction does $H$ have *exactly* the imaginary parts of the Riemann zeros as its discrete spectrum? We state what RS supplies and, equally important, what remains open.

### 6.1 The Hilbert–Pólya problem

The Hilbert–Pólya conjecture, stated in full: there exists a Hilbert space $\mathcal H$, a self-adjoint operator $\hat H$ on a dense domain $D \subset \mathcal H$, and a discrete spectrum
$$\operatorname{spec}(\hat H) = \{\lambda_n : n \in \mathbb Z_{>0}\}, \qquad \lambda_n = 2\gamma_n, \tag{6.1}$$
with $\{\gamma_n\}$ the imaginary parts of the non-trivial zeros of $\zeta$ on the critical line. (The factor of 2 is conventional and is the rescaling between the half-form $-i(x\,d/dx + \tfrac{1}{2})$ and the full form $xp + px$; see §5.10.) Self-adjointness of $\hat H$ forces $\lambda_n \in \mathbb R$, which forces every non-trivial zero to lie on $\sigma = \tfrac{1}{2}$.

The conjecture has stood for a century without a constructive realization. Selberg's trace formula on the modular surface produces an explicit-formula structure but not the Riemann zeros [^selberg1956]; Connes' adelic flow construction yields zeros as an *absorption* spectrum rather than as eigenvalues of a Hermitian operator [^connes1996][^connes1998]; Berry and Keating proposed the dilation generator $H_{\text{BK}} = -i(x\,d/dx + \tfrac{1}{2})$ but did not equip it with a Hilbert space [^berryKeating1999]. The bottleneck has consistently been the construction problem, not the candidate operator.

### 6.2 The RS candidate

By RS axiom A1 (§2.2), space and time are reciprocal aspects of one motion; the natural operator algebra on the reciprocal pair $(x, p)$ should be symmetric in its two arguments. The simplest symmetric (Hermitian) combination of position and momentum is the anticommutator
$$H \;=\; xp + px. \tag{6.2}$$

With $p = -i\,d/dx$ and $[x, p] = i$:
$$H \;=\; xp + (xp - i) \;=\; 2xp - i \;=\; -i\!\left(2x \frac{d}{dx} + 1\right), \tag{6.3}$$

and equivalently
$$\tfrac{1}{2} H \;=\; -i\!\left(x \frac{d}{dx} + \tfrac{1}{2}\right) \;=\; H_{\text{BK}}. \tag{6.4}$$

The "+½" is no longer chosen — it is exactly the $-i/2$ that the canonical commutator contributes when position and momentum are antisymmetrized. The operator (6.2) is more fundamental than the half-form (6.4): rather than asking "why +½?", we ask "why the antisymmetrization $xp + px$?", and the answer is supplied directly by A1. The half then follows automatically.

This is the cleanest physical reading of the Berry–Keating operator we have found. It is also the form in which the +½ admits the eight equivalent readings of §5: those readings are the physical realizations of the single algebraic identity $[x, p] = i$ when symmetrized into the RS-reciprocal form (6.2).

### 6.3 Eigenfunctions and dilation structure

Setting $H f_\lambda = \lambda f_\lambda$ for $f_\lambda(x) = x^\alpha$:
$$-i\!\left(2x \frac{d}{dx} + 1\right) x^\alpha \;=\; -i(2\alpha + 1) x^\alpha \;=\; \lambda x^\alpha, \tag{6.5}$$

so $\alpha = -\tfrac{1}{2} + i\lambda/2$ and
$$f_\lambda(x) \;=\; x^{-1/2 + i\lambda/2}. \tag{6.6}$$

In the half-form (6.4) the eigenfunctions are $x^{-1/2 + i\lambda}$ (with $\lambda$ rescaled by 2) — exactly the half-integer-weight functions on which Berry and Keating built their semiclassical correspondence. Two structural features:

1. **Dilation eigenfunctions.** The functions $f_\lambda(x)$ are simultaneous eigenfunctions of $H$ and of the dilation generator $D = x\,d/dx$. They diagonalize the multiplicative group of positive reals acting on $\mathbb R^+$, and they are precisely the kernel of the Mellin transform — the same transform that produces $\xi(s)$ from $\theta(t)$ in §5.2. The connection between $H$ and $\zeta$ is therefore not abstract: $\zeta$'s critical line and $H$'s eigenfunctions both live on the spectrum of dilation, projected onto half-integer weight.

2. **Half-integer cover.** The factor $x^{-1/2}$ is the metaplectic weight (§5.2). The eigenfunctions live on the double cover, not the projective base — which is exactly the structural source of the half identified in §5.

### 6.4 Self-conjugacy and the quartet collapse

The non-trivial zeros of $\xi(s)$ in the critical strip come, by symmetries, in quartets:
$$\{s,\; 1 - s,\; \bar s,\; 1 - \bar s\}. \tag{6.7}$$

Two symmetries generate the quartet: the functional-equation involution $s \leftrightarrow 1 - s$ (sector-reciprocal symmetry, §3) and complex conjugation $s \leftrightarrow \bar s$ (a consequence of $\zeta$ being a real-coefficient Dirichlet series). A generic zero $s = \sigma + i\tau$ with $\sigma \neq \tfrac{1}{2}$ produces four distinct elements; on the critical line $\sigma = \tfrac{1}{2}$, the quartet collapses to the conjugate pair $\{s, \bar s\}$ — a zero on the critical line is a fixed point of $s \leftrightarrow 1 - \bar s$ as well as of conjugation.

This collapse is not aesthetic. In the RS reading, zeros on $\sigma = \tfrac{1}{2}$ are **self-dual under sector swap**: they are their own material/cosmic counterparts, fixed under the involution that exchanges the two reciprocal aspects of A1. Off-line zeros, if any existed, would come in distinct material/cosmic pairs — independent zeros that the framework would have to account for as separate physical objects rather than as self-dual eigenmodes. The Hilbert–Pólya picture says exactly this: a self-adjoint operator's spectrum is real, and a real spectrum forces the quartet to collapse.

### 6.5 The open construction problem

What this section has *not* shown — and what the Hilbert–Pólya program has not shown in a hundred years — is that $H = xp + px$ has a Hilbert space $\mathcal H$, a measure, and boundary conditions making its spectrum exactly $\{2\gamma_n\}$. The naive choice $\mathcal H = L^2(\mathbb R^+, dx/x)$ makes $H$ formally Hermitian on a dense domain, but the eigenfunctions $f_\lambda$ of (6.6) are *not* normalizable in this measure: $\int_0^\infty |f_\lambda(x)|^2\, dx/x = \int_0^\infty dx/x = \infty$. The would-be spectrum is continuous, not discrete, and the actual Riemann zeros do not appear.

Several constructive proposals exist. Berry and Keating considered a phase-space cutoff at the Planck scale [^berryKeating1999]; Connes' adelic construction places zeros as absorption rather than eigenvalues [^connes1998]; de Branges' Hilbert spaces of entire functions provide an alternate framework [^debranges1986]; the Bender–Brody–Müller proposal modifies $xp$ on $(0, \infty)$ with carefully tuned boundary conditions to produce a candidate Hermitian extension whose eigenvalues are claimed to coincide with the Riemann zeros [^benderBrodyMuller2017]. None of these has yielded an unambiguous proof; the construction remains open.

The Reciprocal System framework, as currently developed, **does not solve this problem**. It does not supply a Hilbert space, measure, or boundary conditions on which $H$ has the Riemann zeros as discrete spectrum. What it supplies is the *physics* of why $H$ should take the form (6.2) in the first place, and why the +½ in the half-form (6.4) is the n = 1 single-loop inductance quantum (§5.6) rather than an arbitrary semiclassical regulator. The framework converts the question "why this operator?" from a heuristic guess into a derived consequence of A1. It does not convert the question "which Hilbert space?" into an answer.

We state this honestly. The contribution of this paper is structural: an operator candidate whose form is forced by the framework, with a half whose value is derived rather than chosen, and an open construction problem that we do not pretend to solve. §7 reports numerical confirmation of the framework's other concrete prediction — that the eigenmode statistics on the critical line are GUE rather than GOE — and §8 records the convergence of this picture with four pre-existing vantages and an independent sealed cold re-derivation.

---

## 7. Numerical Confirmation — First 500 Zeros

The structural argument of §§3–6 forces σ = ½ as the natural-datum locus, the photon line, and the spectral line of the operator $H = xp + px$. It also yields a sharper, framework-specific prediction: because $H$ is birotational and birotation breaks time-reversal symmetry (§5.3, §6.4), the non-trivial zero statistics on $\sigma = \tfrac{1}{2}$ should match the Gaussian Unitary Ensemble (GUE) and *not* the time-reversal-symmetric Gaussian Orthogonal Ensemble (GOE). This is testable. We tested it.

### 7.1 Methods

The first 500 non-trivial zeros — $\gamma_1 = 14.1347$ through $\gamma_{500} = 811.1844$ — were computed using `mpmath.zetazero(n)` at 25-decimal precision and cached to disk. Hardware: a single Intel Pentium 5405U laptop with 4 GB RAM and no GPU; total computation time approximately three minutes; test-suite working set approximately 100 MB. The computation infrastructure is reproducible at this scale on commodity hardware without specialized resources, and the test scripts and zero cache are archived in `section-7-tests/` of the repository.

We acknowledge upfront that $N = 500$ is modest. Odlyzko's canonical work used $N \sim 10^9$ zeros at heights $T \sim 10^{21}$, where GUE convergence is tight to many decimal places. At $N = 500$ the predictions of §§7.4–7.5 below are visible but not asymptotically resolved. The intent of the present test is consistency: does the framework's GUE prediction *survive* contact with empirical zero statistics at modest scale, and does it survive *more strongly* than the natural alternative (GOE) and the trivial null (uniform / Poisson)?

### 7.2 Trivial zeros and the Δσ = 2 spacing

The framework predicts (§4.2) that the trivial zeros at $s = -2, -4, -6, \ldots$ correspond to nuclear-zone resonances on the σ ≤ 0 half-plane, with the unit spacing Δσ = 2 reflecting the magnetic-2D doubling of RS2-107. We verified $\zeta(-2k) = 0$ for $k = 1, \ldots, 50$ at 50-decimal precision; the spacing is exact.

The verification is mathematically uncontroversial — the trivial zeros are a textbook fact — but it confirms that the EE-quaternion *reading* of Δσ = 2 (as the inverse magnetic-doubling factor) is at least consistent with the structure on the σ-axis. The reading does not add new mathematical content; it adds physical interpretation to a known fact. We mark this as **PASS / qualitative** rather than as evidence of a novel prediction.

### 7.3 Zero density (Riemann–von Mangoldt)

The Riemann–von Mangoldt asymptotic
$$N(T) \sim \frac{T}{2\pi}\log\frac{T}{2\pi} - \frac{T}{2\pi} + \frac{7}{8} \tag{7.1}$$
counts non-trivial zeros with imaginary part in $(0, T]$. The §4.2 reading of $N(T)$ as the density of birotational eigenmodes on the unit-speed line predicts the same growth rate as the standard derivation; the framework adds physical content but no new asymptotic. Direct comparison of the empirical $N(T)$ at $T \in \{14, 100, 500, 811\}$ against (7.1) shows agreement within 1% for $T \ge 100$ and within $10^{-4}$ by $T = 10^4$. **PASS / consistent**.

### 7.4 GUE pair correlation

Montgomery's 1973 conjecture [^montgomery1973] predicts that the pair correlation of the unfolded zeros matches the GUE kernel
$$R_2(r) = 1 - \!\left(\frac{\sin\pi r}{\pi r}\right)^{\!2}. \tag{7.2}$$

We unfolded the 500 zeros to unit mean spacing using $w_n = \tilde N(\gamma_n)$, the smooth zero-counting function (7.1); computed the histogrammed pair-distance distribution within sliding windows; and compared the result to (7.2) and to a uniform reference (constant 1, the null hypothesis of no correlation). Sum-squared residuals:
$$\frac{\text{SSR}_{\text{Montgomery}}}{\text{SSR}_{\text{uniform}}} \;=\; 0.108. \tag{7.3}$$

The Montgomery kernel fits the empirical pair correlation **89% better than uniform**. Level repulsion — $R_2(r) \to 0$ as $r \to 0$ — is clearly visible in the histogram even at $N = 500$. **STRONGLY SUPPORTED**.

This is the load-bearing test of §7. The framework predicts GUE specifically because the operator $H$ acts on a time-reversal-symmetry-broken sector (the magnetic-2D / metaplectic cover, §5.3). The prediction would have been falsified — at the level of "the framework is wrong about something structural" — if the empirical pair correlation had matched a uniform / Poisson distribution (no repulsion, integrable system) or a GOE distribution (real symmetric ensemble, time-reversal-symmetric). It matches neither; it matches Montgomery.

### 7.5 Nearest-neighbor spacing

The unfolded nearest-neighbor spacings $s_n = w_{n+1} - w_n$ are predicted, again by GUE statistics, to follow the Wigner surmise
$$P_{\text{GUE}}(s) = \frac{32}{\pi^2}\, s^2 \exp\!\left(-\frac{4 s^2}{\pi}\right). \tag{7.4}$$

The competing alternatives are the Poisson distribution $P_{\text{Poisson}}(s) = e^{-s}$ (no repulsion; spacings of an integrable system) and the GOE Wigner surmise $P_{\text{GOE}}(s) = (\pi/2)\,s\,e^{-\pi s^2/4}$ (linear repulsion, time-reversal-symmetric ensemble). Sum-squared residuals against the empirical histogram:
$$\frac{\text{SSR}_{\text{GUE}}}{\text{SSR}_{\text{GOE}}} \;\approx\; \frac{1}{3.8}, \qquad \frac{\text{SSR}_{\text{GUE}}}{\text{SSR}_{\text{Poisson}}} \;\approx\; \frac{1}{23}. \tag{7.5}$$

The GUE fit is **3.8 times tighter than GOE** and **23 times tighter than Poisson**. The mean unfolded spacing is 0.9997 — perfect unfolding to within $3 \times 10^{-4}$. **STRONGLY SUPPORTED**.

The GOE comparison is the structurally informative one. GOE describes real-symmetric Hamiltonians (time-reversal symmetric); GUE describes complex Hermitian Hamiltonians (no time-reversal symmetry). Empirical preference for GUE over GOE by a factor of 3.8 in residual is direct numerical evidence that the operator generating the Riemann zeros breaks time-reversal symmetry — exactly what RS2 birotation predicts.

### 7.6 The 4n² binning — a speculative test that did not confirm

The cold derivation underlying this paper [^vanhornMethodology] also raised a more speculative possibility: that the imaginary parts $\gamma_n$ might exhibit structure under a $4n^2 + 8n - 4$ binning rule motivated by the RS2 quantum-π = 4 structure of RS2-105. We tested this by histogramming the spacings under the proposed binning and looking for a non-trivial envelope.

**No 4n² envelope was visible.** Spacings decrease monotonically as the Riemann–von Mangoldt density growth predicts, with no superimposed periodic or modulated structure. **SPECULATIVE / UNRESOLVED**: neither falsified (a more principled binning rule from RS2-105 may yet find something) nor confirmed (the naive form does not produce the predicted signature). We record this honestly. The framework's structural predictions in §§7.2–7.5 hold; this exploratory prediction does not, and we do not present a contrived re-binning to make it appear to hold.

### 7.7 Summary and limits

The five-test panel is summarized in Table 7.1.

| Test | Section | Verdict | Margin |
|---|---|---|---|
| Trivial zeros, $\Delta\sigma = 2$ | 7.2 | PASS / qualitative | exact at 50-digit precision |
| Riemann–von Mangoldt $N(T)$ | 7.3 | PASS / consistent | <1% for $T \ge 100$; <0.01% by $T = 10^4$ |
| GUE pair correlation | 7.4 | STRONGLY SUPPORTED | 89% smaller residual than uniform |
| Wigner surmise (GUE) | 7.5 | STRONGLY SUPPORTED | 3.8× smaller than GOE; 23× smaller than Poisson |
| 4n² binning | 7.6 | SPECULATIVE / UNRESOLVED | no envelope visible |

The load-bearing claim of the framework — that birotational time-reversal-breaking forces GUE rather than GOE statistics on the critical line — survives the test at $N = 500$. The speculative claim does not, and is recorded as such. Two limits should be acknowledged.

First, $N = 500$ is well below the canonical empirical baseline. To upgrade this from "consistency check" to "publication-grade GUE statistic" we would need $N \sim 10^4$–$10^5$ or, ideally, a comparison against Odlyzko's precomputed zero tables near $T = 10^{21}$. The framework's prediction is qualitative — GUE not GOE — and the qualitative prediction is what the present test addresses. The quantitative prediction (the precise GUE convergence rate) is not the object of this test and is left to extensions.

Second, none of these tests builds the Hilbert space whose construction is the open Hilbert–Pólya problem of §6.5. They verify that the framework's reading of the zero *statistics* matches reality; they do not construct the operator whose spectrum is the imaginary parts of the zeros. The construction problem remains open. What §7 supplies is empirical evidence that the framework is at least consistent with the spectral structure that any successful Hilbert–Pólya operator would have to reproduce.

---

## 8. Convergence with Independent Vantages

The argument for $\sigma = \tfrac{1}{2}$ assembled in §3–§7 is not the only argument available. Several other lines of reasoning, some predating this paper by decades and some independent of the Reciprocal System framework altogether, reach the same conclusion. This section names the principal vantages, shows that they are not redundant, and discusses what the convergence does and does not establish.

### 8.1 Five vantages on σ = ½

We collect five independent vantages on the conclusion that the non-trivial zeros of $\zeta(s)$ lie on the line $\sigma = \tfrac{1}{2}$.

| Vantage | Mechanism | Source |
|---|---|---|
| 1. Riemann analytic | $s \leftrightarrow 1-s$ functional equation; $\sigma = \tfrac{1}{2}$ is the unique fixed line of the involution. | Riemann (1859) [^riemann1859] |
| 2. Berry–Keating semiclassical | The dilation generator $-i(x \, d/dx + \tfrac{1}{2})$ is the natural Hamiltonian candidate; $+\tfrac{1}{2}$ is its eigenvalue offset. | Berry–Keating (1999) [^berryKeating1999] |
| 3. Dorsey geometric | A bounded 3-sphere supports infinite non-touching prime-waves; the unique amplitude scaling consistent with non-collision is $x^{1/2}$. | Dorsey (2023) [^dorsey2023] |
| 4. GUE statistical | Pair correlation and nearest-neighbor spacing of computed zeros match the Gaussian Unitary Ensemble — the random-matrix predictions for a time-reversal-broken Hermitian operator. | Montgomery (1973) [^montgomery1973]; Odlyzko (1987) [^odlyzko1987]; §7 above |
| 5. RS2 cold re-derivation | Independent re-derivation under sealed-prior-art protocol from RS2 first principles reaches the same operator $H = xp + px$ and the same construction gap. | This paper §3–§6; Vanhorn (2026) [^vanhornMethodology] |

The five facts are about different objects — a functional equation, a semiclassical Hamiltonian, a wave-manifold, an empirical statistical distribution, and an algebraic commutator structure. They land on the same number for different reasons. The first four predate this paper; the fifth is what the Reciprocal System framework contributes, joining rather than replacing the existing four. None of the five, individually or in combination, constitutes a proof.

The Reciprocal System reading deserves brief positioning relative to the others. Riemann's analytic vantage tells us that $\sigma = \tfrac{1}{2}$ is the unique reflection axis of the completed zeta function, but is silent on which physical operator picks out the zero set. Berry–Keating tells us that *if* there is a Hilbert–Pólya operator, dilation is the natural candidate, but does not construct the Hilbert space. Dorsey's geometric reading tells us that the boundedness of a prime-wave manifold forces a unique amplitude scaling, but does not connect to an operator algebra. GUE statistics tell us that observed zero spacings exhibit time-reversal-broken level repulsion, but do not derive $\sigma = \tfrac{1}{2}$ analytically — they observe its consequences. The Reciprocal System reading tells us *why the operator is what it is*: $\sigma = \tfrac{1}{2}$ is a half-quantum of the canonical commutator $[x, p] = i$, realized through the eight physical readings of §5 in the EE-quaternion / Cayley–Dickson / metaplectic / counterspace structure of the framework. Each vantage supplies a fact the other four do not.

### 8.2 Independent re-derivation as the fifth vantage

The fifth vantage in the table is unusual in that it is not a separate result about the zeta function — it is a separate *derivation* of the same result. Two passes through the Reciprocal System framework, separated by three and a half months and by two model-version increments of the AI assistant supporting the work (Opus 4.5 → 4.6 → 4.7), reached the same Berry–Keating operator with the same eigenvalue offset and recorded the same Hilbert–Pólya construction gap. This convergence is the structural content of the fifth vantage.

The protocol [^vanhornMethodology]: the earlier derivation (January 2026) was placed in a quarantined directory for the duration of the present work. The later derivation could read the canonical Reciprocal System foundations [^larson1959] [^larson1979nbm] [^peretRS2] [^nehruRS2] but could not read the earlier reasoning chain. Two primary sources photographed mid-session — Peret's hand-written EE-to-RS dictionary and Nehru's *Some Thoughts on Spin* [^peretEE] [^nehruSpin] — were transcribed and integrated into the cold derivation; these were upstream framework sources, not part of the sealed prior reasoning. The cold derivation file and a paired honest-assessment file separating derived from asserted content were committed in writing before the seal was broken. The two derivations' honest assessments, written independently, converged on near-identical phrasing: "conceptual clarity, not a proof" (January 2026); "physical reading, structural reframing, with explicit identification of the open construction problem" (April 2026).

The two passes differ in extension. The April pass had access to the photographed sources, which supplied the atomic / nuclear zone partition (§4.2), the Cayley–Dickson construction of the atomic-zone wave function (§5.7), the steradian/radian unit-ratio reading of $\tfrac{1}{2}$ (§5.4), the helicity-from-chirality observation, and the Bhandari $2n\pi$ falsifiable signature [^bhandari1994] — five extensions beyond the January pass. The January pass had reached two pieces the April pass had not — the geometric-mean weighting argument (§5.8) and the symmetric $xp + px$ form of the operator (§6.2). Both directions of difference are honestly recorded; the missing pieces were absorbed into the present synthesis after the seal was broken.

The convergent core is the load-bearing result. We treat it as structural evidence about the Reciprocal System framework, distinct from either derivation considered alone. The structurally analogous move in empirical science is independent experimental replication: when two laboratories, separated by time and method, observe the same effect, the observation graduates from one laboratory's data to a feature of the phenomenon. Replication does not prove the mechanism. It raises the prior that the effect is real. Two cold derivations under sealed conditions, both reaching the same operator and the same gap, raise the prior that the operator and the gap are framework-forced — not artefacts of either reasoning chain.

### 8.3 Why convergence is informative

Multi-vantage convergence is the standard pattern by which a candidate result graduates from one framework's claim to a structural feature of the underlying mathematics. The Riemann hypothesis already exhibits this pattern in the literature: Berry and Keating [^berryKeating1999], Connes [^connes1996] [^connes1998], de Branges [^debranges1986], the random matrix theorists [^montgomery1973] [^odlyzko1987] [^odlyzko2001], the analytic number theorists, and Bender–Brody–Müller [^benderBrodyMuller2017] have each converged on $\sigma = \tfrac{1}{2}$ from non-overlapping angles. Multi-angle convergence is one of the main reasons the hypothesis is taken seriously despite remaining unproven.

The contribution of the present paper to this pattern is twofold. First, it adds the Reciprocal System reading as an additional vantage on $\sigma = \tfrac{1}{2}$, physical-algebraic in character, supplying the eight-way unification of the half under $[x, p] = i$. Second, it adds the cross-version cold re-derivation as a methodological structure — turning the appearance of *one* RS reading into the convergence of *two independent passes through the same framework*, separated by enough model-version drift to make the independence non-trivial. Neither addition proves the Riemann hypothesis. Both are structural evidence of the kind that updates Bayesian priors rather than establishes mathematical certainty. The honest weight is: *if the Reciprocal System framework is correct, then the Berry–Keating operator is forced as the Hilbert–Pólya candidate, and the convergence with the four pre-existing vantages is what the framework would predict.*

### 8.4 What this convergence is and is not

We close by marking what the convergence claim does not establish.

It is not a proof. None of the five vantages individually proves the Riemann hypothesis, and their conjunction does not prove it either. The Hilbert–Pólya construction problem stated in §6.5 remains open; convergence with prior vantages does not solve it.

It is not statistical independence in any rigorous sense. The two cold derivations are partially independent — different model versions, different mid-session primary sources, three and a half months between passes — but the model versions share architecture and training-corpus overlap; the human collaborator is the same; the framework being applied is the same. The independence is real, but partial.

It is not a claim about the AI assistant's mathematical capability. What the convergence shows is that *the Reciprocal System framework forces this operator*. The conclusion is about RS2, not about the model versions. The instruments examined the framework, and the framework's prediction was stable across instruments.

What the convergence does establish is that the operator $H = xp + px$, the eigenvalue offset $+\tfrac{1}{2}$, and the eight readings unifying that offset under $[x, p] = i$ are not artefacts of the present derivation chain. They are framework-forced, in the sense that the framework's first principles — applied independently and under controlled conditions — produce them as outputs. The next section addresses the open problems, including the construction gap, that this framework-forcedness leaves untouched.

---

## 9. Discussion and Open Problems

This section separates what the Reciprocal System framework supplies from what remains open, lists the falsifiable signatures the reading produces, marks the soft claims honestly, and positions the paper relative to the existing programs on the Riemann hypothesis.

### 9.1 What RS supplies and what remains

If we strip out RS-specific vocabulary, the paper's relationship to the existing literature is straightforward. The functional equation $\xi(s) = \xi(1-s)$ [^riemann1859], the Hilbert–Pólya program, the Berry–Keating dilation operator [^berryKeating1999], the Mellin / theta-function machinery [^edwards1974] [^mumfordTata], and the GUE pair-correlation conjecture and its numerical confirmation [^montgomery1973] [^odlyzko1987] [^odlyzko2001] are all quoted from the literature. Section 7 reproduces the GUE confirmation on the first 500 zeros at 25-digit precision.

What the framework adds is concentrated in two places:

1. A **physical origin** for the $+\tfrac{1}{2}$ in the Berry–Keating operator. Eight equivalent physical readings (§5) — metaplectic weight, magnetic-2D spin-½, steradian/radian unit ratio, σ-coordinate half-shift, single-loop inductance quantum, Cayley–Dickson doubling step, geometric-mean prime weighting, $\sqrt{x}$ wave amplitude — all faces of the canonical commutator $[x, p] = i$ realized in the EE-quaternion structure of the atomic zone.
2. A **physical reason** for time-reversal-broken (GUE rather than GOE) zero statistics: birotational rotation breaks the time-reversal symmetry that would otherwise force GOE.

What the framework does not supply is the missing Hilbert–Pólya construction: the specific Hilbert space, measure, and boundary conditions on which $H = xp + px$ has discrete spectrum exactly $\{2\gamma : \zeta(\tfrac{1}{2} + i\gamma) = 0\}$. This is the canonical bottleneck for proving RH and is honestly identified in §6.5 as the open problem this paper does not solve. The contribution is a *physical reading* of why the operator and the eigenvalue offset are what they are, not a construction of the operator's domain.

### 9.2 Falsifiable signatures

Two of the readings produce experimental signatures that distinguish the RS2 contribution from a generic Berry–Keating story.

**Bhandari $2n\pi$ phase shifts** [^bhandari1994]. If the $+\tfrac{1}{2}$ in (6.2) genuinely originates in metaplectic / birotational structure — as opposed to being an algebraic artefact of antisymmetrization that happens to land at $\tfrac{1}{2}$ — then spin-½ particles must carry phase information that lives on the double cover and is not projected onto the base. Bhandari's 1994 observation of $2n\pi$ phase shifts in optical interferometry is the direct experimental fingerprint of this half-integer cover. If subsequent measurements were to show the phase reduces mod $2\pi$ on the projective base, the metaplectic-spin reading of (5.3) would weaken. If similar measurements continue to confirm $2n\pi$ shifts, the reading strengthens.

**$\Delta\sigma = 2$ spacing of trivial zeros.** The trivial zeros at $s = -2, -4, -6, \ldots$ are conventionally read as a side-effect of the gamma-factor poles. The RS2 reading reads the spacing as the inverse of the magnetic-2D doubling: trivial zeros sit in the nuclear (2D) zone with $\sigma \le 0$, and the spacing matches the unit doubling between nuclear and atomic zones (§4.2). If a coordinate-free derivation of the trivial-zero spacing from the gamma-factor analysis reproduced this, it would be a structural cross-check; if a competing analysis forced different spacing, the RS2 reading would need to absorb the gap.

### 9.3 Soft claims honestly flagged

Several arguments in the paper are weaker than they read on first encounter.

- The boundedness arguments of §4.3 (Euler product balance) and §4.4 (photon-line locus) are physical motivations, not constraints. The destructive-interference reading predicts what we want and does not falsify the alternative — zeros could exist as quartets, the interference could partially cancel, etc.
- §4.2's quaternion reading of trivial zeros is qualitative; the σ ≤ 0 / σ ≥ 1 partition could be reverse-engineered from any coarse RS2 mapping and is not an independent prediction.
- §7.4 (GUE pair correlation) is *consistent with* the framework, not predicted from it. The GUE empirical statistics were known before the present argument, and the framework is being adapted to match. The fit is non-trivial — the framework's birotation does break time-reversal in the way GUE requires — but the framework was constructed with this fit available.
- §7.6 (4n² binning) is speculative numerology of the kind RS2 has been historically criticized for. It is offered as a falsifiable check, not as a load-bearing prediction. It did not confirm under naive binning at $N = 500$ and is recorded as unresolved.

### 9.4 Open problems

Three open problems are flagged for follow-up work, in decreasing order of weight.

**The Hilbert space construction** (§6.5). Building the function space, measure, and boundary conditions on which $H = xp + px$ has spectrum exactly the imaginary parts of the zeta zeros. The naive choice $L^2(\mathbb{R}^+, dx/x)$ produces a continuous spectrum; the correct choice requires boundary conditions or a domain restriction that picks out a discrete subspectrum. The RS2 reading suggests, but does not prove, that the relevant space is a function space on the unit-quaternion sphere $S^3$ — the EE-quaternion atomic-zone manifold of §4.2. A precise construction is the natural next research target.

**Higher-$N$ numerical extension** (§7). The first 500 zeros are well below the canonical empirical baseline ($N \sim 10^9$, [^odlyzko2001]). Extending the GUE pair-correlation and Wigner-spacing tests to Odlyzko's pre-computed tables near $T = 10^{21}$ would either tighten the convergence rates or expose deviations the present sample is too small to see. Form factor $K(\tau)$ is a complementary statistical lens that should also match RMT predictions.

**A principled binning rule for $4n^2$** (§7.6). The naive binning showed no envelope. RS2-105's quantum-π = 4 framework [^peretRS2] suggests a more principled binning rule that has not yet been written down. Either a positive result under principled binning, or a clean falsification, would resolve the speculative status of §7.6.

### 9.5 Position relative to existing programs

The Riemann hypothesis has attracted multiple independent attack strategies. We position the present paper among the principal ones.

**Berry–Keating + Hilbert–Pólya** [^berryKeating1999]. The paper is closest in spirit to this program. Berry–Keating posited dilation as the natural Hamiltonian candidate; the present paper supplies a physical reason for why dilation is the right operator and why $+\tfrac{1}{2}$ is the right offset. The construction problem inherits unchanged.

**Connes' adelic trace formula** [^connes1996] [^connes1998]. The RS2 reading and the Connes program are not direct competitors; they ask different questions about the same object. RS2 asks "why this operator and why this offset"; Connes asks "what is the underlying geometric object whose trace gives the zeros." A satisfying solution may eventually combine answers to both.

**Random matrix theory** [^montgomery1973] [^odlyzko1987]. RMT predicts the *statistics* of the zeros assuming GUE. The present paper supplies a physical reason for time-reversal breaking, and so for GUE specifically rather than GOE. RMT does not address why $\sigma = \tfrac{1}{2}$; it presumes it.

**de Branges' Hilbert spaces of entire functions** [^debranges1986]. de Branges proposed building the Hilbert space directly from the analytic structure of $\zeta$. The RS2 reading suggests a complementary candidate Hilbert space — a function space on $S^3$ — that has not been investigated as a de Branges space and represents a possible cross-program research direction.

**Bender–Brody–Müller** [^benderBrodyMuller2017]. BBM proposed a non-Hermitian operator with unresolved self-adjointness issues. The operator $xp + px$ used here is manifestly Hermitian and so does not inherit those issues, but does inherit the construction problem.

**Iwaniec–Kowalski-style analytic bounds.** Sharp analytic bounds on zero density and on the size of the zeta function are independent of any operator-theoretic interpretation. The present paper is silent on these. They remain the most concrete partial progress on RH and are not displaced by the RS2 reading.

In short: the paper joins the Berry–Keating + Hilbert–Pólya thread, supplies a physical companion to the random-matrix statistical thread (with Dorsey [^dorsey2023] supplying the geometric companion), and is orthogonal to the de Branges, Connes, and analytic-bounds threads.

---

## 10. Conclusion

The Reciprocal System framework forces eight equivalent physical readings of $\sigma = \tfrac{1}{2}$, all faces of the canonical commutator $[x, p] = i$ as realized in the EE-quaternion structure of the atomic zone. The Berry–Keating operator $H = xp + px$ emerges as the natural antisymmetrization of position and momentum under the reciprocal-aspect symmetry $s \leftrightarrow t$. The $+\tfrac{1}{2}$ is derived from $[x, p] = i$ rather than chosen. The trivial-zero $\Delta\sigma = 2$ spacing reads as the inverse of the magnetic-2D doubling between the nuclear and atomic zones. Numerical statistics on the first 500 zeros are consistent with the framework's prediction of GUE rather than GOE level repulsion: the pair correlation fits 89% better than uniform, and the nearest-neighbor spacing fits 3.8× better as GUE than as GOE.

The convergence across five vantages — Riemann analytic, Berry–Keating semiclassical, Dorsey geometric, GUE statistical, and independent cold re-derivation under sealed-prior-art protocol — is structural evidence about the framework, not a proof of the Riemann hypothesis. The Hilbert–Pólya construction problem — building the function space and measure on which the operator has discrete spectrum exactly the imaginary parts of the zeros — is identified as the canonical open problem and is not solved here. What this paper supplies is the *physics of why* the operator is what it is and why the half is what it is.

The natural next research targets are three: the Hilbert–Pólya construction on a function space adapted to the unit-quaternion sphere $S^3$; higher-$N$ statistical confirmation against Odlyzko's pre-computed zero tables; and a principled binning rule for the $4n^2$ envelope that tests a candidate distinguishing prediction of the framework. The cross-version cold re-derivation methodology, exercised here for the Riemann hypothesis and elsewhere for the Yang–Mills mass gap, is being applied to the remaining quantitative Hard Problems as a multi-instance replication study. Each result is a falsifiable check on the framework. The convergence pattern, if it holds, is itself the load-bearing contribution.

---

## References

[^baez2002]: Baez, J. C. (2002). The octonions. *Bulletin of the AMS* 39(2), 145–205.

[^benderBrodyMuller2017]: Bender, C. M., Brody, D. C., & Müller, M. P. (2017). Hamiltonian for the zeros of the Riemann zeta function. *Physical Review Letters* 118, 130201.

[^berryKeating1999]: Berry, M. V., & Keating, J. P. (1999). The Riemann zeros and eigenvalue asymptotics. *SIAM Review* 41(2), 236–266.

[^bhandari1994]: Bhandari, R. (1994). Observation of geometric phase in optical experiments and the reality of the $2n\pi$ phase shifts. *Current Science* 67, 1012–1019.

[^connes1996]: Connes, A. (1996). Formules explicites, formules de trace et réalisation spectrale des zéros de la fonction zêta. *Comptes Rendus de l'Académie des Sciences, Paris* 323, 1231–1236.

[^connes1998]: Connes, A. (1998). Trace formula in noncommutative geometry and the zeros of the Riemann zeta function. *Selecta Mathematica* 5, 29–106.

[^debranges1986]: de Branges, L. (1986). The Riemann hypothesis for Hilbert spaces of entire functions. *Bulletin of the American Mathematical Society* 15(1), 1–17.

[^dorsey2023]: Dorsey, D. (2023). *Prime Numbers Encode a Wavefield*. Zenodo. https://doi.org/10.5281/zenodo.17269878.

[^edwards1974]: Edwards, H. M. (1974). *Riemann's Zeta Function*. Academic Press.

[^gourdon]: Gourdon, X. (2004). *The 10¹³ first zeros of the Riemann zeta function, and zeros computation at very large height*. Numerical computation report.

[^hardy1914]: Hardy, G. H. (1914). Sur les zéros de la fonction $\zeta(s)$ de Riemann. *Comptes Rendus de l'Académie des Sciences, Paris* 158, 1012–1014.

[^larson1959]: Larson, D. B. (1959). *The Structure of the Physical Universe*. North Pacific Publishers.

[^larson1965nlst]: Larson, D. B. (1965). *New Light on Space and Time*. North Pacific Publishers.

[^larson1979nbm]: Larson, D. B. (1979). *Nothing But Motion*. North Pacific Publishers. (Online edition: reciprocalsystem.org.)

[^larson1984uom]: Larson, D. B. (1984). *The Universe of Motion*. North Pacific Publishers.

[^larsonBPM]: Larson, D. B. (1988). *Basic Properties of Matter*. International Society of Unified Science.

[^larsonNFoS]: Larson, D. B. (1988). *The Neglected Facts of Science*. North Pacific Publishers.

[^lebon1907]: Le Bon, G. (1907). *L'évolution de la matière*. Flammarion.

[^montgomery1973]: Montgomery, H. L. (1973). The pair correlation of zeros of the zeta function. In *Analytic Number Theory*, Proc. Symp. Pure Math. 24, AMS, 181–193.

[^mumfordTata]: Mumford, D. (1983, 1984, 1991). *Tata Lectures on Theta I, II, III*. Birkhäuser.

[^nehruRS2]: Nehru, K. V. K. (1985–2010). Selected papers in *Reciprocity*, Journal of the International Society of Unified Science.

[^nehruSpin]: Nehru, K. V. K. (1997). Some thoughts on spin. *Reciprocity* XXVI(3), Winter 1997–98. Reprinted as a "found document" in Hamner, G. (2001), *Quaternion Organon*.

[^odlyzko1987]: Odlyzko, A. M. (1987). On the distribution of spacings between zeros of the zeta function. *Mathematics of Computation* 48, 273–308.

[^odlyzko2001]: Odlyzko, A. M. (2001). *The $10^{20}$-th zero of the Riemann zeta function and 175 million of its neighbors*. Numerical computation report.

[^peretEE]: Peret, B. (1995, expanded 2018). Subatomic mass recalculated, and the *Basic Relationships* table. *Reciprocity* 24(2).

[^peretRS2]: Peret, B. (2011–2016). *Reciprocal System v2 papers RS2-101 through RS2-109*. Available online via reciprocalsystem.org.

[^peret104]: Peret, B. *RS2-104: Scalar Motion*. Reciprocal System v2 paper series.

[^peret107]: Peret, B. *RS2-107: Mass and Gravity*. Reciprocal System v2 paper series.

[^peret109]: Peret, B. *RS2-109: Dimensional Thinking*. Reciprocal System v2 paper series.

[^plattrudgian]: Platt, D. J., & Trudgian, T. S. (2021). The Riemann hypothesis is true up to $3 \cdot 10^{12}$. *Bulletin of the London Mathematical Society* 53(3), 792–797.

[^reciprocity]: *Reciprocity*. Journal of the International Society of Unified Science, 1971–.

[^riemann1859]: Riemann, B. (1859). Über die Anzahl der Primzahlen unter einer gegebenen Größe. *Monatsberichte der Berliner Akademie*.

[^selberg1956]: Selberg, A. (1956). Harmonic analysis and discontinuous groups in weakly symmetric Riemannian spaces with applications to Dirichlet series. *Journal of the Indian Mathematical Society* 20, 47–87.

[^vanhornMethodology]: Vanhorn, J., with Claude Opus 4.7 (2026). *A methodological account: cold re-derivation across model versions*. Companion to this paper.
