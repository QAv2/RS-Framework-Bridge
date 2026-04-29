---
paper: RS2-108
title: The Lorentz Factor
author: Bruce Peret
date_authored: 2015 (Rev. 26)
drive_fileId: 1tlWjZZMo4Sf28V5n9uRSMeMonrDM8TSy
read_date: 2026-04-28
---

# RS2-108: The Lorentz Factor

## Core thesis

108 reinterprets the Lorentz factor as the equation of a *unit circle* (Pythagorean), with both ± roots considered, in which v is the x-axis and 1/γ the y-axis (mapped to an imaginary axis in the complex plane). The unit circle ranges from v = +1 (outward progression) through v = 0 (direction reversal / birotation, ±γ) to v = -1 (gravity, "speed of light running backwards"). Larson's direction-reversal + rotational-base machinery is shown to be a 1-dimensional kludge for what RS2 upgrades to a 3D quaternion rotation. The c-progression is thus the *fulcrum* between motion in space and motion in time.

## Key constructions

- **Lorentz factor as unit-circle equation.** When v is normalized by c (which is unity in natural units), `γ² = c² - v²` becomes `1 = γ² + v²` — the Pythagorean / unit-circle equation. The square root has two solutions; conventional physics drops the negative.
- **β eliminated.** β = v/c is unnecessary in natural units because c = 1; v is already normalized.
- **Both roots considered.** Honoring the ± sign reveals the negative-velocity branch (v → -1, i.e., the speed of light running backwards = gravity). Suppressing the negative root creates a hidden 0/0 problem ("classic division by zero" allowing 2 = 1 paradox).
- **c-progression as fulcrum.** Unit speed (c) is the boundary between motion-in-space (s/t < 1) and motion-in-time (s/t > 1, equivalently t/s < 1). It is the *maximum* speed of the universe; one can only slow down from it. Photons have no resistance and are carried at this speed.
- **Direction reversal in the complex plane.** Larson's "direction reversal" moves the system from the (+1, 0) progression coordinate to (0, ±γ) — the imaginary axis. This produces a birotation, a counter-rotating pair, whose sum is `(e^(+iγ) + e^(-iγ)) / 2 = cos(γ)`. Larson identifies this cosine as a photon.
- **Rotational base.** Inward scalar rotation moves net motion to (-1, 0) — gravity. The rotational base is therefore the gravitational opposition to the c-progression, at the same magnitude but inward.
- **"Increasing mass" reinterpreted.** Mass remains constant; what increases is the gravitational pull (a "heavier" object resists motion more). The Lorentz Factor is a "kludge hiding the use of imaginary quantities to describe a gravitational field structure."
- **Quaternion upgrade in RS2.** The Lorentz unit circle is 1D, the rotational base is 2D, but RS2 demands 3D rotation, which (per Hamilton) requires 4 quaternion units {1, i, j, k}. RS2 replaces the complex plane with a quaternion, restructuring the photon as a quaternion rotation: 1D electric rotation (k) + 2D magnetic rotation (i·j). Since i·j = k, a birotation can be formed as i·j·(-k).
- **Acceleration ≠ adding velocity.** A rocket does not accelerate by adding speed; it accelerates by *neutralizing the inward motion of gravity*, allowing the object to return to the default speed of unity (c). Particle accelerators reduce resistance; they cannot push past unit speed because there is no further resistance to reduce.
- **Faster-than-light in RS.** Translational velocities are always ≤ unit speed, but other faster-than-light motions exist (manifesting differently from "warp drive").

## Equations / formal relations

Lorentz factor (conventional):

```
γ = 1 / sqrt(1 - (v/c)²)        with β ≡ v/c, so γ = 1 / sqrt(1 - β²)
```

In natural units (c = 1), Peret's reformulation:

```
1/γ = sqrt(c² - v²)             →   sqrt(1 - v²)
(1/γ)² = c² - v²                →   1 = γ² + v²
```

(Gloss: this is a unit circle r² = x² + y² with r = c = 1.)

Right-triangle / unit-circle identification:

```
c² = γ² + v²        →    1 = (1/γ)² + v²
                          (sin² + cos² = 1, with v on x and 1/γ on y)
```

Validity range:

```
-1 ≤ v ≤ +1         (any |v| > 1 is undefined: faster than the universe's max)
```

Birotation as Euler / cosine:

```
( e^(+iγ) + e^(-iγ) ) / 2  =  cos(γ)
```

Quaternion identification of the photon:

```
electric rotation (1D) :  k
magnetic rotation (2D) :  i · j
identity                :  i · j  =  k     (so birotation as i·j·(−k))
```

The three special points on the unit-circle motion diagram:

```
(+1, 0)   :  outward progression of natural reference system (default)
(0, ±γ)   :  direction reversal / birotation / photon (cosine waveform)
(-1, 0)   :  inward rotational base (gravity; "speed of light running backwards")
```

The "2 = 1" paradox derivation footnoted:

```
Let a = b. Then a² = ab; a² + a² = a² + ab; 2a² = a² + ab;
2a² - 2ab = a² + ab - 2ab; 2a² - 2ab = a² - ab; 2(a² - ab) = 1·(a² - ab);
cancel (a² - ab) ⇒ 2 = 1.   (Hidden division by zero.)
```

## Verbatim quotes (definitions)

Lorentz factor as unit circle:

> One should also note that the equation for a right triangle is also the equation for a circle: r² = x² + y², where r is the radius. Because r = c = 1, this is a unit circle, with the velocity on the x axis and the Lorentz factor being the corresponding value on the y axis.

c as fulcrum:

> In the Reciprocal System, the speed of light (unit speed) is the fulcrum between motion in space and motion in time. As such, it is the upper limit of both of those motions, essentially being the maximum speed of the universe, which is referred to as the progression of the natural reference system. You can only slow down from this speed.

Suppressed negative root:

> science also tends to overlook one of the more interesting properties of the square root—that the function returns two solutions, a positive one and a negative one. The negative one is ignored (though the absolute value is never included in the Lorentz equation), because it would indicate that time, length and relativistic mass could also be negative. But if you consider both solutions simultaneously, then a bigger problem arises… they cancel each other out and you end up with the classic "division by zero" problem [...]

Direction reversal in complex plane:

> the gamma function represents the imaginary axis (1/γ = -γ). By default, the Universe is expanding at unit speed, having the coordinates of (+1,0) on the diagram. Larson then introduces the concept of a direction reversal, which results in a linear vibration. This is moving inward (left on the v axis) to the coordinates (0,±1). The progression velocity appears to stop (v=0), but there is now a split across the gamma axis, which is "imaginary" and rotational, creating the two, oppositely-directed rotations that are known as a birotation.

Rotational base = gravity:

> Larson adds an inward scalar rotation to the photon, moving the net motion to the (-1,0) coordinate with a single speed solution, creating the rotational base, whose net motion opposes the progression at the same velocity, the speed of light running backwards that we call gravity, a very familiar concept.

"Lorentz Fudge" assessment:

> Essentially, the Lorentz Factor is just a kludge hiding the use of imaginary quantities to describe a gravitational field structure, in a fashion similar to the imaginary quantities used to describe electric and magnetic fields. This gravitational opposition to the progression is what gives the appearance of increasing mass—even though mass remains constant—since a "heavier" object must have more gravitational pull and be harder to move.

Quaternion upgrade:

> The Lorentz Fudge is a 1-dimensional solution to a 2-dimensional problem, as is Larson's definition of the rotational base. However, the Universe is 3-dimensional and as William Hamilton discovered, it takes 4 dimensions to solve a 3-dimensional rotation: the quaternion.
>
> The RS2 solution was to upgrade the complex plane of the corrected Lorentz Factor and replace it with a quaternion. This, however, changes Larson's 2-unit approach of speed and energy into a 4-unit system of +1, i, i.j and i.j.k=-1. [...] Since i·j = k, a birotation can be formed along electromagnetic lines, using i·j·(-k), providing similar behavior to Prof. KVK Nehru's original birotation model.

The 8-point summary, verbatim items 3-6:

> 3. Unit speed is the maximum speed the physical universe is capable of, expressed in the Reciprocal System as the outward progression of the natural reference system.
>
> 4. The minimum speed is negative unity, the inward motion expressed by gravitation.
>
> 5. The default speed of the Universe is unity. When a conventional object "at rest" is accelerated, what is actually happening is that the inward motion of gravity is being neutralized. A rocket isn't increasing its speed by thrust—the thrust is reducing the effect gravitation is having upon it, allowing it to return to the default speed of unity (the speed of light).
>
> 6. It is impossible to accelerate an object past the speed of light in space, because you are not adding velocity—you are reducing resistance and once that resistance is gone, you are done.

## Connections to other papers

- 105: Quantum-π = 4 vs. analog π is recapped here; the unit-circle of Lorentz uses the analog π.
- 106: Speed ranges and the c-fulcrum interpretation. Direction reversal connects to the displacement notation.
- 107: Mass-gravity reciprocity is reformulated geometrically here as the (+1, 0) ↔ (-1, 0) reflection across v = 0.
- Forward: 109 (Dimensional Thinking) extends the quaternion upgrade. The "i·j·(-k)" birotation model is promised "in a future paper."
- External: KVK Nehru, "The Law of Conservation of Direction," *Reciprocity* 18 № 3 p. 3. Wikipedia (Simple English) for the Lorentz factor formula.

## Riemann-relevance flags

The unit-circle structure with v ∈ [-1, +1] and 1/γ on the imaginary axis is the most directly suggestive structure for a critical-line argument:

- **±γ involution.** The Lorentz factor has a ± symmetry that conventional physics suppresses. RS2 keeps both roots; the symmetry is `γ → -γ` (i.e., reflection across the v-axis, which is reflection across the real axis in the complex plane).
- **Reciprocal involution.** The yin/yang reciprocity s/t ↔ t/s with c as fulcrum is structurally `s ↦ 1/s` reflection, which has fixed points only at |s| = 1.
- **Three points on unit circle.** (+1, 0), (0, ±γ), (-1, 0). The middle point on the imaginary axis is where the cosine waveform / photon lives — this is the structural analogue of the critical line if one identifies "real axis = velocity axis = direction" and "imaginary axis = γ axis."
- **No explicit ½.** The number "½" does not appear in this paper. The spin-½ identification is in 107.
