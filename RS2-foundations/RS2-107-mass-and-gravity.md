---
paper: RS2-107
title: Mass and Gravity
author: Bruce Peret
date_authored: 2014 (Rev. 22)
drive_fileId: 1uEJJYSLdEaic_DMSEPcVD_8muptP3lrm
read_date: 2026-04-28
---

# RS2-107: Mass and Gravity

## Core thesis

107 redefines mass via Gustave LeBon's older "weight / velocity" formula (matching the dimensions of Larson's A-B-C atomic notation), shows mass and gravity are reciprocals (mass = outward angular velocity in time over a volume; gravity = inward linear velocity in space at a point), and explains why electrons/positrons/photons/neutrinos are "massless" (their net temporal displacement is below 3, the threshold where ln(Δt) ≥ 1 unit). RS2 promotes angular velocity to a *primary* motion (the "tao of motion": yang = linear, yin = angular), eliminating the need for Larson's direction-reversal/rotational-base machinery.

## Key constructions

- **LeBon mass = weight / velocity.** The 19th-century definition. Dimensionally matches RS atomic structure: t/s² (weight) divided by s/t (velocity) = t³/s³·t.
- **Modern mass = force / acceleration.** Dimensions: t/s³·t² — also recoverable in RS but a worse match for the A-B-C structure.
- **RS mass from atomic structure.** Mass = magnetic rotation × electric rotation = (t/s × t/s) × (s/t) = t³/s³·t. Magnetic rotation gives primary mass; electric rotation gives a small modifier in equivalent time.
- **Mass / gravity reciprocity (full inversion).** Reciprocal of mass:
  - time → space
  - outward → inward
  - angular (circumferential) → linear (radial)
  - volume → point location
  Yields gravity: an inward, linear velocity in space concentrated at a point — the "center of gravity." Mass and gravity are the *same thing from inverse perspectives*.
- **Step measure vs. growth measure (geometric reciprocity).** Step measure: linear count, range 0 → 1, paired with the speed unit of motion. Growth measure: integral / infinitesimal counting, range 1 → ∞, paired with the energy unit. They convert via the natural logarithm: `Δs = ln(Δt)`.
- **Massless threshold (Δt = 3).** Gravity in space appears only when net temporal displacement Δt is large enough to produce ≥ 1 spatial unit: `gravity = floor(ln(Δt))`. Δt ∈ {0, 1, 2} ⇒ no observable mass. Photons, positrons, electrons, neutrinos are below this threshold. Proton (Δt = 3, including rotational base) is the first massive particle.
- **Free-dimension rule for c-progression carry.** A particle is carried at unit speed by the progression iff it has at least one *free* scalar dimension. Charge occupies 2 dimensions; charged particles cannot be carried. Photons (birotation, occupy 2 of 3) are carried in the third dimension. Uncharged electrons (1 dim) leave 2 free.
- **Birotation and Euler's formula.** Vibration is *not* a primary motion (contra Larson). It arises as shear strain between oppositely directed motions: `e^(ix) + e^(-ix) = 2 cos(x)`.
- **Direction reversal and rotational base (Larson construction, retired in RS2).** Larson, with only linear primary motion, needed direction reversal as a "diameter" to spin around (creating the rotational base). RS2, with both linear and angular as primary, makes every location a potential rotational base and direction reversal becomes unnecessary.
- **Tao of motion.** Linear (yang) translation + angular (yin) rotation, both primary. "Spin is yin." This completes the symmetry.
- **Rotational visualizations.**
  - 2D magnetic rotation: cone whose wide end sweeps the surface of a sphere; takes 720° (4π rad). In physics: spin-½ (appears half-speed because it takes two 360° rotations to complete).
  - 1D electric rotation: spinning disc; takes 360° (2π rad). In physics: integer spin (spin-1).
  - 1D vibration: two opposing electric rotations producing a cosine waveform. (Sine traced out by a rod with a flexible elbow rotated opposite-handedly at each end.)
  - 1D rotational vibration = electric charge / electric field.
  - 2D rotational vibration = magnetic charge / magnetism.
- **Equivalent space / equivalent time.** A 2nd unit of motion expresses temporal motion in equivalent space; a 3rd unit cannot, because that "space" is already occupied. The electric rotation modifies mass via "equivalent time" (the reciprocal concept).

## Equations / formal relations

Modern mass:

```
mass = force / acceleration  =  (t/s³) / (s/t²)  =  t³/s³ · t²
                            =  t³ · t² / s³        (Eq. 1)
```

LeBon mass:

```
mass = weight / velocity     =  (t/s²) / (s/t)   =  t³/s³ · t      (Eq. 2)
```

RS mass from atomic A-B-C:

```
mass = magnetic_rotation × electric_rotation
     = (t/s · t/s) · (s/t)                      (Eq. 3)
     = t³/s³ · t
```

Step / growth measure conversion:

```
Δs = ln(Δt)                  (transform growth-measure in equivalent space → step-measure in linear space)
```

Mass-to-gravity threshold:

```
gravity = floor( ln(Δt) )

Δt = 0, 1, 2   ⇒  gravity = 0   (massless: photon, positron, electron, neutrino)
Δt ≥ 3         ⇒  gravity ≥ 1   (proton is first massive particle; ln(3) ≈ 1.1)
```

Birotation as Euler's identity:

```
e^(ix) + e^(-ix) = 2 cos(x)
```

(Two counter-rotations summing to a vibration / cosine in one dimension.)

Rotational sweeps:

```
2D magnetic (solid) rotation : period 4π rad   (spin-½ in conventional physics)
1D electric (disc) rotation  : period 2π rad   (spin-1  in conventional physics)
```

## Verbatim quotes (definitions)

Mass / gravity reciprocity (full block):

> So, we have mass defined as an outward, angular velocity in time, defining a volume. Let's take a complete reciprocal of mass and see what we have as a natural consequence:
>
> • The aspect of time becomes space.
>
> • Outward motion becomes inward motion.
>
> • Angular (circumferential) velocity becomes linear (radial) velocity.
>
> • Volume becomes a point location.
>
> The reciprocal of mass is therefore an inward, linear velocity in space that can be expressed through a single point. That is the definition of gravity, where the "point" is the "center of gravity." Mass and gravity are the same thing, from inverse perspectives.

Step / growth measure:

> Step measure is the conventional method of measuring finite quantities, just like pacing off steps to measure distance. This is associated with the first unit of motion, speed, with the range of 0→1. [...] Growth measure is associated with the second unit of motion, energy, with the range of 1→∞. Since we cannot do a finite count to infinity, growth measure is done with the Calculus concept of infinitesimals, the integral. To transform this growth measure in equivalent space to a step measure in linear space, the natural logarithm must be used: Δs = ln(Δt).

Massless threshold:

> until the magnitude of a temporal rotation, mass, becomes high enough to produce a single unit of inward, spatial magnitude, gravity does not exist in space. And that occurs with a net temporal speed of 3 displacement units, since ln(3) = 1.1.

> If you are a computer/math person, gravity = floor(ln(Δt)). When Δt = 0, 1 or 2, gravity = 0 = massless.

Free-dimension carry rule:

> In order for a particle to be carried, there needs to be a free dimension, a dimension at unit speed in one of the three scalar dimensions of motion for the progression to have effect.

Tao of motion:

> In RS2, the reevaluation of the Reciprocal System, we assume that the yin, angular velocity is a primary motion along with the yang, linear speed, completing the "tao of motion."

> every location is potentially a "rotational base" and the concept of a "direction reversal" is unnecessary, because rotation is primary and RS2 does not require "something to rotate." This infers that the concept of vibration, which Larson associates with his direction reversal, is not a primary motion but only arises as shear strain from oppositely directed motions, such as the counter-rotations of a birotation as expressed in Euler's formula, e^ix + e^-ix = 2 cos(x).

Spin classification:

> 2-dimensional magnetic rotation: a cone with the wide end expanding across the surface of a sphere. This is known as a solid rotation that takes 720 degrees, or 4π radians to complete. In physics, this is measured as a particle with spin-½, because it appears to take two, 360-degree rotations to complete (they assume it is going at half speed).

## Connections to other papers

- 106: Builds directly on the dimensions/displacements (A-B-C, magnetic + electric rotations) and speed-range structure.
- 105: Quantum-π = 4 underlies the 2D magnetic = "area" / 1D electric = "perimeter" identification used here.
- Forward: 108 (Lorentz Factor) addresses time-dilation; 109 (Dimensional Thinking) systematizes the rotational visualizations.
- External: Gustave LeBon, *The Evolution of Forces* (1908). Larson, *Basic Properties of Matter*, ISUS 1988, p. 7 (Eq. 1-1, "Solid Cohesion"). Self-reference to RS2-106.

## Note on Riemann relevance

The "spin-½" identification (2D magnetic / solid rotation, period 4π) appears here. Conventional physics interprets the half because the rotation takes two full 360° passes — an involution after one cycle, identity after two. This is structurally close to a fixed-point-of-an-involution argument, but the paper does not connect it to any zeta-function structure.
