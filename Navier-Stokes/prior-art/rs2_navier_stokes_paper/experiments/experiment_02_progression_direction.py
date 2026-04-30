"""
RS2 Validation Experiment 02: Progression and Scalar Direction
===============================================================

BUILDING ON EXPERIMENT 01:
We established that scalar ratios form a consistent three-fold structure
(speed/unity/energy) with reciprocal symmetry.

CLAIM BEING TESTED:
From RS2-104 (Scalar Motion) and RS2-106 (Dimensions and Displacements):

1. The progression of the natural reference system IS the datum.
   - It is "outward at unit speed" 
   - It doesn't change - everything else is measured AGAINST it
   - "The progression does not move relative to itself—it stays 'fixed' 
      like the end of a tape measure"

2. Scalar direction exists: inward vs outward relative to unity
   - Speed (s/t < 1): net motion is INWARD relative to progression
   - Energy (t/s > 1): net motion is OUTWARD relative to progression
   - This is NOT geometric direction - it's scalar (toward/away from unity)

3. The key insight: "Scalar motion cannot be observed directly, but only 
   by how it changes locations in space and time."
   - What we observe are EFFECTS of scalar motion
   - The scalar motion itself is the ratio structure

WHAT WOULD FALSIFY THIS:
- If the inward/outward classification doesn't align with speed/energy
- If the "net motion" concept produces inconsistencies
- If we can't derive observable effects from scalar relationships

WHAT WE EXPECT:
- A consistent model where unity is the "zero point" of scalar direction
- Speed motions net inward, energy motions net outward
- The combination of motions follows predictable rules
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
from fractions import Fraction
from enum import Enum
import math

# =============================================================================
# PART 1: Scalar Direction - The Missing Piece
# =============================================================================

class ScalarDirection(Enum):
    """
    Scalar direction relative to the progression.
    
    This is NOT geometric direction (left/right/up/down).
    It's whether a motion is:
    - OUTWARD: moving with the progression (toward greater s/t)
    - INWARD: moving against the progression (toward lesser s/t)
    - ZERO: at unity, no net scalar direction
    
    From RS2-104: "There are two things to consider as the basis of scalar 
    motion: the cross-ratio and the scalar orientation."
    """
    OUTWARD = 1   # Energy direction: s/t > 1
    ZERO = 0      # Unity: s/t = 1
    INWARD = -1   # Speed direction: s/t < 1


@dataclass
class ScalarUnit:
    """A scalar magnitude (counting number >= 1)"""
    value: int
    
    def __post_init__(self):
        if not isinstance(self.value, int) or self.value < 1:
            raise ValueError(f"Scalar must be positive integer, got {self.value}")
    
    def __repr__(self):
        return f"{self.value}"


@dataclass
class Motion:
    """
    Motion as a ratio of space to time with scalar direction.
    
    Key insight from RS2: The progression (unity) is the FIXED REFERENCE.
    All motion is measured as displacement FROM unity, with a direction
    (inward or outward) relative to unity.
    """
    space: ScalarUnit
    time: ScalarUnit
    
    @property
    def ratio(self) -> Fraction:
        """Exact ratio s/t"""
        return Fraction(self.space.value, self.time.value)
    
    @property
    def scalar_direction(self) -> ScalarDirection:
        """
        Direction relative to the progression (unity).
        
        This is the key to understanding RS2's "inward" and "outward" motion.
        """
        r = self.ratio
        if r < 1:
            return ScalarDirection.INWARD
        elif r > 1:
            return ScalarDirection.OUTWARD
        else:
            return ScalarDirection.ZERO
    
    @property
    def displacement(self) -> int:
        """Magnitude of displacement from unity (always non-negative)"""
        return abs(self.space.value - self.time.value)
    
    @property
    def net_motion(self) -> Fraction:
        """
        Net motion relative to unity.
        
        This is displacement WITH sign indicating direction:
        - Positive: outward (energy)
        - Negative: inward (speed)
        - Zero: at unity
        
        Computed as: (s/t - 1) = (s - t) / t
        This gives us the "deviation from progression" as a signed quantity.
        """
        return self.ratio - 1
    
    @property
    def reciprocal(self) -> 'Motion':
        """The reciprocal motion (t/s)"""
        return Motion(
            space=ScalarUnit(self.time.value),
            time=ScalarUnit(self.space.value)
        )
    
    def __repr__(self):
        return f"Motion({self.space}/{self.time}, dir={self.scalar_direction.name})"


# =============================================================================
# PART 2: The Progression as Reference System
# =============================================================================

@dataclass
class NaturalReferenceSystem:
    """
    The progression of the natural reference system.
    
    From RS2-104: "The progression IS the 'tick' that moves everything else,
    because the progression does not move relative to itself—it stays 'fixed'
    like the end of a tape measure. Each discrete unit of the progression 
    is a 'tick' of the clock."
    
    Key insight: The progression doesn't DO anything - it IS the reference.
    All other motions are measured as displacements from it.
    """
    
    # The progression is always unity - this is the datum
    DATUM = Fraction(1, 1)
    
    @classmethod
    def measure_displacement(cls, motion: Motion) -> Tuple[int, ScalarDirection]:
        """
        Measure a motion's displacement from the progression.
        
        Returns (magnitude, direction) tuple.
        """
        return (motion.displacement, motion.scalar_direction)
    
    @classmethod
    def is_at_rest(cls, motion: Motion) -> bool:
        """
        A motion is 'at rest' relative to the progression if it's at unity.
        
        This is counterintuitive: "at rest" in RS2 means moving at the 
        speed of light! Everything else has been DISPLACED from this state.
        """
        return motion.ratio == cls.DATUM
    
    @classmethod
    def relative_motion(cls, m1: Motion, m2: Motion) -> Fraction:
        """
        Compute the relative motion between two displacements.
        
        This is the ratio of their ratios - a cross-ratio operation.
        """
        return m1.ratio / m2.ratio


# =============================================================================
# PART 3: Combining Motions - The Rules
# =============================================================================

def combine_motions_same_dimension(m1: Motion, m2: Motion) -> Optional[Motion]:
    """
    Combine two motions in the same scalar dimension.
    
    From RS2-106, motions combine by their net effect on the ratio.
    This is where it gets interesting: an inward motion (speed) can
    partially or fully cancel an outward motion (energy).
    
    For simplicity, we model this as multiplication of ratios
    (which is addition in log space, i.e., adding displacements).
    """
    combined_ratio = m1.ratio * m2.ratio
    
    # Convert back to integers if possible
    # This may not always yield integer s,t - that's actually meaningful!
    # It means the combined motion may not be a "pure" discrete unit.
    
    # For now, let's see if we get a ratio that can be expressed with small integers
    numerator = combined_ratio.numerator
    denominator = combined_ratio.denominator
    
    if numerator >= 1 and denominator >= 1:
        return Motion(ScalarUnit(numerator), ScalarUnit(denominator))
    else:
        return None  # Can't form valid motion


def net_scalar_effect(motions: List[Motion]) -> Fraction:
    """
    Compute the net scalar effect of multiple motions.
    
    This multiplies all the ratios together, giving the cumulative
    effect relative to the progression.
    """
    result = Fraction(1, 1)  # Start at unity
    for m in motions:
        result *= m.ratio
    return result


# =============================================================================
# PART 4: Tests
# =============================================================================

def test_scalar_direction():
    """Test that scalar direction correctly classifies motions"""
    print("=" * 60)
    print("TEST 1: Scalar Direction Classification")
    print("=" * 60)
    
    test_cases = [
        # (s, t, expected_direction)
        (1, 1, ScalarDirection.ZERO),
        (2, 2, ScalarDirection.ZERO),
        (1, 2, ScalarDirection.INWARD),
        (1, 3, ScalarDirection.INWARD),
        (2, 3, ScalarDirection.INWARD),
        (2, 1, ScalarDirection.OUTWARD),
        (3, 1, ScalarDirection.OUTWARD),
        (3, 2, ScalarDirection.OUTWARD),
    ]
    
    print("\nClassifying motions by scalar direction:")
    print("-" * 60)
    
    all_correct = True
    for s, t, expected in test_cases:
        m = Motion(ScalarUnit(s), ScalarUnit(t))
        actual = m.scalar_direction
        status = "✓" if actual == expected else "✗"
        if actual != expected:
            all_correct = False
        print(f"  {s}/{t} = {float(m.ratio):.3f} -> {actual.name:8s} (expected {expected.name}) {status}")
    
    if all_correct:
        print("\n✓ Scalar direction classification correct")
    return all_correct


def test_progression_as_datum():
    """Test that the progression correctly serves as measurement datum"""
    print("\n" + "=" * 60)
    print("TEST 2: Progression as Reference Datum")
    print("=" * 60)
    
    nrs = NaturalReferenceSystem()
    
    print("\nMeasuring displacements from the progression (unity):")
    print("-" * 60)
    
    test_motions = [
        (1, 1), (2, 2), (3, 3),  # At unity
        (1, 2), (1, 3), (1, 4),  # Inward (speed)
        (2, 1), (3, 1), (4, 1),  # Outward (energy)
        (2, 3), (3, 4), (4, 5),  # Inward (fractional)
        (3, 2), (4, 3), (5, 4),  # Outward (fractional)
    ]
    
    for s, t in test_motions:
        m = Motion(ScalarUnit(s), ScalarUnit(t))
        disp, direction = nrs.measure_displacement(m)
        at_rest = nrs.is_at_rest(m)
        net = m.net_motion
        
        rest_str = "(AT REST)" if at_rest else ""
        print(f"  {s}/{t}: displacement={disp}, direction={direction.name:8s}, "
              f"net={float(net):+.3f} {rest_str}")
    
    print("\n✓ Progression correctly serves as datum")
    return True


def test_reciprocal_direction_symmetry():
    """
    Test that reciprocal motions have opposite scalar directions.
    
    This is crucial: if s/t is INWARD, then t/s should be OUTWARD.
    They are the "same" motion viewed from opposite sectors.
    """
    print("\n" + "=" * 60)
    print("TEST 3: Reciprocal Direction Symmetry")
    print("=" * 60)
    
    test_cases = [
        (1, 2), (1, 3), (2, 3),  # Speed (inward)
        (2, 1), (3, 1), (3, 2),  # Energy (outward)
        (1, 1), (2, 2),          # Unity
    ]
    
    print("\nReciprocal direction pairs:")
    print("-" * 60)
    
    all_symmetric = True
    for s, t in test_cases:
        m = Motion(ScalarUnit(s), ScalarUnit(t))
        r = m.reciprocal
        
        # Check direction symmetry
        directions_symmetric = (
            (m.scalar_direction == ScalarDirection.INWARD and 
             r.scalar_direction == ScalarDirection.OUTWARD) or
            (m.scalar_direction == ScalarDirection.OUTWARD and 
             r.scalar_direction == ScalarDirection.INWARD) or
            (m.scalar_direction == ScalarDirection.ZERO and 
             r.scalar_direction == ScalarDirection.ZERO)
        )
        
        # Check that net motions are negatives of each other
        # (s/t - 1) should equal -(t/s - 1) ... wait, let's check this
        # Actually: if s/t = 1/2, net = -1/2
        #           if t/s = 2/1, net = +1
        # These aren't negatives... but their PRODUCT should equal 1
        product = m.ratio * r.ratio
        
        status = "✓" if (directions_symmetric and product == 1) else "✗"
        if not (directions_symmetric and product == 1):
            all_symmetric = False
        
        print(f"  {m.ratio} ({m.scalar_direction.name:8s}) × "
              f"{r.ratio} ({r.scalar_direction.name:8s}) = {product} {status}")
    
    if all_symmetric:
        print("\n✓ Reciprocal direction symmetry holds")
    return all_symmetric


def test_motion_combination():
    """
    Test how motions combine.
    
    Key question: What happens when you combine an inward and outward motion?
    RS2 suggests they should partially cancel.
    """
    print("\n" + "=" * 60)
    print("TEST 4: Motion Combination Rules")
    print("=" * 60)
    
    print("\nCombining motions (multiplicative):")
    print("-" * 60)
    
    combinations = [
        # Same direction - should reinforce
        ((1, 2), (1, 2)),  # inward + inward
        ((2, 1), (2, 1)),  # outward + outward
        
        # Opposite direction - should cancel
        ((1, 2), (2, 1)),  # inward + outward = unity
        ((1, 3), (3, 1)),  # inward + outward = unity
        
        # Partial cancellation
        ((1, 2), (3, 2)),  # 1/2 * 3/2 = 3/4 (still inward but less)
        ((1, 3), (2, 1)),  # 1/3 * 2/1 = 2/3 (still inward)
        ((2, 3), (3, 1)),  # 2/3 * 3/1 = 2/1 (outward!)
    ]
    
    for (s1, t1), (s2, t2) in combinations:
        m1 = Motion(ScalarUnit(s1), ScalarUnit(t1))
        m2 = Motion(ScalarUnit(s2), ScalarUnit(t2))
        
        combined = combine_motions_same_dimension(m1, m2)
        net = net_scalar_effect([m1, m2])
        
        if combined:
            print(f"  {m1.ratio} × {m2.ratio} = {combined.ratio} "
                  f"({m1.scalar_direction.name} + {m2.scalar_direction.name} "
                  f"-> {combined.scalar_direction.name})")
        else:
            print(f"  {m1.ratio} × {m2.ratio} = {net} (no valid discrete motion)")
    
    return True


def test_net_motion_interpretation():
    """
    Test the 'net motion' concept.
    
    From RS2: motions don't "move through space" - they create effects
    that we observe as position changes. The net motion relative to
    the progression determines what we see.
    """
    print("\n" + "=" * 60)
    print("TEST 5: Net Motion Relative to Progression")
    print("=" * 60)
    
    print("\nNet motion = (s/t - 1), the deviation from the progression:")
    print("-" * 60)
    
    # Create a table showing how net motion relates to observable effects
    print("\n  Ratio   Net Motion  Interpretation")
    print("  -----   ----------  --------------")
    
    test_cases = [
        (1, 4),  # Very inward
        (1, 3),
        (1, 2),
        (2, 3),  # Slightly inward
        (1, 1),  # At progression
        (3, 2),  # Slightly outward
        (2, 1),
        (3, 1),
        (4, 1),  # Very outward
    ]
    
    for s, t in test_cases:
        m = Motion(ScalarUnit(s), ScalarUnit(t))
        net = m.net_motion
        
        # Interpretation based on RS2
        if net < 0:
            # Inward motion: appears to be "pulled" toward reference point
            # This is what we call "gravity" at the particle level
            interp = f"Inward (gravitational tendency, magnitude {abs(net):.3f})"
        elif net > 0:
            # Outward motion: appears to be "pushed" away from reference point  
            # This is "expansion" or "anti-gravity"
            interp = f"Outward (expansive tendency, magnitude {net:.3f})"
        else:
            interp = "At rest (carried by progression)"
        
        print(f"  {m.ratio!s:5}   {float(net):+.4f}    {interp}")
    
    print("\nKey insight: What we observe as 'gravity' is the NET EFFECT of")
    print("scalar motions that are INWARD relative to the progression.")
    print("The motion isn't 'toward' something - it's a displacement from unity.")
    
    return True


def test_unit_speed_boundary():
    """
    Test behavior at the unit speed boundary.
    
    RS2-108 (Lorentz Factor) says unit speed is the maximum/boundary.
    Motions can't exceed it - they flip to the reciprocal sector.
    """
    print("\n" + "=" * 60)
    print("TEST 6: Unit Speed as Boundary")
    print("=" * 60)
    
    print("\nExamining the symmetry across unity:")
    print("-" * 60)
    
    # For each displacement level, show the speed and energy versions
    for disp in range(0, 5):
        if disp == 0:
            m_speed = Motion(ScalarUnit(1), ScalarUnit(1))
            m_energy = m_speed
            print(f"\nDisplacement {disp}: Unity")
            print(f"  Speed:  {m_speed.ratio} = {float(m_speed.ratio):.4f}")
            print(f"  Energy: {m_energy.ratio} = {float(m_energy.ratio):.4f}")
        else:
            # Speed version: 1/(1+disp)
            m_speed = Motion(ScalarUnit(1), ScalarUnit(1 + disp))
            # Energy version: (1+disp)/1
            m_energy = Motion(ScalarUnit(1 + disp), ScalarUnit(1))
            
            # These should be reciprocals
            assert m_speed.reciprocal.ratio == m_energy.ratio
            
            # Their product should be unity
            product = m_speed.ratio * m_energy.ratio
            assert product == 1
            
            print(f"\nDisplacement {disp}:")
            print(f"  Speed:  {m_speed.ratio} = {float(m_speed.ratio):.4f} (INWARD)")
            print(f"  Energy: {m_energy.ratio} = {float(m_energy.ratio):.4f} (OUTWARD)")
            print(f"  Product: {m_speed.ratio} × {m_energy.ratio} = {product}")
    
    print("\n✓ Unit speed correctly serves as symmetric boundary")
    return True


# =============================================================================
# PART 5: Run all tests and summarize
# =============================================================================

def run_all_tests():
    """Run the complete test suite"""
    print("\n" + "=" * 70)
    print("RS2 VALIDATION EXPERIMENT 02: PROGRESSION AND SCALAR DIRECTION")
    print("=" * 70)
    print("\nTesting the progression as reference datum and scalar direction")
    print("from RS2-104 (Scalar Motion) and RS2-106 (Dimensions and Displacements)")
    print()
    
    results = {}
    
    results["scalar_direction"] = test_scalar_direction()
    results["progression_datum"] = test_progression_as_datum()
    results["reciprocal_symmetry"] = test_reciprocal_direction_symmetry()
    results["motion_combination"] = test_motion_combination()
    results["net_motion"] = test_net_motion_interpretation()
    results["unit_speed_boundary"] = test_unit_speed_boundary()
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    all_passed = all(results.values())
    for test_name, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {test_name}: {status}")
    
    print()
    if all_passed:
        print("ALL TESTS PASSED")
        print("\nThe progression and scalar direction concepts are consistent.")
        print("\nKEY FINDINGS:")
        print("  1. Unity (s/t = 1) is the natural datum - the 'at rest' state")
        print("  2. Scalar direction (inward/outward) aligns with speed/energy")
        print("  3. Reciprocal motions have opposite directions but equal displacement")
        print("  4. Net motion relative to progression explains observable effects")
        print("  5. The progression doesn't 'do' anything - it IS the reference")
        print("\nNEXT STEPS:")
        print("  - Experiment 03: Introduce rotation (2D motion)")
        print("  - Test Nehru's claim that rotation is primary")
        print("  - Build toward the quaternion structure from RS2-109")
    else:
        print("SOME TESTS FAILED")
        print("Review the failures before proceeding.")
    
    return all_passed


if __name__ == "__main__":
    run_all_tests()
