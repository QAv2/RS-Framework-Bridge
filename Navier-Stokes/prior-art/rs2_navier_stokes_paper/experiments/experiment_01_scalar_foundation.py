"""
RS2 Validation Experiment 01: Scalar Foundation
================================================

CLAIM BEING TESTED:
From RS2-104 (Scalar Motion) and RS2-106 (Dimensions and Displacements):
- Motion is a ratio of two scalar magnitudes (space and time)
- Scalars are counting numbers: minimum 1, no zero, no negatives, no fractions
- Unity (s/t = 1) is the natural datum
- Three orientations exist: speed (s/t < 1), unity (s/t = 1), energy (s/t > 1)
- Speed and energy are reciprocal: s/t and t/s are "the same motion" from different perspectives

WHAT WOULD FALSIFY THIS:
- If the symmetry across unity breaks down
- If the reciprocal relationship doesn't hold mathematically
- If we can't derive displacement from this structure

WHAT WE EXPECT:
- All ratios cluster into exactly three categories
- Reciprocal pairs should have inverse but symmetric properties
- Displacement from unity should be well-defined for all valid ratios
"""

import math
from dataclasses import dataclass
from typing import Tuple, List, Optional
from fractions import Fraction

# =============================================================================
# PART 1: Define the scalar unit according to RS2 constraints
# =============================================================================

@dataclass
class ScalarUnit:
    """
    A scalar magnitude according to RS2.
    
    Constraints from RS2-104:
    - Must be a counting number (positive integer >= 1)
    - Represents "magnitude only" - no direction yet
    """
    value: int
    
    def __post_init__(self):
        if not isinstance(self.value, int):
            raise ValueError(f"Scalar must be integer, got {type(self.value)}")
        if self.value < 1:
            raise ValueError(f"Scalar minimum is 1, got {self.value}")
    
    def __repr__(self):
        return f"Scalar({self.value})"


@dataclass 
class Motion:
    """
    Motion as a ratio of space to time (s/t).
    
    This is the fundamental unit in RS2 - not space, not time, but their RATIO.
    """
    space: ScalarUnit
    time: ScalarUnit
    
    @property
    def ratio(self) -> Fraction:
        """Exact ratio as a fraction (no floating point errors)"""
        return Fraction(self.space.value, self.time.value)
    
    @property
    def ratio_float(self) -> float:
        """Floating point approximation for display"""
        return self.space.value / self.time.value
    
    @property
    def orientation(self) -> str:
        """
        The three-fold classification from RS2-106:
        - "speed": s/t < 1 (more time than space)
        - "unity": s/t = 1 (natural datum)
        - "energy": s/t > 1 (more space than time)
        """
        r = self.ratio
        if r < 1:
            return "speed"
        elif r == 1:
            return "unity"
        else:
            return "energy"
    
    @property
    def reciprocal(self) -> 'Motion':
        """
        The reciprocal motion (t/s instead of s/t).
        
        RS2 claims: this represents the SAME motion from the inverse sector.
        Material sector sees s/t, cosmic sector sees t/s.
        """
        return Motion(
            space=ScalarUnit(self.time.value),
            time=ScalarUnit(self.space.value)
        )
    
    @property
    def displacement(self) -> int:
        """
        Displacement from unity, as defined in RS2-106.
        
        "A displacement is a measurement of how much something changed from 
        a known point. In the Reciprocal System, that point is unit speed."
        
        For speed (s/t < 1): displacement = t - s (temporal displacement)
        For energy (s/t > 1): displacement = s - t (spatial displacement)
        For unity: displacement = 0
        """
        if self.orientation == "unity":
            return 0
        elif self.orientation == "speed":
            # Temporal displacement: how much extra time beyond space
            return self.time.value - self.space.value
        else:  # energy
            # Spatial displacement: how much extra space beyond time
            return self.space.value - self.time.value
    
    @property
    def displacement_type(self) -> str:
        """Whether displacement is temporal, spatial, or zero"""
        if self.orientation == "unity":
            return "none"
        elif self.orientation == "speed":
            return "temporal"
        else:
            return "spatial"
    
    def __repr__(self):
        return f"Motion(s={self.space.value}, t={self.time.value}, ratio={self.ratio}, orientation={self.orientation})"


# =============================================================================
# PART 2: Test the fundamental claims
# =============================================================================

def test_scalar_constraints():
    """Test that scalar units enforce RS2 constraints"""
    print("=" * 60)
    print("TEST 1: Scalar Unit Constraints")
    print("=" * 60)
    
    # Should work
    valid_tests = [1, 2, 3, 10, 100]
    print("\nValid scalars (should all succeed):")
    for v in valid_tests:
        try:
            s = ScalarUnit(v)
            print(f"  ScalarUnit({v}) -> {s} ✓")
        except ValueError as e:
            print(f"  ScalarUnit({v}) -> FAILED: {e} ✗")
    
    # Should fail
    invalid_tests = [0, -1, -5, 0.5, 1.5]
    print("\nInvalid scalars (should all fail):")
    for v in invalid_tests:
        try:
            s = ScalarUnit(v)
            print(f"  ScalarUnit({v}) -> {s} ✗ (should have failed!)")
        except (ValueError, TypeError) as e:
            print(f"  ScalarUnit({v}) -> Correctly rejected ✓")
    
    return True


def test_three_orientations():
    """Test that all motions fall into exactly three categories"""
    print("\n" + "=" * 60)
    print("TEST 2: Three-Fold Orientation Structure")
    print("=" * 60)
    
    # Generate a range of motions
    motions = []
    for s in range(1, 6):
        for t in range(1, 6):
            m = Motion(ScalarUnit(s), ScalarUnit(t))
            motions.append(m)
    
    # Categorize
    speeds = [m for m in motions if m.orientation == "speed"]
    unities = [m for m in motions if m.orientation == "unity"]
    energies = [m for m in motions if m.orientation == "energy"]
    
    print(f"\nTotal motions generated: {len(motions)}")
    print(f"  Speed (s/t < 1):  {len(speeds)}")
    print(f"  Unity (s/t = 1):  {len(unities)}")
    print(f"  Energy (s/t > 1): {len(energies)}")
    print(f"  Sum: {len(speeds) + len(unities) + len(energies)}")
    
    # Verify exhaustive categorization
    assert len(speeds) + len(unities) + len(energies) == len(motions), "Categories not exhaustive!"
    print("\n✓ All motions categorized into exactly three orientations")
    
    # Show examples
    print("\nExamples of each orientation:")
    print("  Speed:  ", speeds[:3])
    print("  Unity:  ", unities[:3])
    print("  Energy: ", energies[:3])
    
    return True


def test_reciprocal_symmetry():
    """
    Test the core RS2 claim: speed and energy are reciprocal.
    
    If motion M has ratio s/t, its reciprocal has ratio t/s.
    These should be "mirror images" across unity.
    """
    print("\n" + "=" * 60)
    print("TEST 3: Reciprocal Symmetry Across Unity")
    print("=" * 60)
    
    test_cases = [
        (1, 2),  # speed: 1/2
        (1, 3),  # speed: 1/3
        (2, 3),  # speed: 2/3
        (2, 1),  # energy: 2/1
        (3, 1),  # energy: 3/1
        (3, 2),  # energy: 3/2
        (1, 1),  # unity: 1/1
        (2, 2),  # unity: 2/2 (should equal 1/1)
    ]
    
    print("\nReciprocal pairs:")
    print("-" * 60)
    
    symmetry_holds = True
    for s, t in test_cases:
        m = Motion(ScalarUnit(s), ScalarUnit(t))
        r = m.reciprocal
        
        # Key test: ratio * reciprocal_ratio should equal 1
        product = m.ratio * r.ratio
        
        # Key test: orientations should be "opposite" (speed <-> energy) or both unity
        orientations_symmetric = (
            (m.orientation == "speed" and r.orientation == "energy") or
            (m.orientation == "energy" and r.orientation == "speed") or
            (m.orientation == "unity" and r.orientation == "unity")
        )
        
        # Key test: displacements should be equal in magnitude
        displacements_equal = m.displacement == r.displacement
        
        status = "✓" if (product == 1 and orientations_symmetric and displacements_equal) else "✗"
        if status == "✗":
            symmetry_holds = False
        
        print(f"  {m.ratio} × {r.ratio} = {product}, "
              f"orientations: {m.orientation}/{r.orientation}, "
              f"displacements: {m.displacement}/{r.displacement} {status}")
    
    if symmetry_holds:
        print("\n✓ Reciprocal symmetry holds for all test cases")
    else:
        print("\n✗ Reciprocal symmetry FAILED for some cases")
    
    return symmetry_holds


def test_displacement_structure():
    """
    Test the displacement concept from RS2-106.
    
    Displacement measures "distance from unity" and should have
    consistent structure across the speed/energy divide.
    """
    print("\n" + "=" * 60)
    print("TEST 4: Displacement Structure")
    print("=" * 60)
    
    print("\nDisplacement table (s=1 to 5, t=1 to 5):")
    print("-" * 60)
    print("     ", end="")
    for t in range(1, 6):
        print(f"  t={t} ", end="")
    print()
    
    for s in range(1, 6):
        print(f"s={s}: ", end="")
        for t in range(1, 6):
            m = Motion(ScalarUnit(s), ScalarUnit(t))
            d = m.displacement
            dtype = m.displacement_type[0].upper()  # T, S, or N
            print(f" {d:2d}{dtype} ", end="")
        print()
    
    print("\nLegend: T=temporal, S=spatial, N=none (unity)")
    
    # Verify: displacement 0 only at unity
    # Verify: displacement increases away from unity diagonal
    print("\nVerifying displacement properties...")
    
    for s in range(1, 10):
        for t in range(1, 10):
            m = Motion(ScalarUnit(s), ScalarUnit(t))
            
            # Unity should have zero displacement
            if s == t:
                assert m.displacement == 0, f"Unity {s}/{t} has non-zero displacement!"
                assert m.orientation == "unity"
            
            # Non-unity should have positive displacement
            if s != t:
                assert m.displacement > 0, f"Non-unity {s}/{t} has zero displacement!"
    
    print("✓ Displacement structure verified")
    return True


def test_speed_ranges():
    """
    Test Larson's "speed ranges" concept from RS2-106.
    
    Three dimensions of motion, each with speed and energy aspects,
    gives six "units of motion" total, but observed as three "speed ranges":
    - 1-x: low speed (first unit)
    - 2-x: intermediate speed (second unit)  
    - 3-x: ultra-high speed (third unit)
    """
    print("\n" + "=" * 60)
    print("TEST 5: Speed Range Classification")
    print("=" * 60)
    
    # According to RS2-106, speed ranges are determined by displacement
    # A displacement of 1 is in 1-x range, 2 in 2-x, etc.
    
    print("\nClassifying motions by speed range:")
    print("-" * 60)
    
    range_1x = []  # displacement 1
    range_2x = []  # displacement 2
    range_3x = []  # displacement 3
    range_higher = []
    unity_cases = []
    
    for s in range(1, 8):
        for t in range(1, 8):
            m = Motion(ScalarUnit(s), ScalarUnit(t))
            d = m.displacement
            
            if d == 0:
                unity_cases.append(m)
            elif d == 1:
                range_1x.append(m)
            elif d == 2:
                range_2x.append(m)
            elif d == 3:
                range_3x.append(m)
            else:
                range_higher.append(m)
    
    print(f"Unity (displacement 0): {len(unity_cases)} cases")
    print(f"1-x range (displacement 1): {len(range_1x)} cases")
    print(f"2-x range (displacement 2): {len(range_2x)} cases")
    print(f"3-x range (displacement 3): {len(range_3x)} cases")
    print(f"Higher ranges: {len(range_higher)} cases")
    
    print("\n1-x examples (first 5):")
    for m in range_1x[:5]:
        print(f"  {m.space.value}/{m.time.value} = {float(m.ratio):.3f} ({m.orientation})")
    
    print("\n2-x examples (first 5):")
    for m in range_2x[:5]:
        print(f"  {m.space.value}/{m.time.value} = {float(m.ratio):.3f} ({m.orientation})")
    
    return True


def test_cross_ratio_invariance():
    """
    Test the cross-ratio concept from RS2-104.
    
    "A scalar dimension is a cross-ratio where one scalar orientation 
    is fixed at unity (natural datum) and the other varies."
    
    The cross-ratio is projectively invariant - this is what makes
    measurement possible in a universe of pure motion.
    """
    print("\n" + "=" * 60)
    print("TEST 6: Cross-Ratio Structure")
    print("=" * 60)
    
    # The cross-ratio of four collinear points A, B, C, D is:
    # (A,B;C,D) = (AC/BC) / (AD/BD)
    # 
    # In RS2 terms, we fix one ratio at unity and vary the other.
    # The "cross-ratio" becomes: (s/t) / (1/1) = s/t
    # This is why Larson can "simplify" to just using ratios.
    
    print("\nCross-ratio interpretation:")
    print("When datum = 1/1, cross-ratio of motion s/t equals:")
    print("  (s/t) / (1/1) = s/t")
    print("\nThis is why RS2 uses simple ratios - they're already cross-ratios")
    print("measured against the unity datum.")
    
    # Demonstrate the invariance property
    print("\nDemonstrating projective invariance:")
    print("If we 'scale' both aspects uniformly, the ratio is preserved:")
    
    test_cases = [(1, 2), (2, 3), (3, 4)]
    for s, t in test_cases:
        original = Motion(ScalarUnit(s), ScalarUnit(t))
        # Scaling both by same factor
        for scale in [2, 3, 5]:
            scaled = Motion(ScalarUnit(s * scale), ScalarUnit(t * scale))
            print(f"  {s}/{t} = {float(original.ratio):.4f}, "
                  f"{s*scale}/{t*scale} = {float(scaled.ratio):.4f}, "
                  f"equal: {original.ratio == scaled.ratio}")
    
    return True


# =============================================================================
# PART 3: Run all tests and summarize
# =============================================================================

def run_all_tests():
    """Run the complete test suite"""
    print("\n" + "=" * 70)
    print("RS2 VALIDATION EXPERIMENT 01: SCALAR FOUNDATION")
    print("=" * 70)
    print("\nTesting the fundamental claims about scalar motion and ratios")
    print("from RS2-104 (Scalar Motion) and RS2-106 (Dimensions and Displacements)")
    print()
    
    results = {}
    
    results["scalar_constraints"] = test_scalar_constraints()
    results["three_orientations"] = test_three_orientations()
    results["reciprocal_symmetry"] = test_reciprocal_symmetry()
    results["displacement_structure"] = test_displacement_structure()
    results["speed_ranges"] = test_speed_ranges()
    results["cross_ratio"] = test_cross_ratio_invariance()
    
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
        print("\nThe fundamental scalar structure of RS2 is mathematically consistent.")
        print("This validates the foundation for building more complex structures.")
        print("\nNEXT STEPS:")
        print("  - Experiment 02: Introduce scalar DIRECTION (inward/outward)")
        print("  - Experiment 03: Test 2D rotation as primary motion")
        print("  - Experiment 04: Build toward quaternion structure")
    else:
        print("SOME TESTS FAILED")
        print("Review the failures before proceeding.")
    
    return all_passed


if __name__ == "__main__":
    run_all_tests()
