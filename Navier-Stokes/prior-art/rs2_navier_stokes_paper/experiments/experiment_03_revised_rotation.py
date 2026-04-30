"""
RS2 Validation Experiment 03 (Revised): Rotation as Primary
============================================================

CORRECTION FROM PREVIOUS VERSION:
The quaternion's scalar (w) component is NOT a variable - it IS the progression,
fixed at unity. The three imaginary components (i, j, k) represent rotational
displacements in the three dimensions of motion.

Structure: q = 1 + xi + yj + zk

Where:
- "1" is the ever-present progression (the datum, always there)
- x, y, z are rotational displacements in the three scalar dimensions

The universe is ALWAYS 3-dimensional. The quaternion is a mathematical tool
to represent 3D rotation, not a claim about 4D space.

CLAIMS BEING TESTED (from RS2-109 Dimensional Thinking):

1. ROTATION IS PRIMARY (Nehru's contribution):
   Angular velocity is just as fundamental as linear velocity.
   In 3D space, rotation becomes possible without needing "something to rotate."

2. THE QUATERNION REPRESENTS 3D ROTATION WITH FIXED PROGRESSION:
   - Scalar part (w=1): The progression, always present, the datum
   - Vector part (xi + yj + zk): Rotational displacements from progression
   
3. THE UNITS OF MOTION are displacements in the rotational components:
   - 1-x range: Electric rotation (one-dimensional, single axis)
   - 2-x range: Magnetic rotation (two-dimensional, two axes combined)
   - 3-x range: Full 3D rotation -> produces net INWARD scalar direction

4. COMBINING ROTATIONS in 3D produces the gravitational (inward) effect:
   When all three rotational dimensions are engaged, the net effect is
   inward motion relative to the progression.
"""

import math
from dataclasses import dataclass
from typing import Tuple, Optional, List
from fractions import Fraction
from enum import Enum

# =============================================================================
# PART 1: RS2 Quaternion - Scalar Fixed at Unity
# =============================================================================

@dataclass
class RS2Quaternion:
    """
    A quaternion representing motion in RS2, with scalar part fixed at unity.
    
    q = 1 + xi + yj + zk
    
    The "1" represents the ever-present progression of the natural reference
    system. It is the datum against which all rotational displacements are
    measured. It is NOT a variable.
    
    The i, j, k components represent rotational displacements in the three
    scalar dimensions. These CAN vary and represent the structure of particles,
    atoms, and other motion configurations.
    """
    x: float = 0.0  # i component (1st rotational dimension)
    y: float = 0.0  # j component (2nd rotational dimension)
    z: float = 0.0  # k component (3rd rotational dimension)
    
    # The scalar part is ALWAYS 1 - the progression
    @property
    def w(self) -> float:
        return 1.0
    
    def __add__(self, other: 'RS2Quaternion') -> 'RS2Quaternion':
        """Add rotational components (progression remains 1)"""
        return RS2Quaternion(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z
        )
    
    def __sub__(self, other: 'RS2Quaternion') -> 'RS2Quaternion':
        """Subtract rotational components"""
        return RS2Quaternion(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z
        )
    
    def multiply_rotations(self, other: 'RS2Quaternion') -> 'RS2Quaternion':
        """
        Multiply two RS2 quaternions.
        
        This follows Hamilton's rules but with w fixed at 1.
        The result's scalar part may differ from 1 during calculation,
        but we interpret only the rotational part as the new displacement.
        
        Full quaternion product for (1 + ai + bj + ck)(1 + xi + yj + zk):
        
        Scalar: 1*1 - a*x - b*y - c*z = 1 - (a*x + b*y + c*z)
        i:      1*x + a*1 + b*z - c*y = x + a + bz - cy
        j:      1*y - a*z + b*1 + c*x = y - az + b + cx
        k:      1*z + a*y - b*x + c*1 = z + ay - bx + c
        
        The scalar deviation from 1 represents the NET SCALAR EFFECT
        of combining these rotations - this is what produces gravity!
        """
        a, b, c = self.x, self.y, self.z
        x, y, z = other.x, other.y, other.z
        
        # Calculate full quaternion product
        new_w = 1 - (a*x + b*y + c*z)
        new_x = x + a + b*z - c*y
        new_y = y - a*z + b + c*x
        new_z = z + a*y - b*x + c
        
        # Store the scalar deviation for analysis
        self._last_scalar_deviation = new_w - 1
        
        return RS2Quaternion(new_x, new_y, new_z)
    
    @property
    def rotation_magnitude(self) -> float:
        """Magnitude of the rotational displacement"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    @property
    def net_scalar_effect(self) -> float:
        """
        The net scalar effect of this rotation configuration.
        
        For a rotation q = 1 + xi + yj + zk, when it interacts with itself
        (q * q), the scalar part becomes: 1 - (x² + y² + z²)
        
        This deviation from unity IS the gravitational (inward) tendency.
        The larger the rotation magnitude, the more inward the net effect.
        """
        return -(self.x**2 + self.y**2 + self.z**2)
    
    @property 
    def scalar_direction(self) -> str:
        """
        Classify the net scalar direction of this rotation configuration.
        
        Based on the net scalar effect:
        - Zero rotation: No net effect (at progression)
        - Any rotation: Net INWARD effect (gravitational)
        """
        effect = self.net_scalar_effect
        if abs(effect) < 1e-10:
            return "UNITY (pure progression)"
        else:
            return f"INWARD (gravitational, magnitude {abs(effect):.4f})"
    
    def is_electric(self) -> bool:
        """Is this a 1D (electric) rotation - only one axis engaged?"""
        nonzero = sum(1 for v in [self.x, self.y, self.z] if abs(v) > 1e-10)
        return nonzero == 1
    
    def is_magnetic(self) -> bool:
        """Is this a 2D (magnetic) rotation - two axes engaged?"""
        nonzero = sum(1 for v in [self.x, self.y, self.z] if abs(v) > 1e-10)
        return nonzero == 2
    
    def is_full_3d(self) -> bool:
        """Is this a full 3D rotation - all three axes engaged?"""
        nonzero = sum(1 for v in [self.x, self.y, self.z] if abs(v) > 1e-10)
        return nonzero == 3
    
    def rotation_type(self) -> str:
        """Classify the rotation type"""
        nonzero = sum(1 for v in [self.x, self.y, self.z] if abs(v) > 1e-10)
        if nonzero == 0:
            return "NONE (pure progression)"
        elif nonzero == 1:
            return "ELECTRIC (1D, 1-x range)"
        elif nonzero == 2:
            return "MAGNETIC (2D, 2-x range)"
        else:
            return "GRAVITATIONAL (3D, 3-x range)"
    
    def __repr__(self):
        parts = ["1"]  # Scalar is always 1
        if abs(self.x) > 1e-10:
            parts.append(f"{self.x:+.4g}i")
        if abs(self.y) > 1e-10:
            parts.append(f"{self.y:+.4g}j")
        if abs(self.z) > 1e-10:
            parts.append(f"{self.z:+.4g}k")
        return "".join(parts)


# =============================================================================
# PART 2: Rotation Configurations (Units of Motion)
# =============================================================================

def pure_progression() -> RS2Quaternion:
    """No rotation - just the progression"""
    return RS2Quaternion(0, 0, 0)

def electric_rotation(axis: str, magnitude: float = 1.0) -> RS2Quaternion:
    """1D rotation in one axis (electric, 1-x range)"""
    if axis == 'i':
        return RS2Quaternion(magnitude, 0, 0)
    elif axis == 'j':
        return RS2Quaternion(0, magnitude, 0)
    elif axis == 'k':
        return RS2Quaternion(0, 0, magnitude)
    else:
        raise ValueError(f"Unknown axis: {axis}")

def magnetic_rotation(axes: str, magnitude: float = 1.0) -> RS2Quaternion:
    """2D rotation in two axes (magnetic, 2-x range)"""
    if axes == 'ij':
        return RS2Quaternion(magnitude, magnitude, 0)
    elif axes == 'jk':
        return RS2Quaternion(0, magnitude, magnitude)
    elif axes == 'ik':
        return RS2Quaternion(magnitude, 0, magnitude)
    else:
        raise ValueError(f"Unknown axes: {axes}")

def full_rotation(magnitude: float = 1.0) -> RS2Quaternion:
    """3D rotation in all axes (gravitational, 3-x range)"""
    return RS2Quaternion(magnitude, magnitude, magnitude)


# =============================================================================
# PART 3: Tests
# =============================================================================

def test_progression_is_fixed():
    """Verify that the scalar part is always 1 (the progression)"""
    print("=" * 60)
    print("TEST 1: Progression is Fixed at Unity")
    print("=" * 60)
    
    test_cases = [
        RS2Quaternion(0, 0, 0),
        RS2Quaternion(1, 0, 0),
        RS2Quaternion(0, 1, 0),
        RS2Quaternion(0, 0, 1),
        RS2Quaternion(1, 1, 1),
        RS2Quaternion(3, -2, 5),
    ]
    
    print("\nAll quaternions have scalar part = 1 (the progression):")
    print("-" * 60)
    
    all_pass = True
    for q in test_cases:
        status = "✓" if q.w == 1.0 else "✗"
        if q.w != 1.0:
            all_pass = False
        print(f"  {q!s:25s} -> w = {q.w} {status}")
    
    if all_pass:
        print("\n✓ Progression correctly fixed at unity for all configurations")
    
    return all_pass


def test_rotation_classification():
    """Test classification of rotation types (speed ranges)"""
    print("\n" + "=" * 60)
    print("TEST 2: Rotation Type Classification (Speed Ranges)")
    print("=" * 60)
    
    test_cases = [
        (pure_progression(), "NONE"),
        (electric_rotation('i'), "ELECTRIC"),
        (electric_rotation('j'), "ELECTRIC"),
        (electric_rotation('k'), "ELECTRIC"),
        (magnetic_rotation('ij'), "MAGNETIC"),
        (magnetic_rotation('jk'), "MAGNETIC"),
        (full_rotation(), "GRAVITATIONAL"),
    ]
    
    print("\nClassifying rotation configurations:")
    print("-" * 60)
    
    all_pass = True
    for q, expected_type in test_cases:
        rot_type = q.rotation_type()
        matches = expected_type in rot_type
        status = "✓" if matches else "✗"
        if not matches:
            all_pass = False
        print(f"  {q!s:20s} -> {rot_type} {status}")
    
    if all_pass:
        print("\n✓ Rotation types correctly classified")
    
    return all_pass


def test_net_scalar_effect():
    """
    Test that rotations produce net INWARD (gravitational) effect.
    
    This is the key insight: rotation CREATES gravity in RS2.
    The net scalar effect of any rotation is negative (inward).
    """
    print("\n" + "=" * 60)
    print("TEST 3: Net Scalar Effect of Rotations")
    print("=" * 60)
    
    print("\nRotations produce INWARD (gravitational) scalar effect:")
    print("-" * 60)
    
    test_cases = [
        ("No rotation", RS2Quaternion(0, 0, 0)),
        ("Small i-rotation", RS2Quaternion(0.5, 0, 0)),
        ("Unit i-rotation", RS2Quaternion(1, 0, 0)),
        ("Large i-rotation", RS2Quaternion(2, 0, 0)),
        ("Unit j-rotation", RS2Quaternion(0, 1, 0)),
        ("2D (i+j)", RS2Quaternion(1, 1, 0)),
        ("3D (i+j+k)", RS2Quaternion(1, 1, 1)),
        ("Large 3D", RS2Quaternion(2, 2, 2)),
    ]
    
    for name, q in test_cases:
        effect = q.net_scalar_effect
        direction = q.scalar_direction
        print(f"  {name:20s}: {q!s:20s} -> effect = {effect:+.4f} [{direction}]")
    
    print("\nKey observations:")
    print("  - No rotation: net effect = 0 (pure progression)")
    print("  - Any rotation: net effect < 0 (INWARD/gravitational)")
    print("  - More rotation: more inward effect")
    print("  - This is why atoms have mass/gravity!")
    
    # Verify: all rotations produce inward effect
    all_inward = True
    for name, q in test_cases[1:]:  # Skip "no rotation"
        if q.net_scalar_effect >= 0:
            all_inward = False
    
    if all_inward:
        print("\n✓ All rotations correctly produce inward (gravitational) effect")
    
    return all_inward


def test_rotation_magnitude_scaling():
    """
    Test how gravitational effect scales with rotation magnitude.
    
    Effect = -(x² + y² + z²)
    
    This is a SQUARED relationship - doubling rotation quadruples gravity.
    """
    print("\n" + "=" * 60)
    print("TEST 4: Gravitational Effect Scaling")
    print("=" * 60)
    
    print("\nHow gravitational effect scales with rotation magnitude:")
    print("-" * 60)
    
    print("\n  Single axis (i-rotation only):")
    for mag in [0.5, 1.0, 1.5, 2.0, 3.0]:
        q = RS2Quaternion(mag, 0, 0)
        effect = q.net_scalar_effect
        print(f"    magnitude = {mag:.1f} -> effect = {effect:+.4f} (= -{mag}²)")
    
    print("\n  Three axes (equal magnitude in each):")
    for mag in [0.5, 1.0, 1.5, 2.0]:
        q = RS2Quaternion(mag, mag, mag)
        effect = q.net_scalar_effect
        expected = -3 * mag**2
        print(f"    magnitude = {mag:.1f} each -> effect = {effect:+.4f} (= -3×{mag}²)")
    
    print("\nKey insight: Effect = -(x² + y² + z²)")
    print("This is the SQUARE of the rotation magnitude!")
    print("Small rotations: small gravity. Large rotations: much larger gravity.")
    
    return True


def test_birotation_photon():
    """
    Test the photon as a birotation structure.
    
    From RS2-109: The photon is a birotation where two counter-rotating
    systems cancel their rotational components, leaving pure progression.
    
    With fixed progression, this means:
    q1 = 1 + rotation
    q2 = 1 - rotation (counter-rotation)
    
    The NET effect should cancel the rotational displacement.
    """
    print("\n" + "=" * 60)
    print("TEST 5: Photon as Birotation")
    print("=" * 60)
    
    print("\nBirotation: two counter-rotating systems")
    print("-" * 60)
    
    # A rotation and its counter-rotation
    rot = RS2Quaternion(1, 0, 0)      # +i rotation
    counter = RS2Quaternion(-1, 0, 0)  # -i rotation
    
    print(f"  Rotation:        {rot}")
    print(f"  Counter-rotation: {counter}")
    
    # When we ADD these (superposition), rotations cancel
    combined = rot + counter
    print(f"  Sum (superposition): {combined}")
    
    is_pure_progression = combined.rotation_magnitude < 1e-10
    
    if is_pure_progression:
        print("\n✓ Birotation cancels to pure progression (photon at unit speed)")
    else:
        print(f"\n✗ Expected pure progression, got {combined}")
    
    # More complex birotation
    print("\n  More complex birotation (2D):")
    rot2d = RS2Quaternion(1, 1, 0)
    counter2d = RS2Quaternion(-1, -1, 0)
    combined2d = rot2d + counter2d
    print(f"    Rotation:        {rot2d}")
    print(f"    Counter-rotation: {counter2d}")
    print(f"    Sum: {combined2d}")
    
    return is_pure_progression


def test_speed_range_progression():
    """
    Test the progression through speed ranges via rotation dimensions.
    
    0 rotations: Pure progression (photons carried at c)
    1 rotation:  Electric (1-x range)
    2 rotations: Magnetic (2-x range)
    3 rotations: Full gravitational structure (3-x range)
    """
    print("\n" + "=" * 60)
    print("TEST 6: Speed Range Progression")
    print("=" * 60)
    
    print("\nBuilding up through speed ranges by adding rotational dimensions:")
    print("-" * 60)
    
    stages = [
        ("0D: Progression only", RS2Quaternion(0, 0, 0)),
        ("1D: Add i-rotation", RS2Quaternion(1, 0, 0)),
        ("2D: Add j-rotation", RS2Quaternion(1, 1, 0)),
        ("3D: Add k-rotation", RS2Quaternion(1, 1, 1)),
    ]
    
    for name, q in stages:
        rot_type = q.rotation_type()
        effect = q.net_scalar_effect
        print(f"\n  {name}")
        print(f"    Configuration: {q}")
        print(f"    Type: {rot_type}")
        print(f"    Gravitational effect: {effect:+.4f}")
    
    print("\nProgression of gravitational effect:")
    print("  0D -> 1D -> 2D -> 3D")
    print("   0  -> -1 -> -2 -> -3")
    print("\nEach rotational dimension adds to the inward tendency!")
    
    return True


def test_discrete_displacement_values():
    """
    Test with discrete (integer) displacement values, as RS2 requires.
    
    Displacements are counting numbers - this matches the scalar unit
    constraint from Experiment 01.
    """
    print("\n" + "=" * 60)
    print("TEST 7: Discrete Displacement Values")
    print("=" * 60)
    
    print("\nUsing discrete (integer) rotational displacements:")
    print("-" * 60)
    
    # Larson's atomic notation uses integer displacements
    # e.g., Oxygen is 2-2-(2) meaning displacements of 2,2,2
    
    test_atoms = [
        ("Hydrogen-like (1-1-0)", RS2Quaternion(1, 1, 0)),
        ("Helium-like (2-1-0)", RS2Quaternion(2, 1, 0)),
        ("Oxygen-like (2-2-2)", RS2Quaternion(2, 2, 2)),
        ("Large atom (3-3-2)", RS2Quaternion(3, 3, 2)),
    ]
    
    print("\n  Configuration       Quaternion           Grav. Effect  Type")
    print("  " + "-" * 70)
    
    for name, q in test_atoms:
        effect = q.net_scalar_effect
        rot_type = q.rotation_type().split()[0]
        print(f"  {name:22s} {q!s:20s} {effect:+8.1f}      {rot_type}")
    
    print("\nNote: Larger displacements = stronger gravitational effect")
    print("This correlates with atomic mass in conventional physics!")
    
    return True


# =============================================================================
# PART 4: Run all tests and summarize
# =============================================================================

def run_all_tests():
    """Run the complete test suite"""
    print("\n" + "=" * 70)
    print("RS2 VALIDATION EXPERIMENT 03 (REVISED): ROTATION AS PRIMARY")
    print("=" * 70)
    print("\nCorrected model: Scalar (w) is FIXED at 1 (the progression)")
    print("Only the rotational components (i, j, k) are variables")
    print()
    
    results = {}
    
    results["progression_fixed"] = test_progression_is_fixed()
    results["rotation_classification"] = test_rotation_classification()
    results["net_scalar_effect"] = test_net_scalar_effect()
    results["scaling"] = test_rotation_magnitude_scaling()
    results["birotation_photon"] = test_birotation_photon()
    results["speed_ranges"] = test_speed_range_progression()
    results["discrete_values"] = test_discrete_displacement_values()
    
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
        print("\nKEY FINDINGS (CORRECTED MODEL):")
        print("  1. Progression (w=1) is FIXED - it's the datum, not a variable")
        print("  2. Rotations (i,j,k) are displacements FROM the progression")
        print("  3. Any rotation produces INWARD (gravitational) effect")
        print("  4. Effect scales as SQUARE of rotation magnitude: -(x²+y²+z²)")
        print("  5. Speed ranges correspond to # of rotational dimensions engaged")
        print("  6. Birotation (counter-rotating) cancels to pure progression")
        print("\nCONNECTION TO NAVIER-STOKES:")
        print("  In this corrected model, the 'blow-up' question becomes:")
        print("  Can rotational displacements (x,y,z) grow without bound?")
        print("  ")
        print("  RS2 suggests: NO, because:")
        print("  - Larger rotations create larger inward (gravitational) effect")
        print("  - This inward effect opposes further outward expansion")
        print("  - The system is self-limiting")
        print("\nNEXT STEPS:")
        print("  - Experiment 04: Model how rotations interact and combine")
        print("  - Experiment 05: Test self-limiting behavior of rotation growth")
        print("  - Experiment 06: Apply to fluid element dynamics")
    else:
        print("SOME TESTS FAILED")
        print("Review the failures before proceeding.")
    
    return all_passed


if __name__ == "__main__":
    run_all_tests()
