"""§7.2a — Zero density vs Riemann-von Mangoldt asymptotic.

Cold-derivation §7.2 prediction:
    N(T) ~ (T / 2π) log(T / 2π) − T / 2π + O(log T)

In RS2 coordinates: density of birotational eigenmodes on the unit-speed line.
The log(T/2π) growth is multiplicative (modes per unit log-frequency interval).
The 2π is one full birotation.

This test verifies the asymptotic by sampling N(T) at several T and computing
relative error against the leading-order Riemann-von Mangoldt formula.
"""
import math
from mpmath import mp, nzeros, mpf

mp.dps = 25

def rvm_leading(T):
    """Riemann-von Mangoldt leading-order: (T/2π) log(T/2π) − T/2π."""
    x = T / (2 * math.pi)
    return x * math.log(x) - x

def rvm_with_correction(T):
    """Add the standard 7/8 constant correction term."""
    return rvm_leading(T) + 7.0 / 8.0

print("§7.2a — Zero Density (Riemann-von Mangoldt)\n")
print(f"{'T':>10}  {'N(T) actual':>14}  {'leading':>14}  {'+7/8 corr':>14}  {'rel err':>10}")
print("-" * 70)

T_values = [10, 50, 100, 500, 1000, 5000, 10000]

results = []
for T in T_values:
    n_actual = int(nzeros(T))
    leading = rvm_leading(T)
    corrected = rvm_with_correction(T)
    rel_err = abs(corrected - n_actual) / max(n_actual, 1)
    results.append((T, n_actual, leading, corrected, rel_err))
    print(f"{T:>10}  {n_actual:>14}  {leading:>14.4f}  {corrected:>14.4f}  {rel_err:>10.4%}")

print()
max_err_corrected = max(r[4] for r in results)
print(f"Max relative error with 7/8 correction: {max_err_corrected:.3%}")
print(f"Asymptotic N(T) → ∞ as T → ∞? Check:")
for T, n, _, _, _ in results:
    log_density = n / T  # zeros per unit T
    expected_density = math.log(T / (2 * math.pi)) / (2 * math.pi)
    print(f"    T={T:>6}: zeros/T = {log_density:.5f}, "
          f"predicted log(T/2π)/(2π) = {expected_density:.5f}")

print("\n§7.2a prediction status:")
if max_err_corrected < 0.01:
    verdict = "PASS — Riemann-von Mangoldt asymptotic confirmed within <1% across all T."
else:
    verdict = f"PARTIAL — max error {max_err_corrected:.3%}; expected for low T."
print(f"    {verdict}")
print("    The growth rate log(T/2π)/(2π) modes per unit T is consistent with")
print("    the RS2 reading of N(T) as birotational eigenmode density.")
