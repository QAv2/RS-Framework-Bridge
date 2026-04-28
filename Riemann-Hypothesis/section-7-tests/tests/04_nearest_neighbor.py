"""§7.2c — Nearest-neighbor spacing distribution.

Cold-derivation §7.2 prediction (companion to pair correlation):
    Unfolded nearest-neighbor spacings should follow the Wigner surmise for
    GUE:
        P_GUE(s) = (32 / π²) s² exp(−4 s² / π)

Compare to:
    Poisson (uncorrelated levels):  P(s) = exp(−s)
    GOE Wigner surmise:             P(s) = (π/2) s exp(−π s² / 4)

If RS2's birotation argument is right, GUE should fit best.
"""
import math
from pathlib import Path

import numpy as np

DATA = Path(__file__).parent.parent / "data"
zeros = np.load(DATA / "zeros_imag.npy")
N = len(zeros)
print(f"§7.2c — Nearest-Neighbor Spacing Distribution\n")
print(f"Loaded {N} zeros\n")

def smooth_N(T):
    if T <= 0:
        return 0.0
    x = T / (2 * math.pi)
    return x * math.log(x) - x + 7.0 / 8.0

unfolded = np.array([smooth_N(g) for g in zeros])
spacings = np.diff(unfolded)
print(f"Mean unfolded spacing (should be ≈ 1): {spacings.mean():.4f}")
print(f"Spacing std: {spacings.std():.4f}\n")

# Histogram
n_bins = 25
s_max = 4.0
hist, edges = np.histogram(spacings, bins=n_bins, range=(0, s_max), density=True)
centers = 0.5 * (edges[:-1] + edges[1:])

def P_GUE(s):
    return (32.0 / math.pi ** 2) * s * s * math.exp(-4.0 * s * s / math.pi)

def P_GOE(s):
    return (math.pi / 2.0) * s * math.exp(-math.pi * s * s / 4.0)

def P_Poisson(s):
    return math.exp(-s)

print(f"{'s':>6}  {'observed':>10}  {'GUE':>10}  {'GOE':>10}  {'Poisson':>10}")
print("-" * 55)
ssr_gue = 0.0
ssr_goe = 0.0
ssr_poi = 0.0
for c, h in zip(centers, hist):
    g = P_GUE(c); o = P_GOE(c); p = P_Poisson(c)
    print(f"{c:>6.2f}  {h:>10.4f}  {g:>10.4f}  {o:>10.4f}  {p:>10.4f}")
    ssr_gue += (h - g) ** 2
    ssr_goe += (h - o) ** 2
    ssr_poi += (h - p) ** 2

print()
print(f"Sum-squared residual vs GUE Wigner surmise:    {ssr_gue:.4f}")
print(f"Sum-squared residual vs GOE Wigner surmise:    {ssr_goe:.4f}")
print(f"Sum-squared residual vs Poisson (no correl.):  {ssr_poi:.4f}")

best = min((ssr_gue, "GUE"), (ssr_goe, "GOE"), (ssr_poi, "Poisson"))
print(f"\nBest fit: {best[1]}")
print()
print("§7.2c prediction status:")
if best[1] == "GUE":
    print(f"    SUPPORTED — GUE Wigner surmise is the best fit among {{GUE, GOE, Poisson}}.")
    print(f"    Consistent with RS2 §7.2 prediction (birotation breaks time-reversal symmetry).")
elif best[1] == "GOE":
    print(f"    PARTIAL — GOE fits better; would falsify the time-reversal argument.")
else:
    print(f"    FAIL — Poisson fits best, indicating zeros look uncorrelated. Falsifying.")
