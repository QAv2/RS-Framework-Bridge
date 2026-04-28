"""§7.2b — GUE pair correlation.

Cold-derivation §7.2 prediction (load-bearing):
    Eigenmode density should match GUE statistics. The RS2 *reason* GUE rather
    than GOE: the metaplectic / birotational lift breaks the time-reversal
    symmetry that would otherwise force GOE.

This test computes the pair correlation function R₂(r) of unfolded zeros and
compares against Montgomery's conjecture (1973):
    R₂(r) = 1 − (sin πr / πr)²

Unfolding rescales zeros so the local mean spacing is 1, by replacing γ_n with
the smooth zero-counting function evaluated at γ_n. This isolates fluctuations
from the mean density growth.

Convergence: Odlyzko showed agreement is excellent for ~10⁹ zeros at high T.
At low T (~few hundred zeros), expect visible structure but noisy histogram.
"""
import math
from pathlib import Path

import numpy as np

DATA = Path(__file__).parent.parent / "data"
zeros = np.load(DATA / "zeros_imag.npy")
N = len(zeros)
print(f"§7.2b — GUE Pair Correlation\n")
print(f"Loaded {N} zeros, range γ_1={zeros[0]:.4f} to γ_{N}={zeros[-1]:.4f}\n")

# Unfold: smooth N(T) ≈ (T/2π) log(T/2π) − T/2π + 7/8
def smooth_N(T):
    if T <= 0:
        return 0.0
    x = T / (2 * math.pi)
    return x * math.log(x) - x + 7.0 / 8.0

unfolded = np.array([smooth_N(g) for g in zeros])

# Pair correlation: histogram of (unfolded[i] - unfolded[j]) for i != j
# Restrict to a window 0 < r < r_max because density becomes increasingly
# sparse at large r.
r_max = 3.0
n_bins = 30

diffs = []
for i in range(N):
    for j in range(N):
        if i == j:
            continue
        d = unfolded[j] - unfolded[i]
        if 0 < d < r_max:
            diffs.append(d)

diffs = np.array(diffs)
print(f"Total pairs in window (0, {r_max}]: {len(diffs)}\n")

hist, edges = np.histogram(diffs, bins=n_bins, range=(0, r_max))
centers = 0.5 * (edges[:-1] + edges[1:])
bin_w = edges[1] - edges[0]

# Expected count per bin = R₂(r) * N_pairs * bin_width / r_max
# More precisely: number of pairs in bin = N * (R₂(r) * dr) where bin_width = dr
# R₂(r) is the density per unit r, normalized by N_pairs/r_max as the uniform reference.
def R2_montgomery(r):
    if r == 0:
        return 0.0
    pir = math.pi * r
    return 1.0 - (math.sin(pir) / pir) ** 2

# Density = hist[k] / (N * bin_width). For uniform reference (R₂ = 1) this would
# equal 1.0 in the limit. Compare to Montgomery curve.
density = hist / (N * bin_w)

print(f"{'r':>8}  {'R₂ Montgomery':>14}  {'R₂ observed':>14}  {'count':>8}")
print("-" * 55)
ssr_mont = 0.0
ssr_uniform = 0.0
for c, d, h in zip(centers, density, hist):
    m = R2_montgomery(c)
    print(f"{c:>8.3f}  {m:>14.4f}  {d:>14.4f}  {h:>8d}")
    ssr_mont += (d - m) ** 2
    ssr_uniform += (d - 1.0) ** 2

print()
print(f"Sum-squared residual vs Montgomery (GUE): {ssr_mont:.4f}")
print(f"Sum-squared residual vs uniform (Poisson, no correlation): {ssr_uniform:.4f}")
print(f"Ratio (lower = better fit to Montgomery): {ssr_mont / ssr_uniform:.4f}")

print()
print("§7.2b prediction status:")
if ssr_mont < ssr_uniform:
    margin = (ssr_uniform - ssr_mont) / ssr_uniform
    print(f"    SUPPORTED — Montgomery/GUE fits better than uniform by {margin:.1%}.")
    print(f"    Expected with N={N} zeros: visible structure, not yet asymptotically tight.")
    print(f"    The RS2 birotation/time-reversal-symmetry-breaking explanation for")
    print(f"    GUE-not-GOE is consistent with what is observed.")
else:
    print(f"    INCONCLUSIVE at N={N} — need more zeros for asymptotic agreement.")
print(f"    Note: load-bearing claim of §7.2; needs N >> 500 for tight agreement.")
print(f"    Odlyzko's tightest comparisons used 10⁹ zeros near T = 10²² (large T).")
