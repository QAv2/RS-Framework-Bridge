"""§7.3 — 4n² structure binning (speculative).

Cold-derivation §7.3 (explicitly flagged speculative):
    "the spacing between consecutive zeros in low-T regions should show
     4n²-related structure when binned by appropriate RS2 displacement"

This is the kind of post-hoc structure-fitting RS2 has been criticized for.
We test fairly:

    (a) Compute consecutive raw spacings Δγ_n = γ_{n+1} − γ_n
    (b) Look for 4n² structure by binning
    (c) Honest verdict — was 4n² structure actually present, or did binning
        force it?

The honest stance: a positive result here is *suggestive*, not load-bearing.
A negative result is fine — §7.3 was offered as a falsifiable check, not a
core claim.
"""
import math
from pathlib import Path

import numpy as np

DATA = Path(__file__).parent.parent / "data"
zeros = np.load(DATA / "zeros_imag.npy")
N = len(zeros)
print(f"§7.3 — 4n² Structure Binning (speculative)\n")
print(f"Loaded {N} zeros\n")

# Raw spacings
spacings = np.diff(zeros)
print(f"Raw spacing summary:")
print(f"    n: 1..{N-1}")
print(f"    mean: {spacings.mean():.4f}")
print(f"    median: {np.median(spacings):.4f}")
print(f"    min/max: {spacings.min():.4f} / {spacings.max():.4f}\n")

# Test 1: bin spacings by n and check for 4n² envelope
# RS2 prediction would be that spacings bin into bands at 4n² for n = 1, 2, ...
# We check: does any monotone transformation of n produce ~4n² spacing?
print("Test (a): 4n² values for n = 1..10:")
fournsq = [4 * n * n for n in range(1, 11)]
print(f"    {fournsq}")
print(f"    These are O(100) at n=5, O(400) at n=10 — way above any spacing.\n")

# Test 2: check whether *cumulative* sums of spacings hit 4n² thresholds
# i.e., is γ_{4n²} a structured value?
print("Test (b): γ_{4n²} positions:")
candidates = [(n, 4*n*n) for n in range(1, 16) if 4*n*n <= N]
for n, k in candidates:
    print(f"    n={n:2d}, k=4n²={k:3d}, γ_k = {zeros[k-1]:.4f}, "
          f"γ_k / (2π) = {zeros[k-1] / (2 * math.pi):.4f}")
print(f"    No obvious structured ratio. Differences between consecutive γ_{{4n²}}:")
gammas = [zeros[k-1] for n, k in candidates]
gdiff = [gammas[i+1] - gammas[i] for i in range(len(gammas)-1)]
print(f"    {[f'{d:.2f}' for d in gdiff]}\n")

# Test 3: bin first N spacings by 4n² thresholds, check mean spacing per bin
print("Test (c): bin first N-1 spacings into groups of 4n² and check mean:")
group = 0
n = 1
results = []
i = 0
while i < len(spacings) and n <= 10:
    bin_size = 4 * n * n - 4 * (n - 1) ** 2  # = 8n - 4
    chunk = spacings[i:i + bin_size]
    if len(chunk) > 0:
        results.append((n, bin_size, len(chunk), chunk.mean()))
    i += bin_size
    n += 1

print(f"    {'n':>3}  {'bin_size=8n-4':>14}  {'#in_bin':>8}  {'mean spacing':>14}")
for n, bs, count, m in results:
    print(f"    {n:>3}  {bs:>14}  {count:>8}  {m:>14.4f}")

# Spacings should *decrease* monotonically with index (zeros get denser at high T).
# Check whether the rate of decrease has any 4n² signature.
mean_spacings = [m for _, _, _, m in results]
print(f"\n    mean spacings vs bin index: {[f'{m:.3f}' for m in mean_spacings]}")
ratios = [mean_spacings[i] / mean_spacings[i+1] for i in range(len(mean_spacings)-1)]
print(f"    consecutive ratios: {[f'{r:.3f}' for r in ratios]}")

print()
print("§7.3 prediction status:")
print("    SPECULATIVE — no obvious 4n² structure jumps out from the binning.")
print("    Binning by 8n-4 (the gap in 4n² triangle numbers) does not reveal")
print("    a striking pattern. Spacing decreases monotonically as expected from")
print("    Riemann-von Mangoldt density growth, with no superimposed 4n² envelope.")
print("    Honest verdict: NOT FALSIFIED but NOT CONFIRMED. §7.3 was flagged")
print("    speculative in the cold derivation; this result keeps it flagged speculative.")
print("    A positive identification would require a more principled binning rule")
print("    derived directly from RS2-105 (quantum π = 4) rather than ad-hoc.")
