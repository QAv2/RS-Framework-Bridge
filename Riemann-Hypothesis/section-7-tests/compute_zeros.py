"""Compute imaginary parts of first N non-trivial Riemann zeros and cache to .npy.

Uses mpmath.zetazero (Riemann-Siegel under the hood). One-time cost; reused by
all section-7 tests. Saves incrementally so an interrupted run resumes from
the last checkpoint.
"""
import sys
import time
from pathlib import Path

import numpy as np
from mpmath import mp, zetazero

DATA = Path(__file__).parent / "data"
DATA.mkdir(exist_ok=True)
OUT = DATA / "zeros_imag.npy"
LOG = DATA / "zeros_imag.log"

N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
mp.dps = 25  # 25 decimal digits — overkill for stats, cheap insurance

if OUT.exists():
    existing = np.load(OUT)
    start = len(existing)
    print(f"resuming from n={start}")
    zeros = list(existing)
else:
    start = 0
    zeros = []

t0 = time.time()
last_save = t0
for n in range(start + 1, N + 1):
    z = zetazero(n)
    zeros.append(float(z.imag))
    now = time.time()
    if now - last_save > 30 or n == N:
        np.save(OUT, np.array(zeros))
        elapsed = now - t0
        rate = (n - start) / elapsed if elapsed > 0 else 0
        eta = (N - n) / rate if rate > 0 else float("inf")
        msg = f"n={n}/{N}  γ={zeros[-1]:.6f}  rate={rate:.2f}/s  eta={eta:.0f}s"
        print(msg)
        with open(LOG, "a") as f:
            f.write(msg + "\n")
        last_save = now

np.save(OUT, np.array(zeros))
print(f"done. saved {len(zeros)} zeros to {OUT}")
