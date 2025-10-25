#!/usr/bin/env python3
"""
permutation_tests.py
Reusable circular-shift permutation significance testing.

Usage (example):
    python analysis/scripts/permutation_tests.py --data data/PerBin_Averages_Christian_v2.csv
This will do a tiny self-test by correlating residuals with a dummy signal.
"""

import argparse
import numpy as np
import pandas as pd

def z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x))/s if s > 0 else np.zeros_like(x)

def perm_p(R: np.ndarray, D: np.ndarray, n: int = 2000, seed: int = 123) -> tuple[float, float]:
    """Circular-shift permutation p-value for zero-lag correlation."""
    R = np.asarray(R, float)
    D = np.asarray(D, float)
    r_obs = float(np.corrcoef(R, D)[0, 1])
    rng = np.random.default_rng(seed)
    npts = len(D)
    count = 0
    for _ in range(n):
        k = rng.integers(0, npts)
        rp = float(np.corrcoef(R, np.roll(D, k))[0, 1])
        if abs(rp) >= abs(r_obs):
            count += 1
    p = (count + 1) / (n + 1)
    return r_obs, p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/PerBin_Averages_Christian_v2.csv",
                    help="CSV with columns Start, End, Obs_avg, Exp_best_avg")
    ap.add_argument("--perms", type=int, default=500, help="number of permutations (self-test only)")
    args = ap.parse_args()

    perbin = pd.read_csv(args.data).sort_values("Start").copy()
    perbin["Residual"] = perbin["Obs_avg"] - perbin["Exp_best_avg"]
    R = z(perbin["Residual"].values)

    # Quick self-test: correlate residuals with its own 1-step shift (meaningless but shows API)
    D = np.roll(R, 1)
    r, p = perm_p(R, D, n=args.perms, seed=7)
    print(f"[self-test] r={r:.3f}, p≈{p:.3f} (no scientific meaning—demo only)")

if __name__ == "__main__":
    main()
