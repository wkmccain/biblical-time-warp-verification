#!/usr/bin/env python3
"""
monte_carlo_10k.py
------------------
Generate per-bin observed vs expected averages from 10,000 Monte Carlo draws,
following the preregistered procedure described in the package docs:
- Sample each event's date uniformly within its [start, end] range (per draw).
- Bin sampled mid-years into 25-year bins from MIN_YEAR to MAX_YEAR.
- Fit secular baselines on TRAIN (<= 0 CE) only:
    * constant rate (global mean across training bins)
    * piecewise-constant with 8 segments (equal-size segments over training bins)
    * polynomial (degree=5) on training bins (spline-like, fixed complexity)
- Select the baseline with smallest Poisson deviance on TEST (0..100 CE).
- Average observed counts and best-baseline expectations across draws.
- Write CSV with columns: Start, End, Obs_avg, Exp_best_avg, Bin.

Usage:
    python monte_carlo_10k.py \
        --events context/EventList_Biblical_Encounters.csv \
        --out data/PerBin_Averages_Christian_v2.csv \
        --draws 10000 --seed 2025

Notes:
- This script assumes the event list provides clustered anchors with date ranges
  in columns 'start' and 'end'. All events are weighted equally (1.0). If you
  have a 'weight' column, it will be used automatically.
- Binning defaults to 25-year bins spanning -2000..125 (inclusive of the end
  boundary on the last bin). Adjust with --min_year/--max_year/--bin_size.
- Dependencies: numpy, pandas.
"""

import argparse, os
import numpy as np, pandas as pd

def build_bins(min_year=-2000, max_year=125, bin_size=25):
    edges = np.arange(min_year, max_year + bin_size, bin_size, dtype=float)
    starts = edges[:-1].astype(int)
    ends   = edges[1:].astype(int)
    mids   = (edges[:-1] + edges[1:]) / 2.0
    labels = [f"{s}–{e-1}" for s, e in zip(starts, ends)]
    return starts, ends, mids, edges, labels

def poisson_deviance(y, mu, eps=1e-12):
    mu = np.clip(mu, eps, None)
    y = np.asarray(y, float)
    term = np.zeros_like(y, float)
    nz = y > 0
    term[nz] = y[nz] * np.log(y[nz] / mu[nz]) - (y[nz] - mu[nz])
    term[~nz] = -mu[~nz]
    return 2.0 * np.sum(term)

def fit_constant(train_y):
    lam = float(np.mean(train_y))
    def predict(n):
        return np.full(n, lam, dtype=float)
    return predict

def fit_piecewise8(train_y, n_train_bins):
    # Split the training bins into 8 equal-size contiguous segments.
    seg_sizes = [n_train_bins // 8] * 8
    for i in range(n_train_bins % 8):
        seg_sizes[i] += 1
    seg_bounds = np.cumsum([0] + seg_sizes)  # length 9
    seg_means = []
    for i in range(8):
        lo, hi = seg_bounds[i], seg_bounds[i+1]
        seg_means.append(float(np.mean(train_y[lo:hi])) if hi > lo else 0.0)
    seg_means = np.array(seg_means, float)

    def predict(n_total_bins):
        # For bins beyond the training region, use the last segment's mean.
        out = np.empty(n_total_bins, float)
        for i in range(n_total_bins):
            if i < n_train_bins:
                # Map i to segment index
                seg_idx = np.searchsorted(seg_bounds[1:], i, side='right')
                out[i] = seg_means[seg_idx]
            else:
                out[i] = seg_means[-1]
        return out
    return predict

def fit_poly5(train_x, train_y):
    # Fit degree-5 polynomial on training bin midpoints vs counts
    coeff = np.polyfit(train_x, train_y, 5)
    p = np.poly1d(coeff)
    def predict(x_all):
        yhat = p(x_all)
        return np.clip(yhat, 1e-12, None)
    return predict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="context/EventList_Biblical_Encounters.csv",
                    help="CSV with columns at least: start, end; optional 'weight'")
    ap.add_argument("--out", default="data/PerBin_Averages_Christian_v2.csv",
                    help="Output CSV path (will be created)")
    ap.add_argument("--draws", type=int, default=10000, help="Number of Monte Carlo draws")
    ap.add_argument("--seed", type=int, default=2025, help="RNG seed")
    ap.add_argument("--min_year", type=int, default=-2000)
    ap.add_argument("--max_year", type=int, default=125)
    ap.add_argument("--bin_size", type=int, default=25)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Load clustered events with date ranges
    ev = pd.read_csv(args.events)
    if not {"start", "end"}.issubset(ev.columns):
        raise ValueError("Event CSV must include 'start' and 'end' columns (year bounds).")
    weights = ev["weight"].values.astype(float) if "weight" in ev.columns else np.ones(len(ev), float)
    starts, ends, mids, edges, labels = build_bins(args.min_year, args.max_year, args.bin_size)
    n_bins = len(starts)

    # Identify training and test masks (by bin End year)
    bin_ends = ends  # integer
    train_mask = bin_ends <= 0
    test_mask  = (bin_ends > 0) & (bin_ends <= 100)
    n_train = int(np.sum(train_mask))

    rng = np.random.default_rng(args.seed)

    obs_accum = np.zeros(n_bins, float)
    exp_accum = np.zeros(n_bins, float)
    picks = {"constant": 0, "piece8": 0, "poly5": 0}

    # Precompute bin indices helper
    min_year = args.min_year
    bs = float(args.bin_size)

    # Perform draws
    for d in range(args.draws):
        # Sample a mid-year for each event from Uniform[start, end]
        y = rng.uniform(ev["start"].values, ev["end"].values, size=len(ev))
        # Compute bin index for each sampled year
        idx = np.floor((y - min_year) / bs).astype(int)
        # Clip to valid range
        idx = np.clip(idx, 0, n_bins - 1)
        # Accumulate observed counts per bin (weighted)
        obs = np.zeros(n_bins, float)
        np.add.at(obs, idx, weights)

        # Fit baselines on TRAIN
        train_y = obs[train_mask]
        train_x = mids[train_mask]

        # Constant
        const_pred = fit_constant(train_y)(n_bins)

        # Piecewise-8
        piece_pred = fit_piecewise8(train_y, n_train)(n_bins)

        # Poly (degree 5)
        poly_pred = fit_poly5(train_x, train_y)(mids)

        # Evaluate on TEST
        y_test = obs[test_mask]
        const_dev = poisson_deviance(y_test, const_pred[test_mask])
        piece_dev = poisson_deviance(y_test, piece_pred[test_mask])
        poly_dev  = poisson_deviance(y_test,  poly_pred[test_mask])

        devs = np.array([const_dev, piece_dev, poly_dev])
        choice = int(np.argmin(devs))
        if choice == 0:
            exp = const_pred; picks["constant"] += 1
        elif choice == 1:
            exp = piece_pred; picks["piece8"] += 1
        else:
            exp = poly_pred; picks["poly5"] += 1

        # Accumulate
        obs_accum += obs
        exp_accum += exp

    # Averages across draws
    obs_avg = obs_accum / args.draws
    exp_avg = exp_accum / args.draws

    out = pd.DataFrame({
        "Start": starts,
        "End": ends,
        "Obs_avg": np.round(obs_avg, 4),
        "Exp_best_avg": np.round(exp_avg, 4),
        "Bin": labels,
    })
    out.to_csv(args.out, index=False)

    # Write a small summary next to it
    summary_path = os.path.splitext(args.out)[0] + "_MC_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"# Monte Carlo summary ({args.draws} draws)\n")
        f.write(f"Events: {len(ev)}\n")
        f.write(f"Bins: {n_bins} ({args.min_year}..{args.max_year} step {args.bin_size})\n")
        f.write("Baseline picks:\n")
        for k,v in picks.items():
            f.write(f"  - {k}: {v}\n")

    print(f"Wrote {args.out}")
    print(f"Baseline picks: {picks}")

if __name__ == "__main__":
    main()
