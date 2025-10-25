#!/usr/bin/env python3
"""
warp_grid_search.py
Grid search for the approach/contact/aftershock warp model (no Daniel).

Writes:
  data/AC_WARP_grid.csv   (rows of parameters with r0 = zero-lag correlation)
"""

import argparse
import numpy as np
import pandas as pd

def z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x))/s if s > 0 else np.zeros_like(x)

def speed_w(t, t_star=33.0, a=1.0, k=0.006, B=1.6, p=0.8, c=10.0):
    """Logistic approach + Omori-like aftershock tail."""
    t = np.asarray(t, float)
    logistic = 1.0 / (1.0 + np.exp(-k * (t - t_star)))
    w = 1.0 + a * logistic
    mask = t > t_star
    w = w.copy()
    w[mask] += B / np.power(c + (t[mask] - t_star), p)
    return w

def tau_from_w(t, w, dt):
    tau = np.cumsum(w * dt)
    return tau - np.mean(tau)

def contact_pulse_tau(tau, tau_star, sigma_tau):
    return np.exp(-0.5 * ((tau - tau_star) / sigma_tau) ** 2)

def build_warped_signal(t, dt, params):
    w = speed_w(t, **{k: params[k] for k in ["t_star", "a", "k", "B", "p", "c"]})
    tau = tau_from_w(t, w, dt)
    tau_star = tau[np.argmin(np.abs(t - params["t_star"]))]
    D_tau = contact_pulse_tau(tau, tau_star, params["sigma_tau"])
    return z(D_tau)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/PerBin_Averages_Christian_v2.csv",
                    help="CSV with Start, End, Obs_avg, Exp_best_avg")
    ap.add_argument("--out", default="data/AC_WARP_grid.csv", help="grid CSV output path")
    args = ap.parse_args()

    perbin = pd.read_csv(args.data).sort_values("Start").copy()
    perbin["Residual"] = perbin["Obs_avg"] - perbin["Exp_best_avg"]
    t = (perbin["Start"] + perbin["End"]) / 2.0
    t = t.values.astype(float)
    dt = float(np.median(perbin["End"] - perbin["Start"]))
    R = z(perbin["Residual"].values)

    # Parameter ranges (compact, adjustable)
    grid = dict(
        t_star=list(range(-200, 125, 5)),  # Scan from 200 BCE to 120 CE in 5-year bins
        a=[0.5, 1.0, 2.0],
        k=[0.004, 0.006, 0.02],
        B=[0.0, 0.8, 1.6],
        p=[0.8, 1.0],
        c=[10.0, 25.0],
        sigma_tau=[25.0],
    )
        # You can include more if you want to explore pulse width too

    rows = []
    for t_star in grid["t_star"]:
        for a in grid["a"]:
            for k in grid["k"]:
                for B in grid["B"]:
                    for p in grid["p"]:
                        for c in grid["c"]:
                            for sigma_tau in grid["sigma_tau"]:
                                params = dict(t_star=t_star, a=a, k=k, B=B, p=p, c=c, sigma_tau=sigma_tau)
                                D = build_warped_signal(t, dt, params)
                                r0 = float(np.corrcoef(R, D)[0, 1])
                                rows.append({**params, "r0": r0})

    pd.DataFrame(rows).sort_values("r0", ascending=False)\
        .to_csv(args.out, index=False)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
