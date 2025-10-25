#!/usr/bin/env python3
"""
linear_vs_warp_core.py
Original OT→NT approach/contact/aftershock analysis (no Daniel).

Reads:
  data/PerBin_Averages_Christian_v2.csv

Writes:
  data/AC_WARP_grid.csv
  data/AC_WARP_summary.json
  data/figures/AC_WARP_alignment.png
  data/figures/AC_WARP_speed.png
  data/figures/AC_WARP_tau.png
  data/figures/AC_WARP_comparison.png
"""

import argparse, json, os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

# ---------- utilities ----------
def z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x))/s if s > 0 else np.zeros_like(x)

def perm_p(R: np.ndarray, D: np.ndarray, n: int = 3000, seed: int = 123) -> tuple[float, float]:
    """Circular-shift permutation p-value for zero-lag correlation."""
    R = np.asarray(R, float); D = np.asarray(D, float)
    r_obs = float(np.corrcoef(R, D)[0, 1])
    rng = np.random.default_rng(seed); npts = len(D); count = 0
    for _ in range(n):
        k = rng.integers(0, npts)
        rp = float(np.corrcoef(R, np.roll(D, k))[0, 1])
        if abs(rp) >= abs(r_obs): count += 1
    p = (count + 1) / (n + 1)
    return r_obs, p

def gaussian_linear(t, center=33.0, sigma=25.0):
    g = np.exp(-0.5 * ((t - center) / sigma) ** 2)
    return z(g)

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
    return z(D_tau), tau, w

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/PerBin_Averages_Christian_v2.csv",
                    help="CSV with Start, End, Obs_avg, Exp_best_avg")
    ap.add_argument("--outdir", default="data", help="output directory")
    ap.add_argument("--figdir", default="data/figures", help="figures directory")
    ap.add_argument("--perms", type=int, default=3000, help="permutations for p-values")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    perbin = pd.read_csv(args.data).sort_values("Start").copy()
    perbin["Residual"] = perbin["Obs_avg"] - perbin["Exp_best_avg"]
    t = ((perbin["Start"] + perbin["End"]) / 2.0).values.astype(float)
    dt = float(np.median(perbin["End"] - perbin["Start"]))
    R = z(perbin["Residual"].values)

    # --- Linear baseline (small sigma grid) ---
    sigma_grid = [20.0, 25.0, 40.0, 60.0]
    lin_rows = []
    best_lin = (None, -9.0)  # (sigma, r)
    for s in sigma_grid:
        D_lin = gaussian_linear(t, 33.0, s)
        r = float(np.corrcoef(R, D_lin)[0, 1])
        lin_rows.append({"sigma": s, "r": r})
        if r > best_lin[1]:
            best_lin = (s, r)
            D_lin_best = D_lin
    r_lin, p_lin = perm_p(R, D_lin_best, n=args.perms, seed=9)

    # --- Warp grid (compact) ---
    grid = dict(
        t_star=[33.0],
        a=[0.5, 1.0, 2.0, 3.0],
        k=[0.004, 0.006, 0.02, 0.05],
        B=[0.0, 0.8, 1.6],
        p=[0.8, 1.0, 1.2],
        c=[10.0, 25.0],
        sigma_tau=[20.0, 25.0, 40.0, 60.0],
    )
    rows = []
    best = {"r": -9.0}
    for t_star in grid["t_star"]:
        for a in grid["a"]:
            for k in grid["k"]:
                for B in grid["B"]:
                    for p in grid["p"]:
                        for c in grid["c"]:
                            for sigma_tau in grid["sigma_tau"]:
                                params = dict(t_star=t_star, a=a, k=k, B=B, p=p, c=c, sigma_tau=sigma_tau)
                                D_warp, tau, w = build_warped_signal(t, dt, params)
                                r0 = float(np.corrcoef(R, D_warp)[0, 1])
                                rows.append({**params, "r0": r0})
                                if r0 > best["r"]:
                                    best.update({"r": r0, "params": params, "D": D_warp, "tau": tau, "w": w})

    grid_path = os.path.join(args.outdir, "AC_WARP_grid.csv")
    pd.DataFrame(rows).sort_values("r0", ascending=False).to_csv(grid_path, index=False)

    # Permutation on best warp
    r_warp, p_warp = perm_p(R, best["D"], n=args.perms, seed=11)

    # --- Save summary ---
    summary = {
        "linear": {"best_sigma_calendar": best_lin[0], "r": float(r_lin), "p": float(p_lin)},
        "warp":   {"best_params": best["params"], "r": float(r_warp), "p": float(p_warp)},
        "delta_r": float(r_warp - r_lin)
    }
    with open(os.path.join(args.outdir, "AC_WARP_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("Summary:", json.dumps(summary, indent=2))

    # --- Figures ---
    # Alignment
    plt.figure(figsize=(10, 4))
    plt.plot(t, R, label="Residual (z)", linewidth=2)
    plt.plot(t, best["D"], label="Warp pulse (z)", linewidth=2)
    plt.axvline(best["params"]["t_star"], linestyle="--", alpha=0.6)
    plt.title(f"Alignment (lag 0): warp r={r_warp:.3f}, p≈{p_warp:.3f}; linear r={r_lin:.3f}, p≈{p_lin:.3f}")
    plt.xlabel("Year"); plt.ylabel("z-score"); plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.figdir, "AC_WARP_alignment.png"), dpi=150)
    plt.close()

    # Warp speed w(t)
    plt.figure(figsize=(10, 4))
    plt.plot(t, best["w"], linewidth=2)
    plt.axvline(best["params"]["t_star"], linestyle="--", alpha=0.6)
    plt.title("Warp speed w(t): approach → contact → aftershock tail")
    plt.xlabel("Year"); plt.ylabel("w(t)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.figdir, "AC_WARP_speed.png"), dpi=150)
    plt.close()

    # Tau mapping
    plt.figure(figsize=(10, 4))
    plt.plot(t, best["tau"], linewidth=2)
    plt.axvline(best["params"]["t_star"], linestyle="--", alpha=0.6)
    plt.title("Prophetic time mapping τ(t) (cumulative of w)")
    plt.xlabel("Year"); plt.ylabel("τ(t) [shifted]")
    plt.tight_layout()
    plt.savefig(os.path.join(args.figdir, "AC_WARP_tau.png"), dpi=150)
    plt.close()

    # Comparison bar
    plt.figure(figsize=(7, 4))
    x = np.arange(2)
    vals = [r_lin, r_warp]
    labels = ["Linear (best σ)", "Warp (best)"]
    plt.bar(x, vals)
    plt.xticks(x, labels, rotation=0)
    plt.title("Correlation at lag 0: Linear vs Warp")
    plt.ylabel("Pearson r")
    plt.tight_layout()
    plt.savefig(os.path.join(args.figdir, "AC_WARP_comparison.png"), dpi=150)
    plt.close()

if __name__ == "__main__":
    main()
