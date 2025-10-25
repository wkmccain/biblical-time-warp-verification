#!/usr/bin/env python3
"""
build_baseline_plot.py
Rebuilds the composite baseline figure:
 - Observed (MC avg) vs Expected baseline (MC avg) per 25-yr bin
 - Optional overlays: reconstructed event weights, Daniel per-bin signal
 - Annotates OT pulses and the Christ window

Inputs (required):
  data/PerBin_Averages_Christian_v2.csv     (columns: Start, End, Obs_avg, Exp_best_avg)

Optional inputs (if present):
  context/EventBin_Coverage_vs_Observed.csv (Event_weight_scaled)
  context/Daniel_Signal_perBin.csv          (Daniel_weight_max)

Outputs:
  data/figures/BASELINE_composite.png
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def z(x):
    x = np.asarray(x, float)
    s = np.nanstd(x)
    return (x - np.nanmean(x))/s if s > 0 else np.zeros_like(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perbin", default="data/PerBin_Averages_Christian_v2.csv")
    ap.add_argument("--event_cov", default="context/EventBin_Coverage_vs_Observed.csv")
    ap.add_argument("--daniel", default="context/Daniel_Signal_perBin.csv")
    ap.add_argument("--out", default="data/figures/BASELINE_composite.png")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Load per-bin baseline series
    perbin = pd.read_csv(args.perbin).sort_values("Start").copy()
    perbin["Mid"] = (perbin["Start"] + perbin["End"]) / 2.0

    # Core series
    t   = perbin["Mid"].values
    obs = perbin["Obs_avg"].values
    exp = perbin["Exp_best_avg"].values

    # Optional overlays
    event_z = None
    if os.path.exists(args.event_cov):
        ev = pd.read_csv(args.event_cov).sort_values("Start")
        if "Event_weight_scaled" in ev.columns:
            event_z = z(ev["Event_weight_scaled"].values)

    daniel_z = None
    if os.path.exists(args.daniel):
        ds = pd.read_csv(args.daniel).sort_values("Start") if "Start" in pd.read_csv(args.daniel, nrows=1).columns else pd.read_csv(args.daniel)
        if "Daniel_weight_max" in ds.columns:
            daniel_z = z(ds["Daniel_weight_max"].values)

    # Plot
    plt.figure(figsize=(12, 5))
    plt.plot(t, obs, label="Observed (MC avg)", linewidth=2)
    plt.plot(t, exp, label="Expected baseline (MC avg)", linewidth=2)

    # Optional z-scored overlays on a 2nd axis to avoid scale confusion
    ax = plt.gca()
    if event_z is not None or daniel_z is not None:
        ax2 = ax.twinx()
        if event_z is not None:
            ax2.plot(t, event_z, linestyle="--", alpha=0.6, label="Event weight (z)")
        if daniel_z is not None:
            ax2.plot(t, daniel_z, linestyle=":", alpha=0.7, label="Daniel signal (z)")
        ax2.set_ylabel("z-score (overlay)")

        # Merge legends
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc="upper left")
    else:
        ax.legend(loc="upper left")

    # Annotate OT pulses (approximate bins nearest these centers)
    def annotate_pulse(x0, txt):
        xi = t[np.argmin(np.abs(t - x0))]
        yi = obs[np.argmin(np.abs(t - x0))]
        plt.annotate(txt, xy=(xi, yi), xytext=(xi, yi + (np.nanmax(obs) - np.nanmin(obs))*0.08),
                     arrowprops=dict(arrowstyle="->"), fontsize=9)

    annotate_pulse(-520, "OT pulse: Persian Return (~520 BCE)")
    annotate_pulse(-165, "OT pulse: Maccabean Crisis (~165 BCE)")

    # Highlight Christ window 0–75 CE
    plt.axvspan(0, 75, alpha=0.12, color="grey")
    plt.text(2, np.nanmax(obs), "Christ window (0–75 CE)", fontsize=9, va="top")

    plt.title("Composite View: OT pulses → Christ mainshock → Aftershock")
    plt.xlabel("Year (CE; negative = BCE)")
    plt.ylabel("Weighted events per 25-yr bin (MC avg)")
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    plt.close()
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
