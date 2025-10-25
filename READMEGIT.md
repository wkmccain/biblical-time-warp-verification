# Prophetic Intensity Aligns with a Natural Pulse at ~33 CE

**Author:** Wiley K. McCain  
**Preregistration:** https://doi.org/10.17605/OSF.IO/EXYWU  
**Licensing:** Code — Prosperity Public License 3.0.0 (non-commercial) · Data/Docs — CC BY-NC 4.0  
**Contact:** <your-email-or-site>

---

## Overview

We test whether the biblical timeline shows a **burst pattern**—**approach → contact → aftershock**—similar to foreshock/mainshock/aftershock in seismology.  
A neutral, prereg-style pipeline compares a **no-warp (linear)** baseline to a **time-warp** model that concentrates weight near a candidate contact year and allows a heavy-tailed unwind.

Across 10k Monte Carlo resamplings of event date ranges (25-year bins), the **best alignment centers near ~33 CE** in the biblical set. Running the **same pipeline** on Greek/Roman/Han controls shows **no comparable advantage**. Exact statistics for your run live in `data/AC_WARP_summary.json`.

> Note: The 33 CE center is **discovered**, not hard-coded. Results may vary slightly by dataset revision/seed; verify from the summary JSON rather than a fixed number in this README.

---

## Quickstart (reproduce in minutes)

```bash
conda env create -f environment.yml
conda activate timeline

# 1) Monte Carlo per-bin averages (10,000 draws, fixed seed)
python scripts/monte_carlo_10k.py   --events context/EventList_Biblical_Encounters.csv   --out    data/PerBin_Averages_Christian_v2.csv   --draws  10000 --seed 2025

# 2) Optional baseline figure
python scripts/build_baseline_plot.py   --perbin data/PerBin_Averages_Christian_v2.csv

# 3) Linear vs Warp comparison (writes grid, summary, figures)
python scripts/linear_vs_warp_core.py   --data   data/PerBin_Averages_Christian_v2.csv   --outdir data --figdir data/figures --perms 10000
```

**Outputs to check**
- `data/AC_WARP_summary.json`  ← r, p, and best-fit parameters  
- `data/AC_WARP_grid.csv`  
- `data/figures/AC_WARP_alignment.png`, `..._speed.png`, `..._tau.png`, `..._comparison.png`  
- `data/figures/BASELINE_composite.png` (if you ran step 2)

---

## Controls (recommended)

Use the same pipeline to show the effect is **not** generic:

```bash
# Greek (raw events, year-only OK → start=end=year)
python scripts/monte_carlo_10k.py --events context/Greek_Events.csv   --out data/PerBin_Averages_Greek.csv --draws 10000 --seed 2025
python scripts/linear_vs_warp_core.py --data data/PerBin_Averages_Greek.csv   --outdir data/greek --figdir data/figures/greek --perms 10000

# Roman & Han (if per-bin files are provided instead of raw events)
python scripts/linear_vs_warp_core.py --data data/PerBin_Averages_Roman.csv   --outdir data/roman --figdir data/figures/roman --perms 10000
python scripts/linear_vs_warp_core.py --data data/PerBin_Averages_Han.csv   --outdir data/han --figdir data/figures/han --perms 10000
```

---

## Repo contents

```
/scripts/
  monte_carlo_10k.py         # builds per-bin MC averages
  build_baseline_plot.py     # makes BASELINE_composite.png
  linear_vs_warp_core.py     # runs comparisons, figures, summary/grid
  warp_grid_search.py, permutation_tests.py, make_checksums.py
/context/
  EventList_Biblical_Encounters.csv        # core input (start,end[,weight])
  EventBin_Coverage_vs_Observed.csv        # optional overlay for baseline figure
/data/ (created by scripts)
/data/figures/ (created by scripts)
```

---

## Methods snapshot

- **Uncertainty:** sample each event uniformly over `[start,end]`, 10,000 draws.  
- **Binning:** 25-year bins from −2000..+125.  
- **Baselines (≤0 CE only):** constant, 8-segment piecewise constant, and degree-5 polynomial; pick by **Poisson deviance** on 0–100 CE test window.  
- **Signal:** residual `Obs_MC − Exp_baseline` and its alignment with the warp model centered at candidate years; reported via Pearson **r**, **Δr** vs. linear, and permutation **p** (circular shifts).  
- **Integrity:** optional `make_checksums.py` to write `CHECKSUMS.txt`.

---

## Danielic overlay (exploratory)

Files such as `Daniel_Signal_perBin.csv` and `Daniel_Markers_to_Bins.csv` illustrate a **theological overlay**. They’re **not** part of the core claim and should be treated as **exploratory**; the main conclusion rests on the prereg-style pipeline above.

---

## Licensing

- **Code:** Prosperity Public License 3.0.0 (non-commercial) — see `LICENSE`  
- **Data, figures, docs:** CC BY-NC 4.0 — see `LICENSE-data`  
Commercial licensing available on request.

---

## Citation

> McCain, W. K. (2025). *Prophetic Intensity Aligns with a Natural Pulse at ~33 CE* (software & dataset).  
Include the repo URL/DOI and the preregistration DOI in references.
