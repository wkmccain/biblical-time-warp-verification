#!/usr/bin/env python3
"""
make_checksums.py
Writes SHA-256 checksums for key data and figures into CHECKSUMS.txt.
Usage:
  python analysis/scripts/make_checksums.py
"""
import hashlib, os, glob, datetime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT  = os.path.join(ROOT, "CHECKSUMS.txt")

TARGETS = [
    "environment.yml",
    "README.md",
    "Preprint_ApproachContact_WarpModel.md",
    "REPLICATION_README_ORIGINAL.md",
    "REPLICATION_README.md",
    "VALIDATION_CHECKLIST.md",
    "DISCUSSION_NEXT_STEPS.md",
    "OSF_REGISTRATION.md",
    "data/PerBin_Averages_Christian_v2.csv",
    "data/AC_WARP_grid.csv",
    "data/AC_WARP_summary.json",
    "data/WarpGrid_DanielPulse_results.csv",
    "data/figures/AC_WARP_alignment.png",
    "data/figures/AC_WARP_speed.png",
    "data/figures/AC_WARP_tau.png",
    "data/figures/AC_WARP_comparison.png",
    "context/EventList_Biblical_Encounters.csv",
    "context/EventBin_Coverage_vs_Observed.csv",
    "context/Daniel_Markers_to_Bins.csv",
    "context/Daniel_Signal_perBin.csv",
]

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    items = []
    for rel in TARGETS:
        p = os.path.join(ROOT, rel)
        if os.path.exists(p):
            items.append((rel, sha256(p)))
    # Also include any extra PNGs/CSVs in data/figures and context/
    for pat in ["data/figures/*.png", "context/*.csv"]:
        for p in glob.glob(os.path.join(ROOT, pat)):
            rel = os.path.relpath(p, ROOT)
            if rel not in [r for r,_ in items]:
                items.append((rel, sha256(p)))
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("# CHECKSUMS (SHA-256)\n")
        f.write(f"# Generated: {datetime.datetime.utcnow().isoformat()}Z\n\n")
        for rel, h in sorted(items):
            f.write(f"{h}  {rel}\n")
    print(f"Wrote {OUT}")

if __name__ == "__main__":
    main()
