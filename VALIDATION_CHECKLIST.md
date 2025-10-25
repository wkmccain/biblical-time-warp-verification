
# Validation Checklist (for external reviewers)

Use this checklist to replicate and audit the results.

## Inputs
- [ ] `PerBin_Averages_Christian_v2.csv` present and readable
- [ ] Inspect residual column: `Residual = Obs_avg − Exp_best_avg`

## Baseline (linear)
- [ ] Compute Gaussian at 33 CE, σ=25 y -> z‑score
- [ ] Pearson r with residuals ∈ [0.67, 0.69]
- [ ] Permutation p ∈ [0.02, 0.05]

## Warp model
- [ ] Define speed `w(t)` with parameters (t★=35, a=1.0, k=0.08, tail_B=0.5, tail_p=0.7, c=25.0)
- [ ] Compute τ(t) = cumulative sum of w(t)*Δt, mean‑centered
- [ ] Build pulse in τ around τ(t★) with σ=25, z‑score
- [ ] Pearson r with residuals ∈ [0.72, 0.74]
- [ ] Permutation p ∈ [0.01, 0.02]

## Sensitivity
- [ ] Try σ ∈ {20, 30, 40} for linear; confirm warp remains higher
- [ ] Vary `c` ∈ {10, 25, 50}; confirm r remains ≥ 0.70 for warp
- [ ] Try 10‑year and 50‑year bins (if available); confirm qualitative advantage

## Transparency
- [ ] Post code and the **exact parameters** used
- [ ] Share the produced grid file or parameter sweep results
- [ ] Declare limitations: structure ≠ cause; dependency on chosen residual series; binning effects

If your results fall within the listed bands, you have successfully validated the core claim: a **non‑linear approach/contact time mapping** centered near 33 CE explains the intensity curve **better** than linear time.
