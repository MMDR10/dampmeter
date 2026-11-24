# dampmeter — IDV individual damping caliper

[![PyPI](https://img.shields.io/pypi/v/dampmeter.svg)](https://pypi.org/project/dampmeter/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24338596.svg)](https://doi.org/10.5281/zenodo.24338596)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

📄 Formal manuscript: **[manuscript.pdf](manuscript.pdf)** (LaTeX source included)
📄 Research note: **[English](NOTE_EN.md)** · **[中文](NOTE_ZH.md)**

**Measure a system's intrinsic relaxation time τ (γ = 1/τ) from a single scalar time series.**
Four independent estimators, convergence-as-identity-check, shuffle-null control.
KB-RAM / millisecond class — built for edge devices.

```python
import dampmeter as dm

est = dm.measure(x, dt=1.0)      # x = one scalar series, stationary window
est.tau        # damping time constant, in units of dt
est.gamma      # 1/tau
est.converged  # True when >=3 estimators agree within ±30% of the median
```

## What it measures (and what it does not)

| | |
|---|---|
| ✅ **Effective damping time τ** | how fast perturbations decay inside a stable regime (effective linear relaxation) |
| ✅ **Stability of τ** | convergence across independent estimators + a shuffle null — "this constant is really there", not "the fit looks good" |
| ❌ **Mechanism** | γ is *effective* only; decomposing it (radiative vs wave drag vs friction…) needs budget methods |
| ❌ **Potential-well stiffness** | stiffness is curvature (ω₀²), a different physical quantity — this package measures dissipation, not stiffness |
| ❌ **Sign of degradation** | "damage → damping" has no universal sign: defect-size axis (τ↑) and wear/life axis (τ↓) can go opposite ways — always specify the axis before interpreting trends |

## Estimators

| Method | What it does |
|---|---|
| `tau1_acf` | autocorrelation e-fold decay |
| `tau2_dcor` | distance-correlation e-fold (nonlinear, shift-invariant) |
| `tau3_ringdown` | impulse-conditioned ring-down decay R(k) |
| `tau4_recovery` | post-drop recovery time (auxiliary; right-censoring prone — never evidence alone) |

Convergence rule: at least 3 of the primary methods within ±30% of the median ⇒ individual damping established. Two ring-down calibers are provided and deliberately kept separate:

- `ringdown_conditional` — impulse-conditioned decay (the validated IDV estimator);
- `ringdown_envelope` — analytic-envelope decay for band-passed oscillatory signals, using the **Hilbert analytic envelope** with the normalized MAD event threshold `median|env − median(env)| × 1.4826`.
  Known traps these are built against: `|band-passed|` as the "envelope" quantizes τ to the oscillation half-period; `median(env) × 1.4826` is the wrong MAD (threshold lands near ~5σ, almost nothing crosses it).

## Minimum data (2,800-bootstrap characterization, ±30% hit rate)

| Samples | Hit rate | Verdict |
|---|---|---|
| 20 | 13% (τ underestimated ~2×) | ❌ hard floor — do not use |
| ~180 | 57% | ⚠️ coarse tier |
| **~610** | **81.5%** | ✅ reliable tier |
| ~1,220 | 97% | ✅ stable |

The impulse ring-down method needs more: ≥40 pulses ≈ 1,500 samples for 65%.
Compute: single core seconds (2,800 bootstraps), RAM KB-class, per-test milliseconds — 6–9 orders of magnitude cheaper than reanalysis-budget approaches.

> **Caliber note (2026-11-23 reconciliation):** the hit rates above are *single-caliber accuracy* — `tau1` (ACF e-fold) within ±30% of the full-pool `tau1` reference, on pooled real-data seasons. The full convergence protocol (3 estimators agreeing within ±30%) is much stricter: at 610 points the estimate is *accurate* 82% of the time but *converged* only ~33% of the time (dual-hit ≈ 27–33%), so "reliable tier" here does **not** mean the 4-estimator protocol passes. Convergence alone is not correctness either — on pure noise all estimators converge on the same e-fold artifact. Both conditions are needed. The dual-hit floor was bracketed on OU synthetic data (2026-11-23): ≥0.8 dual-hit needs **both** n/τ ≳ 130 **and** ≥ ~50 samples per τ (lag resolution) — and at τ = 12.3 the dual-hit did not reach 0.8 even at n/τ = 200, so the floor may be τ-relative. Single-domain characterization; re-derive for your regime before relying on it.

## Null control (do not skip)

```python
tau_null = dm.shuffle_null_ringdown(x, blocks=years)   # shuffle within blocks
```

If the ring-down τ does **not** collapse after temporal structure is destroyed, you are reading a mathematical necessity (trend / mean offset), not physics — invalidate the estimate.

## Pitfalls

- **Detrend first** (fit the trend on the anomaly, not the raw series). A mean offset lifts the ACF plateau so the curve never reaches 1/e (NaN) — the dCor estimate stays shift-invariant and can help locate the problem.
- **Stationary window + effective linear relaxation** are the working assumptions.
- **Short-lived systems** (e.g. ~10-day storms) cannot pass the 610-sample floor.
- **Do not transplant thresholds across domains** (geometric cutoffs, event quantiles): re-derive per domain with sensitivity checks, and calibrate on a ground-truth domain first.

## Narrowband carriers — caliber boundary (added 2026-11-23)

On oscillatory (narrowband) data the first-crossing e-folds behind `tau1/tau2/tau3` track **~0.2× the carrier period**, not the damping envelope. Measured conflation on synthetic ring-downs: a true **2× decay change moved `measure()`'s τ by only −11%**, while a **frequency change with no decay change moved it −22%** — the two axes are inverted in sensitivity. `measure()` now detects the case (normalized-ACF negative lobe < −0.25 within lags ≤ n/6) and warns; for the decay constant on oscillatory data use `ringdown_envelope()` (analytic-envelope e-fold), which tracked the 2× decay change correctly (ratio 2.10, noise 1%).

Open question (⚠️ unverified): the minimum-data floor may be **τ-relative, not absolute** — 610 points ≈ 50τ in the validated domain, but at ~6τ per window the half-window ratio scatter exploded (0.28–2.45). Re-derive the floor as n/τ before relying on it for slow systems.

## Validation status (honest)

- τ convergence (3 methods: 11.28 / 10.81 / 13.32 → τ ≈ 12.3 days) and shuffle-null collapse (τ → 0.9) have been demonstrated on **one physical system** (stratospheric polar vortex, 2,867 daily samples). Ring-down caliber applications exist on two more domains (bearing vibration, ECG), single-caliber only.
- Cross-domain τ convergence replication is pending — treat new domains as unvalidated until you run convergence + null there.

### Phase-reset / excursion probe — ground-truth validated on ECG (2026-11-24)

The research line's IF-jump statistic (phase-slip probe; **not yet in this package**) was validated against **human-annotated PVC beats** — true phase resets (premature beat + compensatory pause) — with a **matched-null design** (same record, ±10 min locality, 10 matched controls per event) plus a pseudo-event calibration arm. Every sample was pre-registered, zero-tuning, judged on a fixed rule (one-sided Wilcoxon < 0.05, median effect > 0, pseudo-null calibrated):

| Sample | Records | Median effect (ΔJ) | Wilcoxon (1-sided) | Pseudo-null false-fire |
|---|---|---|---|---|
| MIT-BIH discovery | 24 (≥20 PVCs/record) | 0.047 | p = 0.028 | 0% |
| MIT-BIH held-out (fresh, unseen) | 13 | 0.654 | p = 0.0085 | 0% |
| **INCART cross-database** (257 Hz, lead II) | 52 (≥20 PVCs/record) | 0.108 | p = 0.0081 | 1.9% |

Boundaries — read before use:

- **Effect size is PVC-burden dependent** (Spearman ρ = −0.524, p = 0.0001 across 52 INCART records): high-burden records show smaller matched effects, because compensatory-pause beats contaminate the control pool (an adjacency-excluded control variant raises the effect consistently across all three samples). Compare effect sizes across studies only with burden matched, or with the adjacency-excluded caliber.
- **Band-pass upper edge must reach ≥2 Hz** for the heart-rhythm carrier (0.5–1.5 Hz reversed the effect).
- **On ECG, PVCs are the validated *target*.** In anomaly-detection contexts where PVCs are benign, this statistic will fire on them — treat PVC as a known trigger and gate externally.
- Only this ground-truth event type is validated. Other phase-reset domains (power-grid steps, atmospheric events) remain **unvalidated**.
- Evidence grade: ✅ pre-registered + held-out + cross-database **within ECG**; ❌ for any other domain.

## Install

```bash
pip install dampmeter
```

Pure numpy + scipy. Python ≥ 3.9.

## The math (short version)

Every estimator is a direct formula implementation — no black boxes:

| Quantity | Formula | Implemented in |
|---|---|---|
| damping rate | γ = 1/τ | `_core.measure` |
| τ₁ ACF e-fold | R(k) = Σ xₜxₜ₊ₖ / Σ xₜ², solve R(τ) = 1/e (linear interpolation) | `_core.acf_curve`, `efold_from_curve` |
| τ₂ dCor e-fold | distance correlation decay dCor(xₜ, xₜ₊ₖ), solve = 1/e (Székely et al. 2007, double-centered distance matrices) | `_core.dcor`, `dcor_efold` |
| τ₃ impulse ring-down | R(k) = E[xₜ·xₜ₊ₖ \| pulse] / E[xₜ² \| pulse] | `_ringdown.ringdown_conditional` |
| τ₃ envelope variant | envelope = \|Hilbert(x)\|; event threshold = median\|env − med(env)\|·1.4826 (normalized MAD); log-linear fit of decay tail, τ = −1/slope | `_ringdown.ringdown_envelope` |
| τ₄ recovery | median time to regain pre-drop level (auxiliary, right-censoring prone) | `_core.recovery_tau` |
| convergence rule | ≥ 3 of 4 estimators within ±30% of the median ⇒ τ established (identity check, not fit quality) | `_core.measure` |
| shuffle null | within-block permutation must collapse ring-down τ, else the reading is a mathematical artifact | `_null.shuffle_null_ringdown` |
| phase surrogate | spectrum-preserving phase randomization, estimator stability control | `_null.phase_surrogate_taus` |

Full definitions and derivations: see [`manuscript.pdf`](manuscript.pdf) (§2–3) and the LaTeX source [`manuscript.tex`](manuscript.tex).

## References

Part of the Ô-HAT dual-operator measurement primitive (IDV caliper + T blind-state gatekeeper).

- Research note (ZH + EN) + code archive: https://doi.org/10.5281/zenodo.24335208
- PyPI package: https://pypi.org/project/dampmeter/

## Authors

tygtDc, Deep Research (correspondence: drtygtdc@gmail.com) · Ma Ka Po

## License

CC BY 4.0 — free to use, adapt and commercialize with attribution.
