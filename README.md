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

## Validation status (honest)

- τ convergence (3 methods: 11.28 / 10.81 / 13.32 → τ ≈ 12.3 days) and shuffle-null collapse (τ → 0.9) have been demonstrated on **one physical system** (stratospheric polar vortex, 2,867 daily samples). Ring-down caliber applications exist on two more domains (bearing vibration, ECG), single-caliber only.
- Cross-domain τ convergence replication is pending — treat new domains as unvalidated until you run convergence + null there.

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
