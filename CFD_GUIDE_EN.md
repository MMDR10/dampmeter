# Reading Aerodynamic Damping from Single-Point Probes: A Dual-Operator Guide for Unsteady CFD

> How to use a **single-channel scalar-series dual-operator primitive** in unsteady CFD
> (URANS / DES / LES): the **IDV damping caliper** (measure aerodynamic damping γ = 1/τ)
> and the **T blind-state gatekeeper** (flag when readings cannot be trusted).
> Implementation: `pip install dampmeter` · Code + bilingual research note:
> [GitHub MMDR10/dampmeter](https://github.com/MMDR10/dampmeter) · Archive: [10.5281/zenodo.24335208](https://doi.org/10.5281/zenodo.24335208)
>
> **Evidence statement (up front):** the IDV core (four-estimator convergence + shuffle
> null) is validated on a geophysical system, and the ring-down caliber on mechanical
> vibration data. **The CFD interface itself has ZERO validated cases** — this is a
> forward-looking methodology proposal. What follows is "how to wire it up" and "how to
> break it yourself", not a promise that it works.

---

## 1. Two CFD pain points

**Pain 1 — how do you measure damping?**
An unsteady run leaves you monitor-point time series (wake pressure p(t), velocity u'(t), lift coefficient Cl(t)). To get the *aerodynamic damping* of the operating point (the key quantity for flutter), the classic routes are free-decay experiments or generalized-aerodynamic-force modal projections — expensive, heavy, and not always feasible from single-point data.

**Pain 2 — how do you judge convergence?**
Residuals plateauing ≠ the physical transient is over. Statistics still drifting while you post-process is the classic fake-convergence trap of unsteady CFD.

The dual operators target exactly these: **IDV answers "how much damping, and how stable"; T answers "when will these readings lie".**

## 2. The operators in 30 seconds

- **IDV (damping caliper):** one scalar series in, effective damping time τ (γ = 1/τ) out. Four independent estimators (ACF e-fold / distance-correlation e-fold / impulse ring-down / post-drop recovery) each measure the same τ; at least three must converge within ±30% before any number is returned — **convergence is an identity check ("the constant is really there"), not goodness-of-fit**. A shuffle null is included: destroy temporal structure, re-measure, τ must collapse.
- **T (gatekeeper):** produces no signal; its output is a negative alarm — "when not to trust directional readings". ⚠️ Honest note: T's semantics are validated on 2D trajectory domains; the 1D scalar interface (CFD probes) is **not yet derived** — see usage 2.

Compute class: KB RAM, millisecond tests, 610 samples for a reliable-tier reading.

## 3. Usage 1: aerodynamic damping extraction (flutter-boundary scans) ⭐

Classic aeroelasticity: as flow speed U approaches the flutter critical, effective damping
softens (γ→0). IDV offers a lightweight route — read effective γ straight off a single
probe, sweep operating points, track γ(U).

**Wiring:**

1. **Signal:** pick 1–3 monitors per condition (Cp(t), u'(t), Cl(t) all fine); confirm a **stationary window**.
2. **Detrend:** remove mean + linear trend (fit on the anomaly). Skip this and the ACF plateau never reaches 1/e — τ comes back NaN. (We hit this bug ourselves.)
3. **Sampling:** do NOT use the CFD time step blindly. You need ≥20–30 samples per τ and ≥610 samples total (see Pitfall 1).
4. **Measure:**

```python
import numpy as np
import dampmeter as dm

x = probe_series - np.mean(probe_series)   # detrended single-point series
est = dm.measure(x, dt=dt_sample)          # dt_sample in physical time units

est.tau          # damping time constant
est.gamma        # 1/tau — the softening indicator
est.converged    # True = at least 3 estimators agree within ±30%
est.as_dict()    # per-estimator taus, warnings, data tier
```

5. **Sweep** (U, or angle of attack α, or reduced frequency k):

```python
results = {}
for U in speeds:
    x = detrend(load_probe(case=U))
    est = dm.measure(x, dt=dt_sample)
    if est.converged:
        results[U] = est.gamma
# gamma falling with U = damping softening = the classic flutter precursor
```

6. **Null check (mandatory):**

```python
tau_null = dm.shuffle_null_ringdown(x, blocks=blocks)
# must collapse far below the measured tau — otherwise you are reading
# a mathematical necessity, not physics
```

**Interpretation discipline:**

- IDV measures **effective damping** (dissipative time constant), **not stiffness** (well curvature ω₀²). Flutter involves both; γ softening is an *indicator*, not a sufficient condition.
- "Damage/condition change → how γ moves" has no universal sign (dual-axis damage hypothesis) — calibrate direction on a ground-truth benchmark before extrapolating.
- Effective γ only — no mechanism decomposition.

## 4. Usage 2: true limit-cycle convergence (🔭 forward-looking)

**Idea:** during the initial transient, statistics drift; once a developed state (LCO / developed shedding) is reached, drift extinguishes. Use rolling-window τ stability ("red light off") as a *physical* convergence criterion against residual-plateau fake convergence:

- run `dm.measure()` over consecutive windows; τ converging within ±30% across windows + stable event structure = transient over
- inverse use: τ still drifting = red light on = **don't post-process this segment yet**

⚠️ T's "state detachment" semantics are only validated on 2D trajectory domains; the 1D red-light definition is not yet derived. This usage is a **forward-looking proposal** — run your own sensitivity checks (Section 6) before relying on it.

## 5. Usage 3: grid/scheme physical convergence

Residual convergence ≠ physical convergence. Run the same case across grid densities / numerical schemes and compare γ from the same monitor series — **if γ is grid-insensitive (within ±30%), that is stronger evidence of physical convergence than any residual plot.**

## 6. Five pitfalls (empirical lessons)

1. **610 samples ≠ 610 time steps.** If τ = 10 ms and you sample every 0.01 ms, 610 samples cover 6 ms (0.6 τ) — fatal. The window must span multiple τ (aim ≥30–50 τ); downsample if needed.
2. **No detrend → τ = NaN.** Mean offsets lift the ACF plateau below 1/e never being reached.
3. **Envelope must be the Hilbert analytic envelope.** Using |band-passed| quantizes τ at the oscillation half-period (empirical artifact: 0.078 ms). `dampmeter.ringdown_envelope()` ships the corrected caliber (Hilbert + normalized MAD threshold).
4. **Don't transplant thresholds across cases.** Re-derive per case with sensitivity checks (cross-domain threshold transplant failed 4/4 in our tests).
5. **Never report from a single estimator.** Convergence + null is the bar; `est.converged == False` means "this window yields no number" — respect it.

## 7. How to break it (falsification protocol)

The CFD interface has no validated cases — the most valuable thing you can do is beat it:

1. Pick a benchmark with a published flutter boundary (AGARD 445.6 wing, Goland beam, NACA 0012 section).
2. Sweep URANS/DES across the critical speed, take single-point probe series at each point.
3. **Pre-register the prediction:** γ(U) should soften monotonically and approach zero near the published critical; direction should agree with classical aerodynamic-damping estimates.
4. Score it: if it holds, the interface earns empirical status; if it fails, note which link broke (window too short? stationarity violated? effective γ vs modal damping mismatch?) — a falsification record is valuable either way.

**Weapons, not conclusions.**

## 8. Evidence grades (honest)

| Claim | Grade |
|------|---------|
| IDV four-estimator convergence measures intrinsic τ (+ shuffle null) | ✅ cross-validated (geophysics + mechanical ring-down) |
| 610-point reliable floor, KB-RAM / ms compute | ✅ cross-validated (2,800 bootstraps) |
| Hilbert envelope / MAD caliber necessity | ✅ cross-validated (reverse lesson) |
| γ(U) softening as flutter indicator via IDV | 🔭 theoretically coherent (classical aeroelasticity), IDV route untested |
| Limit-cycle true-convergence criterion | 🔭 forward-looking (T 1D interface not derived) |
| Grid physical-convergence comparison | 🔭 forward-looking (reasonable extension, untested) |

## Links

- Install: `pip install dampmeter`
- Code + bilingual research note: https://github.com/MMDR10/dampmeter
- Archive (ZH + EN + code): https://doi.org/10.5281/zenodo.24335208

*tygtDc, Deep Research · Ma Ka Po (CC BY 4.0)*
*AI contribution disclosure: MKP provided original intuition and research direction; tygtDc, Deep Research performed structural synthesis and writing; Gemini served as dialogue partner on writing strategy.*
