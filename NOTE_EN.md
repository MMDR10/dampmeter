# A Dual-Operator Measurement Primitive for Single-Channel Scalar Time Series: The IDV Individual Damping Caliper and the T Blind-State Gatekeeper

**tygtDc, Deep Research · Ma Ka Po**
(Correspondence: drtygtdc@gmail.com)
License: CC BY 4.0

---

## Problem definition

Traditional nonlinear dynamics reconstructs high-dimensional attractors from full-field grids (reanalysis-scale data) with heavy recomputation, is extremely noise-sensitive, and in most cases yields only *indirect inference* (cross-system convergent cross mapping / teleconnection). This note proposes a lightweight dual-operator measurement primitive that depends on **a single scalar time series only**, distilled from cross-domain experiments (stratospheric polar vortex, bearing run-to-failure, ECG, typhoon tracks):

- **IDV (primary operator) — the individual damping caliper.** It measures a system's *recovery capacity against perturbation* while inside a stable attractor, outputting the effective damping time constant τ (or γ = 1/τ). Four independent estimators converge on the same τ, with a shuffle-null self-check — convergence is not "the fit looks good"; it is an identity check that "this constant is really there".
- **T (secondary operator) — the blind-state gatekeeper.** It produces no signal of its own; its output is a negative alarm by design. It judges *when directional readings cannot be trusted*. Cross-domain evidence (typhoon recurvature, ECG ectopic beats, polar vortex wave-2 splitting) shows the unified physical cause behind the red light is **the system detaching from its previous stable attractor state (state detachment / topological rearrangement)** — not high energy or external noise per se.

Both operators are the same species: **meta-operators** — they do not measure the phenomenon itself but the *trustworthiness of measurement*. IDV answers "how large and how stable is this system's recovery constant"; T answers "when will these readings lie to you".

> **Terminology note.** "Recovery capacity / effective damping" here means the dissipative time constant τ (effective damping). It is **not** potential-well stiffness (curvature, ω₀²), which must be measured independently from oscillation frequency and is outside this primitive. If "rigidity" is used as an analogy, this distinction must be stated.

---

## 1. IDV: direct measurement of attractor recovery damping (the core)

### 1.1 Definition and estimators

Within a stationary window (assumptions: stationary window + effective linear relaxation), four independent estimators each measure the same τ from a system's own scalar series:

| Estimator | Caliber |
|---|---|
| τ₁ | linear autocorrelation (ACF) e-fold |
| τ₂ | nonlinear distance-correlation e-fold (shift-invariant; can cross-locate implementation bugs) |
| τ₃ | transient impulse ring-down (Hilbert analytic envelope + normalized MAD threshold) |
| τ₄ | post-drop recovery (right-censoring prone; auxiliary only) |

**Convergence rule:** at least three methods within ±30% of the median ⇒ individual damping established.
**Shuffle-null control:** destroy temporal structure by within-block shuffling and re-measure; the ring-down τ must collapse — if it does not, you are reading a mathematical necessity, not physics; invalidate.

### 1.2 Empirical baseline (stratospheric polar vortex, 2,867 daily u60 samples)

- Three-method convergence: τ₁ **11.28** / τ₂ **10.81** / τ₃ **13.32** days, all within median **12.3 ±30%** → **γ ≈ 0.081/day**
- Shuffle null: ring-down τ collapses to **0.9 days** ✅ (real physical structure)
- Conclusion: cross-system CCF/CCM indirect inference can retire — **damping can be measured directly**

### 1.3 Minimum-data and compute characterization (2,800 bootstraps)

| Data | ±30% hit rate | Verdict |
|---|---|---|
| 20 points | 13% (τ systematically underestimated ~2×) | ❌ hard floor, unusable |
| ~180 points | 57% | ⚠️ coarse tier |
| **~610 points** | **81.5%** | ✅ **reliable-tier floor** |
| ~1,220 points | 97% | ✅ stable |

- The impulse ring-down method needs more: ≥40 pulses ≈ **1,500 points** for 65%
- Compute: single core **17 s**, RAM **KB-class**, per-test **milliseconds** — 6–9 orders of magnitude below reanalysis budgets; minimal build = ACF + ring-down only (O(n·k)); dCor is the single heavy step and can be dropped
- Suitable for paleoclimate proxies (ice cores / tree rings / speleothems) and edge devices

### 1.4 Dual-axis damage hypothesis — the reading is stable, the interpretation must name the axis

"Damage → how damping changes" has **no universal sign**; at least two axes run in opposite directions (cross-domain measurements):

- **Defect-size axis** (controlled severity, CWRU): larger defects ring down slower = damping drops (τ↑), 11/11 cells ρ>0 (p≈0.001), monotonic 1.1→1.5→2.1 ms
- **Wear / life axis** (run-to-failure, XJTU-SY): τ *shortens* with degradation (13/15 negative, p=0.0074) = damping rises, suspected cause debris/lubricant-degradation as a new damping source; the artifact hypothesis was excluded by strict event-isolation re-measurement (p=0.0037)

**Usage:** the IDV reading is trustworthy, but narratives like "degradation ⇒ damping weakens" must name which axis is being walked, with the direction calibrated on a ground-truth domain first.

### 1.5 Boundaries

- Output is **effective γ**: "how much, how stable" only — no mechanism decomposition (radiative vs wave-drag budgets need other methods)
- Assumptions: stationary window + effective linear relaxation; short-lived systems (~10-day typhoons) cannot pass the 610-point floor
- **τ-convergence evidence exists on one physical domain only** (polar vortex); bearing/ECG are single-caliber ring-down applications — do not generalize before cross-domain replication
- Implementation discipline: ring-down must use the **Hilbert analytic envelope + normalized MAD threshold** (median|env−med|×1.4826). Using |band-passed| as the envelope measures the oscillation half-period (empirical lesson: a 0.078 ms quantization artifact); a wrong MAD formula pushes the threshold so high that almost nothing crosses it

---

## 2. T: attractor-failure / state-detachment gatekeeper (secondary operator)

### 2.1 Definition

T is a **meta-operator / gatekeeper** measuring "unstable directional readings"; its output is a negative alarm (red light) by design. Do not judge it by "does it produce a signal" — **the standard is "does it hold the line when failure happens" (fail-loud)**. Design key: the state marker must use spatial morphology *orthogonal* to directional jumps (e.g. a polar vortex tearing from a single core into two), or the definition becomes circular.

### 2.2 Cross-domain evidence: red-light cause = state detachment

| Domain | Event | Evidence |
|---|---|---|
| Typhoon 2D tracks | recurvature/redirection | blind state 56.3°/38% vs steady 20.0°/89% sharp cliff; catastrophic cell 101.4° |
| ECG (MIT-BIH) | ectopic beats | state detachment ρ=**+0.55** dominant (energy residual only 0.17); normal-beat axis jump 1.41° vs ventricular 25.07° |
| Polar vortex | wave-2 splitting | split state ρ=**+0.402**; controlling for energy, the energy term vanishes ρ=**+0.020** |

Unified physical cause: **"the system detaches from its previous stable attractor state (state detachment / topological rearrangement)"** — the energy hypothesis is downgraded, the noise hypothesis eliminated. The measurement generalizes; the cause must be specified per domain (typhoon = recurvature, ECG = ectopic beats).

### 2.3 Semantic and implementation boundaries (important)

- T semantics are **bound to the 2D trajectory domain** (directional readings / blind state / redirection are 2D concepts); 1D waveguide domains (MJO-like) are upstream background modulators (they modulate genesis 2.5× — confirmed — but not recurvature — double null), not parallel cross-domain targets. Mismatched dimension semantics is a category error
- **Implementation is domain-adapted:** geometric thresholds (45° etc.) transplanted across fields failed 4/4 in practice; every new domain requires re-derivation + sensitivity checks
- Current status: **candidate retained** (semantics generalize across three domains; implementation not finalized). Using T on 1D scalar domains requires re-deriving "state detachment" — a forward-looking interface, not a ready-made plug

---

## 3. Hard mathematical and engineering boundaries (fail-loud guardrails)

1. **610-point (10-season) reliable floor**; 180 points is coarse tier only; **a 20-point short window is fatal** (τ underestimated ~2×; short windows miss the long tail)
2. **Geometric thresholds must not be transplanted across fields** — granularity/thresholds must be re-derived per target domain's dynamic time scale and calibrated on a ground-truth domain first
3. **Only effective total damping is measured; internal physical components are not decomposed**
4. **Sign interpretation must name the axis** (dual-axis damage hypothesis)
5. **Caliber discipline:** envelope = Hilbert analytic envelope, threshold = normalized MAD; code/docs/templates unified to the code as single truth
6. **Preregistration discipline:** lock criteria before looking at data; amendments may fix definition bugs only, never threshold criteria; exploratory results are not validation

---

## 4. Engineering and CFD single-point probe interface table

**Wiring rule in one sentence:** one scalar series in a stationary window → IDV yields γ + robustness (convergence + null verified); if you have directional readings → add T for red lights. Evidence status is marked per row:

| Engineering / simulation domain | Single-channel input x(t) | IDV measure (recovery damping) | T measure (failure / detachment) | Evidence |
|---|---|---|---|---|
| Atmosphere / geophysics (empirical baseline) | zonal wind u60(t), etc. | intrinsic damping τ≈12.3 days, γ≈0.081/day | polar-vortex split red light (wave-2 orthogonal morphology) | ✅ measured |
| Medical signals (empirical baseline) | single-lead ECG | — (damping fingerprint is domain-specific; ECG counter-evidence) | ectopic-beat red light ρ=0.55 | ✅ measured |
| Machinery / rotating machines (empirical baseline) | bearing vibration impulse response | ring-down τ tracking damage evolution (axis-dependent: defect-size axis τ↑ / wear axis τ↓) | — | ✅ measured (CWRU / XJTU) |
| Unsteady CFD | wake / wall single-point pressure p(t) or velocity u(t) | aerodynamic damping τ (flutter / Hopf-bifurcation critical softening early warning) | transient-drift red light extinguishing = true convergence criterion (limit-cycle steady state) | 🔭 forward-looking (physically plausible, untested) |
| Edge sensors (MCU / automotive) | single-axis vibration a_z(t) or yaw rate | structural / damper damping evolution tracking (lifetime direction must be axis-calibrated first) | state-detachment early warning (bearing early pitting, tire slip; naturally rejects high-energy-rough-road false alarms — energy term measured to vanish) | 🔭 forward-looking (mechanical-domain τ has empirical grounding; T single-axis interface pending) |
| Neural networks / training | per-sample loss residual or logits-drift scalar | representational residual capacity (weights "hardened" / critical softening) | state jump on the eve of catastrophic forgetting / representation collapse | ❌ untested hypothesis (zero data, pure speculation — not a promise) |

> Edge-device spec: KB RAM, millisecond tests, 610-point floor, pure numpy/scipy — the open implementation is on PyPI: **`pip install dampmeter`** (v1.0.0).

---

## 5. Open-use statement

This note is **open publication + technology diffusion**: method, caliber and boundaries are all in the open. Anyone may use, package and deploy it (CC BY 4.0). The only requests are to **use the corrected calibers** (the Hilbert analytic envelope + normalized MAD estimator) and to ground-truth-calibrate before entering a new domain — not for our sake, but so your own readings do not produce fake numbers.

---

## Evidence-grade checklist

| Conclusion | Grade | Basis |
|------|---------|------|
| IDV three-method convergence measures intrinsic τ (vortex 12.3 d) | ✅ cross-validated | 3-estimator convergence + shuffle null collapse to 0.9 d |
| 610-point floor, KB-RAM / millisecond compute | ✅ cross-validated | 2,800-bootstrap characterization |
| T red-light cause = state detachment (not energy/noise) | ✅ cross-validated | ECG 0.55 + vortex 0.402; energy term → 0.020 when controlled |
| Dual-axis damage (defect axis τ↑ / wear axis τ↓) | ✅ cross-validated | CWRU 11/11 + XJTU 13/15 + artifact exclusion |
| Geometric thresholds not transplantable across fields | ✅ cross-validated | four domains, four failures |
| IDV cross-domain universal τ convergence | 🔄 to verify | convergence on one domain (vortex); bearing/ECG single-caliber only |
| CFD convergence criterion / flutter early warning | 🔭 forward-looking | physically plausible, zero measurement |
| MCU lifetime-prognosis interface | 🔭 forward-looking | mechanical τ grounded; direction needs axis calibration |
| NN training probe | ❌ unverifiable | zero-data hypothesis |

## Open Questions

1. Can IDV τ convergence be fully replicated on a second physical domain (candidate: Greenland ice-stream velocity, confirmatory)?
2. How should T's "state detachment" be re-defined for 1D scalar domains (edge sensors, NN training dynamics)?
3. Can the two axis mechanisms of dual-axis damage (defect softening vs debris/lubricant new damping sources) be separated by independent observation?
4. Can the CFD limit-cycle convergence criterion (transient-drift red light extinguishing) genuinely replace residual-based convergence criteria?

## References

- Open implementation: `dampmeter` 1.0.0, PyPI — https://pypi.org/project/dampmeter/
- Method details, usage boundaries and the full engineering interface table: the Zenodo deposit accompanying this note (includes all source code)

---

*AI contribution disclosure: MKP provided original intuition and research direction; tygtDc, Deep Research performed structural synthesis and writing; Gemini served as dialogue partner on writing strategy.*
