"""Core estimators: ACF e-fold, distance-correlation e-fold, recovery,
convergence check, and the top-level ``measure`` entry point."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from ._ringdown import ringdown_conditional

#: 1/e crossing point for e-fold decay curves.
EFOLD = float(1.0 / np.e)

# Minimum-data reliability tiers (2,800-bootstrap characterization,
# ±30% hit rate): 20 pts -> 13% (tau underestimated ~2x, hard floor),
# ~180 pts -> 57% (coarse tier), ~610 pts -> 81.5% (reliable tier).
HARD_FLOOR_POINTS = 20
COARSE_MIN_POINTS = 180
RECOMMENDED_MIN_POINTS = 610


def check_data_length(n: int) -> str:
    """Return the reliability tier for a series of ``n`` samples.

    ``"insufficient"`` (<20: do not use), ``"coarse"`` (20..179),
    ``"reliable"`` (>=610), ``"marginal"`` (180..609).
    """
    if n < HARD_FLOOR_POINTS:
        return "insufficient"
    if n < COARSE_MIN_POINTS:
        return "coarse"
    if n < RECOMMENDED_MIN_POINTS:
        return "marginal"
    return "reliable"


def dcor(x: Sequence[float], y: Sequence[float]) -> float:
    """Distance correlation between two 1-D samples (biased estimator,
    O(n^2) memory — subsample upstream for long series)."""
    x = np.asarray(x, float)[:, None]
    y = np.asarray(y, float)[:, None]
    a = np.abs(x - x.T)
    b = np.abs(y - y.T)
    A = a - a.mean(0)[None, :] - a.mean(1)[:, None] + a.mean()
    B = b - b.mean(0)[None, :] - b.mean(1)[:, None] + b.mean()
    dcov2 = (A * B).mean()
    dvx = (A * A).mean()
    dvy = (B * B).mean()
    if dvx <= 0 or dvy <= 0:
        return 0.0
    return float(np.sqrt(max(dcov2, 0.0) / np.sqrt(dvx * dvy)))


def efold_from_curve(curve: Sequence[float]) -> float:
    """Lag (float) where a normalized decay curve first crosses 1/e.

    ``curve[0]`` must be 1.0; NaN entries are tolerated. Linear
    interpolation between bracketing points. Returns NaN when the curve
    never reaches 1/e (e.g. a plateau lifted by a mean offset — check
    your detrending).
    """
    for k in range(1, len(curve)):
        c1 = curve[k]
        if np.isnan(c1):
            continue
        if c1 < EFOLD:
            c0 = curve[k - 1]
            if np.isnan(c0):
                return float(k)
            if c0 == c1:
                return float(k)
            f = (c0 - EFOLD) / (c0 - c1)
            return float(k - 1 + f)
    return float("nan")


def acf_curve(x: Sequence[float], max_lag: int) -> List[float]:
    """Normalized autocorrelation curve, lags 0..max_lag (curve[0] = 1)."""
    x = np.asarray(x, float)
    x = x - x.mean()
    var = float(np.dot(x, x))
    if var <= 0:
        return [1.0] + [float("nan")] * max_lag
    curve = [1.0]
    for k in range(1, max_lag + 1):
        if k >= len(x):
            curve.append(float("nan"))
            continue
        curve.append(float(np.dot(x[:-k], x[k:]) / var))
    return curve


def acf_efold(x: Sequence[float], max_lag: Optional[int] = None) -> float:
    """tau1: ACF e-fold decay time (in samples)."""
    n = len(x)
    max_lag = max_lag or max(10, min(n // 3, 200))
    return efold_from_curve(acf_curve(x, max_lag))


def dcor_efold(
    x: Sequence[float],
    max_lag: Optional[int] = None,
    subsample: int = 1200,
    seed: int = 0,
) -> float:
    """tau2: distance-correlation e-fold decay time (in samples).

    Nonlinear counterpart of the ACF; shift-invariant, useful for
    cross-checking. Pair indices are subsampled beyond ``subsample``
    to keep the O(n^2) kernel tractable.
    """
    x = np.asarray(x, float)
    n = len(x)
    max_lag = max_lag or max(10, min(n // 3, 200))
    rng = np.random.default_rng(seed)
    curve = [1.0]
    for k in range(1, max_lag + 1):
        idx = np.arange(0, n - k)
        if len(idx) < 30:
            curve.append(float("nan"))
            continue
        if len(idx) > subsample:
            idx = rng.choice(idx, subsample, replace=False)
        curve.append(dcor(x[idx], x[idx + k]))
    return efold_from_curve(curve)


def recovery_tau(
    x: Sequence[float],
    drop_mask: Optional[Sequence[bool]] = None,
    threshold: float = 2.0,
    cap: Optional[int] = None,
) -> float:
    """tau4 (auxiliary): median time-to-recover after large drops.

    Right-censoring warning: when the window ends before recovery, the
    capped value is *censoring*, not an estimate. Never treat tau4 alone
    as evidence; use it as a supporting check only.
    """
    x = np.asarray(x, float)
    if drop_mask is None:
        dx = np.diff(x, prepend=x[0])
        sig = float(np.std(dx[1:])) if len(x) > 2 else 0.0
        drop_mask = dx < -threshold * sig
        drop_mask[0] = False
    drop_mask = np.asarray(drop_mask, bool)
    cap = cap or max(10, len(x) // 20)
    recs = []
    for t in np.where(drop_mask)[0]:
        if t < 1 or t + cap >= len(x):
            continue
        lvl = x[t - 1]
        r = cap
        for k in range(1, cap + 1):
            if x[t + k] >= lvl:
                r = k
                break
        recs.append(r)
    return float(np.median(recs)) if recs else float("nan")


@dataclass
class DampingEstimate:
    """Result of :func:`measure`.

    ``tau`` is the converged estimate in *time units* (samples x dt);
    ``gamma = 1/tau``. ``converged`` is True when at least
    ``min_converging`` estimators fall within ±30% of the median —
    convergence is an identity check ("this constant is really there"),
    not a goodness-of-fit score.
    """

    tau: float
    gamma: float
    median_tau: float
    taus: Dict[str, float]
    converged: bool
    within_30pct: int
    n_methods: int
    n_samples: int
    data_tier: str
    dt: float
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return {
            "tau": self.tau,
            "gamma": self.gamma,
            "median_tau": self.median_tau,
            "taus": self.taus,
            "converged": self.converged,
            "within_30pct": f"{self.within_30pct}/{self.n_methods}",
            "n_samples": self.n_samples,
            "data_tier": self.data_tier,
            "dt": self.dt,
            "warnings": self.warnings,
        }


def measure(
    x: Sequence[float],
    dt: float = 1.0,
    max_lag: Optional[int] = None,
    pulse_threshold: float = 2.0,
    min_converging: int = 3,
    dcor_subsample: int = 1200,
    seed: int = 0,
    include_recovery: bool = True,
) -> DampingEstimate:
    """Estimate the intrinsic damping time tau of a scalar time series.

    Runs three independent estimators — tau1 ACF e-fold, tau2 dCor
    e-fold, tau3 impulse-conditioned ring-down — plus the auxiliary tau4
    post-drop recovery. Convergence: >= ``min_converging`` methods within
    ±30% of the median => the constant is considered established.

    Parameters
    ----------
    x : scalar series (assumed stationary window, effectively linear
        relaxation). Detrend before calling — fitting on raw data with a
        mean offset lifts the ACF plateau so the curve never reaches 1/e.
    dt : sample spacing; tau and gamma are returned in these time units.
    """
    x = np.asarray(x, float)
    n = len(x)
    warnings: List[str] = []
    tier = check_data_length(n)
    if tier == "insufficient":
        warnings.append(
            f"n={n} < {HARD_FLOOR_POINTS}: below the hard floor — tau is "
            "systematically underestimated (~2x on short windows). Do not use."
        )
    elif tier == "coarse":
        warnings.append(
            f"n={n}: coarse tier (~57% hit rate at 180 pts). Treat as "
            "indicative only; >=610 pts recommended for reliable tier."
        )
    elif tier == "marginal":
        warnings.append(f"n={n}: between coarse and reliable tiers (~610 pts).")

    tau1 = acf_efold(x, max_lag=max_lag)
    tau2 = dcor_efold(x, max_lag=max_lag, subsample=dcor_subsample, seed=seed)

    dx = np.diff(x, prepend=x[0])
    sig = float(np.std(dx[1:])) if n > 2 else 0.0
    pulse = np.abs(dx) > pulse_threshold * sig
    pulse[0] = False
    tau3 = ringdown_conditional(x, pulse_mask=pulse, max_lag=max_lag)

    taus = {"tau1_acf": tau1, "tau2_dcor": tau2, "tau3_ringdown": tau3}
    if include_recovery:
        taus["tau4_recovery"] = recovery_tau(x, drop_mask=dx < -pulse_threshold * sig)

    valid = {k: v for k, v in taus.items() if not np.isnan(v)}
    if not valid:
        return DampingEstimate(
            tau=float("nan"), gamma=float("nan"), median_tau=float("nan"),
            taus=taus, converged=False, within_30pct=0,
            n_methods=len(taus), n_samples=n, data_tier=tier, dt=dt,
            warnings=warnings + ["no estimator produced a finite tau"],
        )

    # tau4 is auxiliary (right-censoring prone): convergence is judged
    # on tau1/tau2/tau3, matching the validated IDV protocol.
    primary = {k: v for k, v in valid.items() if k != "tau4_recovery"} or valid
    med = float(np.median(list(primary.values())))
    within = sum(1 for v in primary.values() if abs(v - med) <= 0.3 * med)
    converged = within >= min_converging
    if not converged:
        warnings.append(
            f"only {within}/{len(primary)} primary estimators within ±30% of "
            f"median {med:.3g}: individual damping NOT established (either the "
            "window is not a stationary relaxation regime, or data too short)."
        )
    tau = med * dt
    return DampingEstimate(
        tau=tau,
        gamma=(1.0 / tau) if tau > 0 else float("nan"),
        median_tau=med,
        taus={k: (v * dt if not np.isnan(v) else float("nan")) for k, v in taus.items()},
        converged=converged,
        within_30pct=within,
        n_methods=len(primary),
        n_samples=n,
        data_tier=tier,
        dt=dt,
        warnings=warnings,
    )
