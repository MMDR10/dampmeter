"""Ring-down estimators.

Two calibers, deliberately kept separate:

* :func:`ringdown_conditional` — impulse-conditioned decay R(k), the
  estimator validated in the polar-vortex IDV case (tau3).
* :func:`ringdown_envelope` — analytic-envelope decay on band-passed
  oscillatory signals (bearing / vibration style). This is the
  CORRECTED caliber: Hilbert analytic envelope + normalized MAD event
  threshold ``median|env - median(env)| * 1.4826``.

Known failure modes these implementations are built against:

* using |band-passed| as the "envelope" quantizes tau at the oscillation
  half-period (measures ringing, not damping);
* using ``median(env) * 1.4826`` as the MAD is wrong — the threshold
  lands near ~5 sigma under a Rayleigh envelope and almost nothing
  crosses it.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy import signal as sp_signal


def ringdown_conditional(
    x: Sequence[float],
    pulse_mask: Optional[Sequence[bool]] = None,
    threshold: float = 2.0,
    max_lag: Optional[int] = None,
    min_events: int = 10,
) -> float:
    """tau3: impulse-conditioned ring-down e-fold (in samples).

    R(k) = E[x(t)x(t+k) | pulse] / E[x(t)^2 | pulse]: decay of the
    conditional response after large-impulse days. NaN when fewer than
    ``min_events`` usable pulses.
    """
    x = np.asarray(x, float)
    if pulse_mask is None:
        dx = np.diff(x, prepend=x[0])
        sig = float(np.std(dx[1:])) if len(x) > 2 else 0.0
        pulse_mask = np.abs(dx) > threshold * sig
        pulse_mask[0] = False
    pulse_mask = np.asarray(pulse_mask, bool)
    max_lag = max_lag or max(10, min(len(x) // 5, 100))
    ps = np.where(pulse_mask)[0]
    if len(ps) < min_events:
        return float("nan")
    base = float(np.mean(x[pulse_mask] ** 2))
    if base <= 0:
        return float("nan")
    curve = [1.0]
    for k in range(1, max_lag + 1):
        pk = ps[ps + k < len(x)]
        if len(pk) < min_events:
            curve.append(float("nan"))
            continue
        curve.append(float(np.mean(x[pk] * x[pk + k]) / base))
    from ._core import efold_from_curve

    return efold_from_curve(curve)


def _normalized_mad(env: np.ndarray) -> float:
    """CORRECT MAD: median|env - median(env)| * 1.4826."""
    med = float(np.median(env))
    return float(np.median(np.abs(env - med)) * 1.4826)


def ringdown_envelope(
    x: Sequence[float],
    fs: float = 1.0,
    band: Optional[Tuple[float, float]] = None,
    mad_k: float = 1.0,
    min_events: int = 5,
) -> float:
    """Analytic-envelope ring-down decay time (in time units = 1/fs).

    Band-pass (optional), take the Hilbert analytic envelope, smooth it,
    threshold events at ``mad_k * normalized_MAD``, then fit the
    exponential decay of each event envelope tail (peak → thr/e); the
    median event tau is returned. NaN when fewer than ``min_events`` events.

    ``band``: (low, high) in Hz, relative to fs. Default = no band-pass
    (use when the series is already band-limited).

    The one caliber rule: envelope = |hilbert(bandpassed)|. Never
    |bandpassed| — that quantizes tau to the oscillation half-period.
    """
    x = np.asarray(x, float)
    if len(x) < 20:
        return float("nan")
    y = x - x.mean()
    if band is not None:
        lo, hi = band
        nyq = fs / 2.0
        b, a = sp_signal.butter(2, [lo / nyq, hi / nyq], btype="band")
        y = sp_signal.filtfilt(b, a, y)
    env = np.abs(sp_signal.hilbert(y))
    # light smoothing of the envelope (3-point moving average)
    if len(env) >= 3:
        env = np.convolve(env, np.ones(3) / 3.0, mode="same")

    # trim Hilbert/convolution edge artifacts
    edge = min(5, len(env) // 10)
    env = env[edge:len(env) - edge] if edge > 0 else env

    thr = mad_k * _normalized_mad(env)
    if thr <= 0:
        return float("nan")

    # event onsets: upward crossings of the threshold
    above = env > thr
    onsets = np.where(np.diff(above.astype(int)) == 1)[0] + 1
    taus: List[float] = []
    for i, t0 in enumerate(onsets):
        # tail: from local peak to thr/e, cut before the next event rises
        seg = env[t0:]
        pk_rel = int(np.argmax(seg[:min(10, len(seg))]))
        tail = seg[pk_rel:]
        stop = np.where(tail <= thr / np.e)[0]
        end = int(stop[0]) if len(stop) else len(tail)
        if i + 1 < len(onsets):
            end = min(end, int(onsets[i + 1] - t0 - pk_rel))
        if end < 4:
            continue
        tail = tail[:end]
        yy = np.log(np.maximum(tail, 1e-300))
        tt = np.arange(len(tail), dtype=float)
        slope = np.polyfit(tt, yy, 1)[0]
        if slope < 0:
            taus.append(-1.0 / slope)
    if len(taus) < min_events:
        return float("nan")
    return float(np.median(taus))
