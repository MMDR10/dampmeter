"""Null controls.

A damping estimate is only physical if it collapses when temporal
structure is destroyed. Two standard controls:

* :func:`shuffle_null_ringdown` — shuffle within blocks (e.g. within
  each year/season) to destroy temporal structure while preserving the
  marginal distribution; the ring-down tau should collapse.
* :func:`phase_surrogate_taus` — phase-randomized surrogates preserve
  the power spectrum; ACF tau should survive (same spectrum) while the
  point is to check estimator stability, not to falsify physics.

If the null does NOT collapse, you are reading a mathematical necessity
(e.g. a monotone trend or a mean offset), not physics — invalidate.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from ._core import acf_efold
from ._ringdown import ringdown_conditional


def shuffle_null_ringdown(
    x: Sequence[float],
    blocks: Optional[Sequence[int]] = None,
    threshold: float = 2.0,
    seed: int = 0,
) -> float:
    """Ring-down tau after within-block shuffling (NaN if no events).

    ``blocks``: label per sample (e.g. year numbers); shuffling happens
    within each block. Default: one block (whole series).
    """
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    xs = x.copy()
    blocks = np.asarray(blocks) if blocks is not None else np.zeros(len(x), int)
    for b in np.unique(blocks):
        m = np.where(blocks == b)[0]
        xs[m] = rng.permutation(xs[m])
    dx = np.diff(xs, prepend=xs[0])
    sig = float(np.std(dx[1:])) if len(x) > 2 else 0.0
    pulse = np.abs(dx) > threshold * sig
    pulse[0] = False
    return ringdown_conditional(xs, pulse_mask=pulse)


def phase_surrogate_taus(
    x: Sequence[float],
    n: int = 10,
    seed: int = 0,
) -> np.ndarray:
    """ACF tau of ``n`` phase-randomized surrogates (spectrum-preserving)."""
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        X = np.fft.rfft(x - x.mean())
        ph = rng.uniform(0, 2 * np.pi, len(X))
        ph[0] = 0
        xs = np.fft.irfft(np.abs(X) * np.exp(1j * ph), n=len(x))
        out.append(acf_efold(xs))
    return np.asarray(out)
