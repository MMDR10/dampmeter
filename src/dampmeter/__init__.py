"""dampmeter — IDV individual damping caliper.

Measure a system's intrinsic relaxation (damping) time constant ``tau``
from a single scalar time series, with four independent estimators,
convergence-as-identity-check, and a shuffle null control.

Terminology: this package measures the *effective damping* (dissipative
time constant tau, gamma = 1/tau) — NOT the potential-well stiffness
(curvature, omega0^2). The two are different physical quantities.

The output is an *effective* gamma only: it says "how much, and how
stable", never "which mechanism". Sign interpretation of damping change
is axis-dependent (defect-size axis vs wear/life axis can go opposite
directions) — always specify the axis before interpreting trends.

Caliber note: the ring-down envelope estimator uses the Hilbert
analytic envelope with the normalized MAD threshold
median|env - median(env)| * 1.4826. Using |band-passed| as the envelope
measures the oscillation half-period instead of damping decay — do not
swap the calibers.
"""
from ._core import (
    EFOLD,
    HARD_FLOOR_POINTS,
    COARSE_MIN_POINTS,
    RECOMMENDED_MIN_POINTS,
    DampingEstimate,
    acf_curve,
    acf_efold,
    check_data_length,
    dcor,
    dcor_efold,
    efold_from_curve,
    measure,
    recovery_tau,
)
from ._ringdown import ringdown_conditional, ringdown_envelope
from ._null import phase_surrogate_taus, shuffle_null_ringdown

__version__ = "1.0.0"

__all__ = [
    "EFOLD",
    "HARD_FLOOR_POINTS",
    "COARSE_MIN_POINTS",
    "RECOMMENDED_MIN_POINTS",
    "DampingEstimate",
    "acf_curve",
    "acf_efold",
    "check_data_length",
    "dcor",
    "dcor_efold",
    "efold_from_curve",
    "measure",
    "recovery_tau",
    "ringdown_conditional",
    "ringdown_envelope",
    "phase_surrogate_taus",
    "shuffle_null_ringdown",
    "__version__",
]
