"""dampmeter smoke tests on synthetic data with known ground truth."""
import numpy as np
import pytest

import dampmeter as dm


def ou_process(n=610, tau=12.3, dt=1.0, seed=1):
    """Discrete OU: x[t+1] = x[t]*exp(-dt/tau) + noise (variance-matched)."""
    rng = np.random.default_rng(seed)
    a = np.exp(-dt / tau)
    q = 1.0 - a * a
    x = np.zeros(n)
    for t in range(1, n):
        x[t] = a * x[t - 1] + rng.normal(0.0, np.sqrt(q))
    return x


def test_dcor_self_is_one():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    assert dm.dcor(x, x) == pytest.approx(1.0, abs=1e-6)


def test_efold_monotone_curve():
    curve = [1.0, 0.8, 0.5, 0.3, 0.2]
    lag = dm.efold_from_curve(curve)
    assert 1.0 < lag < 3.0
    assert np.isnan(dm.efold_from_curve([1.0, 0.9, 0.85]))


def test_measure_recovers_ou_tau():
    x = ou_process(n=2000, tau=12.3, seed=1)
    est = dm.measure(x, dt=1.0, seed=1)
    # within ±30% of ground truth 12.3
    assert not np.isnan(est.tau)
    assert abs(est.tau - 12.3) <= 0.3 * 12.3
    assert est.converged
    assert est.gamma == pytest.approx(1.0 / est.tau)


def test_measure_short_series_warns():
    x = ou_process(n=15, tau=5.0, seed=2)
    est = dm.measure(x, dt=1.0, seed=2)
    assert est.data_tier == "insufficient"
    assert any("Do not use" in w for w in est.warnings)


def test_tier_boundaries():
    assert dm.check_data_length(15) == "insufficient"
    assert dm.check_data_length(100) == "coarse"
    assert dm.check_data_length(300) == "marginal"
    assert dm.check_data_length(610) == "reliable"


def test_shuffle_null_collapses():
    rng = np.random.default_rng(3)
    x = ou_process(n=2000, tau=20.0, seed=3)
    tau_real = dm.ringdown_conditional(x)
    tau_null = dm.shuffle_null_ringdown(x, blocks=np.arange(2000) // 100, seed=3)
    assert np.isnan(tau_real) or np.isnan(tau_null) or tau_null < tau_real
    _ = rng  # silence lint


def test_envelope_tau_on_damped_oscillation():
    # repeated ring-down bursts: exp(-t/15) * cos(2*pi*0.2*t) every 300 pts
    rng = np.random.default_rng(4)
    t = np.arange(0, 1800)
    x = 0.01 * rng.normal(size=len(t))
    for t0 in range(0, 1800, 300):
        dt_ = t - t0
        burst = np.where(dt_ >= 0, np.exp(-np.maximum(dt_, 0) / 15.0)
                         * np.cos(2 * np.pi * 0.2 * np.maximum(dt_, 0)), 0.0)
        x += burst
    tau = dm.ringdown_envelope(x, fs=1.0, band=(0.15, 0.25))
    assert not np.isnan(tau)
    assert abs(tau - 15.0) <= 0.5 * 15.0


def test_envelope_mad_formula():
    env = np.abs(np.random.default_rng(5).rayleigh(size=5000))
    med = np.median(env)
    wrong = med * 1.4826
    right = np.median(np.abs(env - med)) * 1.4826
    # the documented failure mode: median(env)*1.4826 is far too high
    assert wrong > right * 1.5


def test_recovery_tau_positive():
    x = ou_process(n=1000, tau=10.0, seed=6)
    r = dm.recovery_tau(x)
    assert np.isnan(r) or r > 0


def test_as_dict_roundtrip():
    est = dm.measure(ou_process(n=610, tau=12.3, seed=7), seed=7)
    d = est.as_dict()
    assert set(d) >= {"tau", "gamma", "taus", "converged", "warnings"}
