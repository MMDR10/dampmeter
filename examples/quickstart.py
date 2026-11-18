"""Quickstart: measure the intrinsic damping time of a scalar series."""
import numpy as np

import dampmeter as dm


def main():
    # synthetic Ornstein-Uhlenbeck process with known tau = 12.3
    rng = np.random.default_rng(20261118)
    n, tau_true = 2000, 12.3
    a = np.exp(-1.0 / tau_true)
    x = np.zeros(n)
    for t in range(1, n):
        x[t] = a * x[t - 1] + rng.normal(0.0, np.sqrt(1 - a * a))

    est = dm.measure(x, dt=1.0)
    print("ground truth tau =", tau_true)
    print("estimate         =", round(est.tau, 2))
    print("gamma = 1/tau    =", round(est.gamma, 5))
    print("per-method taus  :", {k: round(v, 2) for k, v in est.taus.items()})
    print("converged        :", est.converged, f"({est.within_30pct}/{est.n_methods} within ±30%)")
    print("data tier        :", est.data_tier)

    # null control: destroy temporal structure within 100-sample blocks
    tau_null = dm.shuffle_null_ringdown(x, blocks=np.arange(n) // 100)
    print("shuffle-null ring-down tau =", round(tau_null, 2), "(should collapse)")


if __name__ == "__main__":
    main()
