#!/usr/bin/env python3
"""openfoam_probes_example.py — feed OpenFOAM postProcessing output into dampmeter.

Zero C++ compilation, zero solver changes: this script only reads files that
OpenFOAM already writes (probes / forceCoeffs / solver logs) and measures the
effective damping time tau (gamma = 1/tau) with ``dampmeter``.

Usage
-----
  python openfoam_probes_example.py probes    postProcessing/probes/0/p [probe_index]
  python openfoam_probes_example.py envelope  postProcessing/forceCoeffs/0/coefficient.dat [column_index]
  python openfoam_probes_example.py residuals log.pimpleFoam [field]
  python openfoam_probes_example.py windows   postProcessing/probes/0/p [probe_index] [window] [stride]
  python openfoam_probes_example.py grid      coarse/p fine/p
  python openfoam_probes_example.py demo

Modes
-----
  probes     monitor-point series (p, T, scalar)  -> dm.measure()  -> tau
             vector fields (U) are reduced to |U| per probe.
  envelope   forceCoeffs time series (Cl of vortex shedding / flutter)
             -> ringdown_envelope() -> amplitude decay time.
  residuals  solver log initial residuals -> log10 -> dm.measure()
             -> "is this convergence real?" check.
  windows    sliding-window tau (default 610 @ stride 305) -> tau stability
             -> the false-convergence fuse: jumping tau = flow not settled.
  grid       same monitor on two meshes / schemes -> tau_eff ratio
             -> quantifies numerical dissipation (grid independence).
  demo       synthetic signals with known tau, full pipeline self-test.

Gotchas (read CFD_GUIDE_ZH.md §4 for the long version)
------------------------------------------------------
  * adjustTimeStep yes -> non-uniform timestamps -> this script resamples
    to uniform dt and tells you so. Check the output.
  * Need >= ~610 samples for the reliable tier; 20 points is a hard floor.
  * Always detrend (this script does a linear detrend for you).
  * tau_eff contains BOTH physical and numerical viscosity -- that is a
    feature in `grid` mode, a trap everywhere else.

Evidence status: validated on atmospheric vortex / ECG / bearing vibration.
The CFD application is forward-looking and has NOT been validated on CFD
data yet. Run it, break it, tell us where.
"""
from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path

import numpy as np

import dampmeter as dm

HARD_FLOOR = 20
RELIABLE_TIER = 610


# ---------------------------------------------------------------- helpers

def detrend_linear(x: np.ndarray) -> np.ndarray:
    """Remove a least-squares linear trend (numpy only)."""
    t = np.arange(len(x), dtype=float)
    k, b = np.polyfit(t, x, 1)
    return x - (k * t + b)


def load_probes(path: str, iprobe: int = 0):
    """Parse an OpenFOAM ``probes`` output file.

    Returns (t, x, nprobes). Handles scalar fields (p, T ...) with one
    value per probe and vector fields (U) with one ``(u v w)`` triple per
    probe (reduced to magnitude).
    """
    t_all, vals_all = [], []
    nprobes = 0
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = re.search(r"Probe\s+(\d+)", s)
                if m:
                    nprobes = max(nprobes, int(m.group(1)) + 1)
                continue
            tok = s.replace("(", " ").replace(")", " ").split()
            t_all.append(float(tok[0]))
            vals_all.append([float(v) for v in tok[1:]])
    t = np.asarray(t_all)
    if not vals_all:
        raise SystemExit(f"no data rows in {path}")
    width = len(vals_all[0])
    if nprobes == 0:  # header missing -> infer
        nprobes = width if width % 3 else width // 3
    if iprobe >= nprobes:
        raise SystemExit(f"probe index {iprobe} out of range (file has {nprobes})")
    if width == nprobes:  # scalar field
        x = np.array([row[iprobe] for row in vals_all])
        kind = "scalar"
    elif width == 3 * nprobes:  # vector field -> magnitude
        x = np.array([
            float(np.linalg.norm(row[3 * iprobe:3 * iprobe + 3]))
            for row in vals_all
        ])
        kind = "vector->|U|"
    else:
        raise SystemExit(
            f"unexpected column count {width} for {nprobes} probe(s) in {path}"
        )
    return t, x, nprobes, kind


def ensure_uniform_dt(t: np.ndarray, x: np.ndarray):
    """Return (x_uniform, dt, resampled: bool). Linear resample if needed."""
    diffs = np.diff(t)
    dt = float(np.median(diffs))
    if dt <= 0:
        raise SystemExit("non-increasing time stamps")
    if np.max(np.abs(diffs - dt)) > 0.01 * dt:
        t_u = np.arange(t[0], t[-1] + dt / 2, dt)
        return np.interp(t_u, t, x), dt, True
    return x, dt, False


def report(x: np.ndarray, dt: float, label: str) -> dm.DampingEstimate | None:
    """Detrend + measure + shuffle-null, printed engineer-style."""
    n = len(x)
    print(f"\n=== {label} ===")
    print(f"samples: {n}  dt: {dt:g}")
    if n < HARD_FLOOR:
        print(f"REFUSED: n={n} < {HARD_FLOOR} hard floor. "
              "Write output more often (writeInterval).")
        return None
    if n < RELIABLE_TIER:
        print(f"WARNING: n={n} < ~{RELIABLE_TIER} reliable tier -- coarse result.")
    xd = detrend_linear(x)
    est = dm.measure(xd, dt=dt)
    print(f"tau_eff = {est.tau:.4g}   gamma = {est.gamma:.4g}")
    print(f"converged: {est.converged}   "
          f"(within +-30%: {est.within_30pct}/{est.n_methods})")
    for name, v in est.taus.items():
        print(f"    {name:<10s} {v:.4g}")
    for w in est.warnings:
        print(f"    ! {w}")

    blocks = np.arange(n) // max(1, n // 10)  # 10 blocks
    # NOTE: shuffle_null_ringdown returns tau in SAMPLES -> x dt for time units
    null = dm.shuffle_null_ringdown(xd, blocks=blocks, seed=0)
    null_t = null * dt
    print(f"shuffle-null tau: {null_t:.4g} (time units)")
    if np.isnan(null):
        print("  null = NaN: no pulse events survived shuffling -- "
              "no clear event structure to falsify.")
    elif est.tau <= 1.5 * dt:
        print("  NO STRUCTURE: tau ~ 1 sample -- the detrended series is "
              "essentially white. Nothing to falsify, nothing to measure.")
    elif est.tau > 0 and null_t < 0.5 * est.tau:
        print("  PASS: tau collapses after shuffling -- "
              "the time structure you measured is real.")
    else:
        print("  FAIL: tau survives shuffling -- "
              "you may be reading a trend/offset, not physics. Do not trust.")
    return est


# ---------------------------------------------------------------- modes

def mode_probes(path: str, iprobe: int) -> None:
    t, x, nprobes, kind = load_probes(path, iprobe)
    print(f"file: {path}  ({nprobes} probe(s), field kind: {kind}, "
          f"using probe #{iprobe})")
    x_u, dt, resampled = ensure_uniform_dt(t, x)
    if resampled:
        print("! timestamps were non-uniform (adjustTimeStep?) -- "
              "resampled to uniform dt via linear interpolation.")
    report(x_u, dt, f"probes probe#{iprobe}")


def mode_envelope(path: str, icol: int) -> None:
    data = np.loadtxt(path, comments="#", ndmin=2)
    t = data[:, 0]
    if icol >= data.shape[1]:
        raise SystemExit(f"column {icol} out of range ({data.shape[1]} columns)")
    x = data[:, icol]
    x_u, dt, resampled = ensure_uniform_dt(t, x)
    if resampled:
        print("! timestamps non-uniform -- resampled to uniform dt.")
    n = len(x_u)
    print(f"file: {path}  column {icol}  samples: {n}")
    if n < HARD_FLOOR:
        print(f"REFUSED: n={n} < {HARD_FLOOR} hard floor.")
        return
    xd = detrend_linear(x_u)
    # NOTE: ringdown_envelope returns tau in SAMPLES -> x dt for time units
    tau_samp = dm.ringdown_envelope(xd, fs=1.0 / dt)
    tau_t = tau_samp * dt
    print(f"envelope ring-down tau: {tau_t:.4g} (time units) "
          f"[{tau_samp:.4g} samples]")
    print("  oscillatory caliber -- amplitude decay of the Hilbert "
          "envelope, not the ACF of the raw signal.")


def mode_residuals(path: str, field: str | None) -> None:
    pat = re.compile(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+)")
    series: dict[str, list[float]] = {}
    with open(path) as f:
        for line in f:
            m = pat.search(line)
            if m:
                series.setdefault(m.group(1), []).append(float(m.group(2)))
    if not series:
        raise SystemExit(f"no 'Initial residual' lines found in {path}")
    fields = [field] if field else sorted(series)
    print(f"log: {path}  fields: {', '.join(sorted(series))}")
    for fld in fields:
        if fld not in series:
            raise SystemExit(f"field '{fld}' not in log (have: {sorted(series)})")
        r = np.asarray(series[fld])
        print(f"\n[{fld}] {len(r)} iterations, "
              f"residual {r[0]:.2e} -> {r[-1]:.2e}")
        # residuals span orders of magnitude -> work on log10
        report(np.log10(r), 1.0, f"residuals {fld} (log10, dt=1 iteration)")
    print("\nhow to read this (false-convergence check):")
    print("  * converged=False -> residual series still carries time "
          "structure (plateaus / oscillations) -> do NOT stop the run yet.")
    print("  * NO STRUCTURE -> decay looks clean (what is left is solver "
          "noise). Consistent with genuine convergence.")
    print("  * PASS -> real relaxation structure in the residuals.")
    print("  * FAIL -> suspicious: apparent structure is a trend artifact.")
    print("  CAUTION: this measures time STRUCTURE, not residual LEVEL. A "
          "perfect exponential decay is white after detrend even mid-run -- "
          "always check the residual level and the flow field yourself.")


def mode_windows(path: str, iprobe: int, window: int, stride: int) -> None:
    """Sliding-window tau -- the false-convergence fuse."""
    t, x, nprobes, kind = load_probes(path, iprobe)
    x_u, dt, resampled = ensure_uniform_dt(t, x)
    if resampled:
        print(f"! {path}: non-uniform timestamps -- resampled.")
    n = len(x_u)
    print(f"file: {path}  ({kind}, probe#{iprobe})  samples: {n}  "
          f"window: {window}  stride: {stride}")
    if n < window:
        raise SystemExit(f"series ({n}) shorter than one window ({window})")
    taus: list[float] = []
    for s in range(0, n - window + 1, stride):
        seg = detrend_linear(x_u[s:s + window])
        est = dm.measure(seg, dt=dt)
        taus.append(est.tau)
        print(f"  window [{s}:{s + window}]  tau = {est.tau:.4g}  "
              f"converged = {est.converged}")
    arr = np.array([v for v in taus if np.isfinite(v)])
    if len(arr) >= 2 and abs(arr.mean()) > 0:
        cv = float(np.std(arr) / abs(arr.mean()))
        print(f"\ntau across windows: mean {arr.mean():.4g}, CV {cv:.2%}")
        if cv > 0.5:
            print("  UNSTABLE: tau jumps between windows -- the flow has NOT "
                  "settled into a stationary regime. Do not stop the run.")
        else:
            print("  STABLE: tau consistent across windows -- consistent "
                  "with a settled regime (verify residual level too).")
    else:
        print("\nnot enough finite windows to judge stability.")


def mode_grid(path_a: str, path_b: str) -> None:
    ests = {}
    for tag, path in (("A", path_a), ("B", path_b)):
        t, x, _, kind = load_probes(path, 0)
        x_u, dt, resampled = ensure_uniform_dt(t, x)
        if resampled:
            print(f"! {path}: non-uniform timestamps -- resampled.")
        ests[tag] = report(x_u, dt, f"grid {tag} ({path}, {kind})")
    a, b = ests.get("A"), ests.get("B")
    if a and b and a.tau > 0 and b.tau > 0:
        ratio = a.tau / b.tau
        print("\n=== grid / scheme comparison ===")
        print(f"tau_A / tau_B = {ratio:.3f}")
        print("tau_eff = physical + numerical viscosity. The smaller tau, "
              "the heavier the damping.")
        if ratio > 1.3:
            print("  -> B damps >30% faster: B's mesh/scheme adds visible "
                  "numerical dissipation. Grid independence NOT reached for B.")
        elif ratio < 1 / 1.3:
            print("  -> A damps >30% faster: A's mesh/scheme adds visible "
                  "numerical dissipation. Grid independence NOT reached for A.")
        else:
            print("  -> tau_eff agrees within 30%: consistent dissipation "
                  "across the two setups (good sign for grid independence).")


def mode_demo() -> None:
    rng = np.random.default_rng(42)

    # --- signal 1: Ornstein-Uhlenbeck relaxation, tau_true = 0.15 s ---
    # The canonical "perturbation decays inside a stable regime" series
    # (stationary, ongoing small perturbations -- like a settled flow with
    # background unsteadiness).
    print("=" * 60)
    print("demo 1: OU relaxation  tau_true = 0.15 s  (measure caliber)")
    tau_true, dt, n = 0.15, 0.005, 12000
    a = np.exp(-dt / tau_true)
    x = np.zeros(n)
    for k in range(1, n):
        x[k] = a * x[k - 1] + rng.standard_normal() * np.sqrt(1 - a * a)
    report(x, dt, "demo 1: OU process -> dm.measure()")

    # --- signal 2: intermittent oscillation bursts, decay tau = 0.12 s ---
    # 30 ring-down bursts of a 25 Hz carrier (like periodic vortex-shedding
    # bursts after small kicks).
    print("=" * 60)
    print("demo 2: burst-decay oscillation  tau_true = 0.12 s  "
          "(ringdown_envelope caliber)")
    fs = 200.0
    t = np.arange(0, 30, 1 / fs)
    env = np.zeros_like(t)
    for t0 in np.arange(0.5, 29, 1.0):
        m = t >= t0
        env[m] += np.exp(-(t[m] - t0) / 0.12)
    y = env * np.cos(2 * np.pi * 25 * t) + 0.01 * rng.standard_normal(len(t))
    tau_samp = dm.ringdown_envelope(y, fs=fs)
    print(f"envelope ring-down tau: {tau_samp / fs:.4g} s "
          f"[{tau_samp:.4g} samples]  (true 0.12)")

    # --- integration test: write a fake probes file and re-parse it ---
    print("=" * 60)
    with tempfile.NamedTemporaryFile(
        "w", suffix="_p", delete=False
    ) as f:
        f.write("# Probe 0 (0.1 0.05 0)\n# Probe 1 (0.2 0.05 0)\n# Time\n")
        t_ou = np.arange(n) * dt
        for ti, xi in zip(t_ou, x):
            f.write(f"{ti:.6g}\t{xi:.6g}\t{xi * 0.5:.6g}\n")
        p = f.name
    t2, x2, npb, kind = load_probes(p, 0)
    x_u, dt2, resampled = ensure_uniform_dt(t2, x2)
    print(f"parsed fake probes file: {len(x2)} rows, {npb} probes, "
          f"kind={kind}, resampled={resampled}, dt={dt2:g}")
    report(x_u, dt2, "demo 1 again via fake probes file (probe #0)")
    Path(p).unlink()
    print("\ndemo done.")


# ---------------------------------------------------------------- cli

def main(argv: list[str]) -> None:
    if len(argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    mode = argv[1]
    args = argv[2:]
    if mode == "probes":
        if not args:
            raise SystemExit("usage: probes <file> [probe_index]")
        mode_probes(args[0], int(args[1]) if len(args) > 1 else 0)
    elif mode == "envelope":
        if not args:
            raise SystemExit("usage: envelope <coefficient.dat> [column_index]")
        mode_envelope(args[0], int(args[1]) if len(args) > 1 else 2)
    elif mode == "residuals":
        if not args:
            raise SystemExit("usage: residuals <log> [field]")
        mode_residuals(args[0], args[1] if len(args) > 1 else None)
    elif mode == "windows":
        if not args:
            raise SystemExit(
                "usage: windows <file> [probe_index] [window] [stride]")
        mode_windows(
            args[0],
            int(args[1]) if len(args) > 1 else 0,
            int(args[2]) if len(args) > 2 else 610,
            int(args[3]) if len(args) > 3 else 305,
        )
    elif mode == "grid":
        if len(args) < 2:
            raise SystemExit("usage: grid <fileA> <fileB>")
        mode_grid(args[0], args[1])
    elif mode == "demo":
        mode_demo()
    else:
        print(__doc__)
        raise SystemExit(2)


if __name__ == "__main__":
    main(sys.argv)
