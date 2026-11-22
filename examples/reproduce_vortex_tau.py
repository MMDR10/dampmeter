#!/usr/bin/env python3
"""Reproduce the polar-vortex validation numbers from the raw CSV in this repo.

Data: data/vortex_daily_hist.csv — 2,867 daily values of 60°N zonal-mean
zonal wind (u60) at 10 hPa, NCEP/NCAR Reanalysis 1 (PSL THREDDS OPeNDAP),
1979-10-01 .. 2025-11-30.

Published result (dual-operator note): tau1=11.28 d / tau2=10.81 d /
tau3=13.32 d (median 12.3 d), shuffle-null collapses to 0.9 d.

Run:  python examples/reproduce_vortex_tau.py   (needs numpy + pandas only)

CALIBER NOTE (important):
The CSV is 46 ragged Oct-Nov segments, NOT a continuous daily series.
The paper's caliber therefore pairs samples *within the same season*
(season_pairs below, vendored from the preregistered script
individual_damping_test.py). Naively passing the concatenated array to a
generic ACF/measure() assumes uniform 1-day spacing across year gaps and
gives a biased ~9 d. The ~12.3 vs ~9 d gap is a caliber artifact, not a
data dispute — we flag it here on purpose.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "data", "vortex_daily_hist.csv")

SEED = 20261117
EFOLD = 1.0 / np.e


# ---- vendored from the preregistered script (exact caliber) ----------------

def dcor(x, y):
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


def efold_from_curve(curve):
    for k in range(1, len(curve)):
        if curve[k] < EFOLD:
            c0, c1 = curve[k - 1], curve[k]
            if c0 == c1:
                return float(k)
            return float(k - 1 + (c0 - EFOLD) / (c0 - c1))
    return float("nan")


def season_pairs(dates, k):
    """(i, i+k) index pairs where both samples fall in the same season year."""
    sid = dates // 10000
    ok = np.array([sid[i] == sid[i + k] for i in range(len(dates) - k)])
    return np.where(ok)[0]


def acf_efold(x, dates, kmax=30):
    var = np.var(x)
    curve = [1.0]
    for k in range(1, kmax + 1):
        ii = season_pairs(dates, k)
        if len(ii) < 30:
            curve.append(np.nan)
            continue
        curve.append(float(np.mean(x[ii] * x[ii + k]) / var))
    return efold_from_curve(curve)


def dcor_efold(x, dates, kmax=30, subsample=1200, rng=None):
    rng = rng or np.random.default_rng(SEED)
    curve = [1.0]
    for k in range(1, kmax + 1):
        ii = season_pairs(dates, k)
        if len(ii) < 30:
            curve.append(np.nan)
            continue
        if len(ii) > subsample:
            ii = rng.choice(ii, subsample, replace=False)
        curve.append(dcor(x[ii], x[ii + k]))
    return efold_from_curve(curve)


def ringdown_efold(x, dates, pulse_mask, kmax=20):
    base = np.mean(x[pulse_mask] ** 2)
    curve = [1.0]
    for k in range(1, kmax + 1):
        ps = np.where(pulse_mask)[0]
        ps = ps[ps + k < len(x)]
        sp = season_pairs(dates, k)
        ps = np.array([t for t in ps if sp.size and
                       (dates[t] // 10000) == (dates[t + k] // 10000)])
        if len(ps) < 10:
            curve.append(np.nan)
            continue
        curve.append(float(np.mean(x[ps] * x[ps + k]) / base))
    return efold_from_curve(curve)


def recovery_days(x, drop_idx, cap=30):
    recs = []
    for t in drop_idx:
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


# ---- load + caliber (identical to the manuscript figures) ------------------

df = pd.read_csv(CSV)
df = df[df["level"] == 10] if "level" in df.columns else df
df["date"] = pd.to_datetime(df["date"])
if df["date"].duplicated().any():
    df = df.groupby("date", as_index=False)["u60"].mean()
df["doy"] = df["date"].dt.strftime("%m-%d")
clim = df.groupby("doy")["u60"].mean()
tnum = (df["date"] - df["date"].min()).dt.days.values.astype(float)
resid = df["u60"].values - df["doy"].map(clim).values     # anomaly
x = resid - np.polyval(np.polyfit(tnum, resid, 1), tnum)  # detrend on anomaly
dates = (df["date"].dt.year * 10000 + df["date"].dt.month * 100
         + df["date"].dt.day).values

dx = np.diff(x, prepend=x[0])
sig = float(np.std(dx[1:]))
pulse = np.abs(dx) > 2 * sig
pulse[0] = False

print(f"series: u60 @10hPa, {len(x)} daily samples "
      f"({df['date'].min().date()} .. {df['date'].max().date()}), NCEP/NCAR R1")
tau1 = acf_efold(x, dates)
tau2 = dcor_efold(x, dates, rng=np.random.default_rng(SEED))
tau3 = ringdown_efold(x, dates, pulse)
drop = dx < -2 * sig
drop[0] = False
tau4 = recovery_days(x, np.where(drop)[0])
print(f"tau1_acf       = {tau1:.2f} d   (primary)")
print(f"tau2_dcor      = {tau2:.2f} d   (primary)")
print(f"tau3_ringdown  = {tau3:.2f} d   (primary)")
print(f"tau4_recovery  = {tau4:.2f} d   (auxiliary, right-censored -- never "
      f"evidence alone)")
vals = [v for v in (tau1, tau2, tau3) if not np.isnan(v)]
med = float(np.median([v for v in (tau1, tau2, tau3, tau4) if not np.isnan(v)]))
print(f"median (all four, matching the note) = {med:.2f} d")
print(f"primary within median +-30%: "
      f"{sum(1 for v in vals if abs(v - med) <= 0.3 * med)}/{len(vals)}")

# shuffle null: destroy time structure within each season year
rng = np.random.default_rng(SEED)
x_sh = x.copy()
for y in np.unique(dates // 10000):
    m = np.where(dates // 10000 == y)[0]
    x_sh[m] = rng.permutation(x_sh[m])
dsh = np.diff(x_sh, prepend=x_sh[0])
null3 = ringdown_efold(x_sh, dates, np.abs(dsh) > 2 * sig)
print(f"shuffle-null tau3 = {null3:.2f} d  (must collapse far below the "
      f"estimates above)")
