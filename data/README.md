# Validation data

## `vortex_daily_hist.csv` — polar vortex u60 series (the τ ≈ 12.3 d validation)

- **What**: daily 60°N zonal-mean zonal wind (`u60`) at 10 hPa, 1979-10-01 .. 2025-11-30, 2,867 samples
- **Source**: NCEP/NCAR Reanalysis 1 (NOAA PSL THREDDS OPeNDAP), public domain data
- **Structure**: 46 ragged Oct–Nov segments (one per year), **not** a continuous series — see the caliber note in `../examples/reproduce_vortex_tau.py`
- **Used for**: the dual-operator note's empirical baseline — τ₁ 11.28 d / τ₂ 10.81 d / τ₃ 13.32 d (median 12.3 d), shuffle-null collapses to 0.9 d

Reproduce: `python examples/reproduce_vortex_tau.py` (numpy + pandas only).

Other columns (`zcap`, `cent_lat`, `D_fold`, …) are polar-vortex structure
diagnostics from the broader research project; not needed for the τ reproduction.

Validation on ECG uses **MIT-BIH Arrhythmia Database** (public, PhysioNet);
bearing vibration data comes from a public bearing test dataset — see
`manuscript.pdf` for the exact references.
