# Source layout

`fastfish` is divided by responsibility. The installed `fish`, `fish-viewer`,
`fastfish`, and `fastfish-interactive` commands remain the public entry points.

| Directory | Responsibility |
| --- | --- |
| `fastfish/commands/` | The CLASS-like shell and command-line preparation tools. |
| `fastfish/viewer/` | Qt/web interfaces, display helpers and the persistent preview worker. |
| `fastfish/processing/` | Shared preview/batch execution, spectral resampling and cube production. |
| `fastfish/calibration/` | ON/OFF temperature calibration and science reference handling. |
| `fastfish/baseline/` | Masked polynomial, sine/cosine and standing-wave baseline models. |
| `fastfish/rfi/` | Shared advisory catalog and observed-spectrum RFI detection. |
| `fastfish/pointing/` | Telemetry, coordinates, Doppler corrections and time references. |
| `fastfish/io/` | FITS readers and observation-project discovery. |
| `fastfish/diagnostics/` | Calibration/product diagnostics and the synthetic self-check. |
| `fastfish/data/` | Shared RFI JSON references and the IERS FITS table. |
| `fastfish/_compat/` | Aliases for the former flat module names; no scientific implementation. |

The main flow is `commands/shell.py` → `processing/stages.py`, using
`processing/preview.py` for the same science calculations as interactive previews.
`processing/pipeline.py` handles cube geometry, gridding, combination and validation.
The Qt interface is `viewer/qt_viewer.py`; `viewer/worker.py` performs the previews
locally or over SSH. Calibration, baseline, RFI, pointing and FITS algorithms live
in their corresponding directories. New code should use those organized imports.

For example, import `fit_components` from `fastfish.baseline.components`.
The older `fastfish.baseline_components` import and `python -m fastfish.qt_viewer`
still work through compatibility aliases. These aliases point to the same module
objects, rather than maintaining duplicate implementations.

`_paths.py` resolves packaged resources centrally. The noise-diode reference remains
in the package's external `Tcal/` directory. Observation inputs, sessions, caches
and cubes belong to each data project's `fish_runs/`, outside this source tree.

The directory reorganization changes source fingerprints. Completed products are
preserved; existing source-fingerprint safeguards continue to apply to new batch
work. It does not authorize reusing products with incompatible fingerprints.
