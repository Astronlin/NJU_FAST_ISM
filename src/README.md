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

## Developer quick start

Run these commands from the package root (the directory containing `pyproject.toml`,
not from `src/`). Use the same Python environment for installation and execution:

```bash
python3 -m pip install -e '.[qt]'
python3 -c "import fastfish; print(fastfish.__file__)"
fish --help
```

Editable installation loads your source edits directly. Reopen the shell/viewer
and restart preview workers after editing; an existing process retains imported
code and cached results. For a separate development checkout, use a separate
Python environment so its `fish` command cannot replace a colleague's installation.
The Qt extra is for the desktop viewer; the web interface uses the separate
`interactive` extra. See `pyproject.toml` for dependencies and entry points.

## Read the code in this order

1. `commands/shell.py`: `FishShell`, `defaults`, `validate_params`, `_run`.
   This translates commands into session changes, previews or batch recipes.
2. `processing/session.py`: session JSON, revision checks and catalogue loading.
3. `processing/preview.py`: `PreviewDataset` and shared calibration, baseline and
   resampling helpers. This is the bridge between selected input rows and science.
4. `processing/stages.py`: `fingerprint`, `run_group`, `run_batch`; restartable
   per-scan/beam records, upstream dependencies and multiprocessing.
5. `processing/pipeline.py`: cube geometry, gridding, weighted combination,
   validation, configuration helpers and source hashing.
6. `viewer/worker.py`, then `viewer/qt_viewer.py`: preview requests/results and
   their presentation. `viewer/usability.py` owns guided stages, draft feedback,
   plot-range dialogs, fit details and copyable run commands.

Normal interactive flow:

```text
project / file in -> metadata -> go where -> selection + pointing cache
                             -> tsyscal -> applycal -> baseline -> table -> xy_map
Qt controls -> preview worker -> PreviewDataset -> shared science helpers -> plots
FISH stage run -> frozen recipe -> stage workers -> compact records/disk arrays
                                                    -> cube -> validation
```

`file in` starts from file metadata, rather than loading every native spectrum.
`go where` calculates/reuses coordinates. The public stage names differ from the
internal batch names: `tsyscal` is `reference`, `applycal` is `calibration`,
`table` is `resample`, and `xy_map` is `cube`.
`processing/correction.py` is a separate, provenance-controlled downstream
comparison path; it is not the normal configurable-component viewer workflow.

## Where to make common changes

| Desired change | Start here | Also check |
| --- | --- | --- |
| Add a terminal command | `commands/shell.py`, `FishShell.do_*` | Help text, errors and command-file behavior. |
| Change labels, layout, dialogs or guidance | `viewer/qt_viewer.py`, `viewer/usability.py` | Small windows, detached settings, both themes. |
| Change zoom keys or footprint selection | `viewer/keyboard.py`, `viewer/footprint.py` | Rectangle/polygon selection and large point sets. |
| Add a preview result or array | `viewer/worker.py::preview_arrays` | Qt receiver, stage-specific requests and transport size. |
| Add a processing parameter | Shell `defaults` / `validate_params` | UI draft/sync, worker application, cache keys and stage fingerprints. |
| Change gain estimation | `processing/preview.py::calibrate`, `processing/spectral.py::frequency_gain` | Compact saved coefficients and reconstruction. |
| Change ON/OFF or Tpsr smoothing | `calibration/reference.py`, `calibration/extracted.py`, `calibration/tsys.py` | Distinct `n_t_cal` and `n_t_apply` controls and timestamp boundaries. |
| Add a baseline function | `baseline/components.py` | Normalization, basis, reconstruction, records and UI editor. |
| Change RFI search | `rfi/detection.py`, `rfi/catalog.py` | Protected H I, absent priors, AA/BB and mask propagation. |
| Change spectral output | `processing/spectral.py`, preview `output_grid` / `resample` | Doppler scales, channel edges, coverage and RFI fractions. |
| Change cube bounds/frame/kernel | `processing/pipeline.py` | Pointing support, WCS, memory estimates and direct/merge validation. |
| Change input discovery or time handling | `io/`, `pointing/` | SCI midpoints, CAL times, telemetry and cache invalidation. |
| Change worker admission | `processing/scheduler.py` | Per-worker memory, I/O limits and extra output fields. |

Before editing a function, search its callers, including reconstruction and batch
paths. For example, from the package root:

```bash
rg -n 'fit_components|reconstruct_components' src/fastfish
rg -n 'n_t_apply|n_t_cal' src/fastfish
rg -n 'fingerprint|source_hash' src/fastfish/processing
```

## Example: customize the baseline without changing Python

The existing component interface already supports ordered sin/cos and polynomial
functions. In FISH, after opening a session:

```text
set baseline_components [{"function":"sincos","range_kms":null,"period_MHz":null},{"function":"polynomial","degree":2,"range_kms":[-1500,1500]}]
baseline
```

This previews an automatic-period sine/cosine followed by a quadratic residual.
Use `baseline run` only when ready to process the loaded selection. The Qt
component table edits the same parameter. A range is both the fitting and
subtraction domain; components do not extrapolate outside it. H I protection,
detected RFI and extra exclusions remove fit anchors, while the model can still
be evaluated inside the protected interval. Components fit sequentially.

To add a genuinely new function in Python, update all of:

- `normalize_components`: validate and serialize the new specification.
- `_basis` and fitting logic: define units, normalization and support requirements.
- `evaluate_record` / reconstruction: reproduce the fitted result from records.
- Qt component editor and its parameter serialization: expose valid controls.
- Tests for known signals, insufficient support and protected-line invariance.

Avoid fitting in the GUI itself. Preview and batch must call the same science
implementation; saved records must reconstruct the model after the fit process
has exited.

## Example: add a parameter safely

Suppose a new baseline option changes science results. Add its default and type/
range checks in the shell, carry it through session parameters into the effective
configuration, then expose it in the Qt draft reader and accepted-state sync.
Include it in relevant preview cache keys and the baseline/downstream stage
fingerprints in `processing/stages.py`. Check any separate correction path that
also uses it. A UI-only display option should stay out of science parameters.

The viewer distinguishes saved settings from draft edits. A delayed response must
not overwrite newer edits. `Update preview` submits validated settings; session
revision checks prevent stale writers replacing a newer session. Preserve these
contracts when adding controls. Changing only a label should not trigger a fit.

## Data and scientific contracts to preserve

- Native frequencies in the science helpers are generally **GHz**; catalogue and
  sinusoidal-period fields use **MHz**; velocities use **km/s**. Check each API
  before converting. Pointing angles are degrees and internal times commonly MJD.
- Typical native polarization arrays are `(timestamps, channels, 2)` for AA/BB.
  Baselines are fitted independently per polarization before their final average.
  Inspect the called reader/helper before assuming an array orientation.
- Science integration midpoints are `UTOBS + EXPOSURE/2`. CAL timestamps already
  represent midpoints. Pointing, interpolation and Doppler use consistent times.
- `n_t_cal` smooths ON−OFF for temperature derivation; `n_t_apply` smooths the
  temperature reference before application. Display averaging changes plots only.
- RFI catalogues provide search priors. Only observed detections become baseline
  exclusions; a catalogue entry alone is not a flag. Preserve protected H I.
  Quality-rejected spectra and RFI flags are separate concepts.
- `n_b` controls overlap rebinning. `n_g` is a target Gaussian FWHM in native-channel
  units; extra smoothing accounts for the assumed native response. Interpolation,
  smoothing and rebinning must not be substituted without checking the response.
- Compact calibration/baseline records are intentional. Reconstruct only needed
  chunks; do not retain full native arrays for an entire observation.
- Unsupported diagnostics must remain explicitly unavailable, not silently zero.
  Arithmetic cube validation does not prove baseline truth or absolute calibration.

## Caches, provenance and memory

Observation inputs live under project `org_spec/`, `org_cal/`, `feed_info/`.
Generated sessions, coordinate caches, stage records and cubes belong inside that
project's `fish_runs/`. Do not write observation caches into this source package.
Use `_paths.py` for packaged resources rather than hard-coded machine paths.

Stage paths include fingerprints of source/settings/input/catalogue information.
Never bypass a mismatch or rename old products to make them appear compatible.
For an intentional reuse, document exactly which records remain compatible and
verify their provenance. Preserve completed reference products during comparisons.
`processing/pipeline.py::source_hash` currently hashes `.py`, `.fits` and `.json`
files beneath the Python package. This README is outside that hash scope; editing
it alone does not invalidate scientific stage fingerprints. Inspect the hash
implementation before changing its scope or resuming a modified recipe.

Process spectra in bounded chunks. Worker count is a ceiling, not a guarantee that
all workers launch: scheduler admission considers memory and I/O limits. When
adding arrays/fields, update memory estimates and propagation through gridding,
combination and validation. Keep Qt callbacks short; expensive work belongs in the
preview worker. Reuse cached reads only when input, parameters and masks match.

## Verification before sharing a change

For a development environment with dependencies installed:

```bash
python3 -m compileall -q src/fastfish
python3 -m fastfish.diagnostics.selfcheck
```

The synthetic self-check covers selected numerical routines; it does not fully
exercise every viewer feature or component-list configuration. Add checks suited
to the changed behavior rather than relying on this alone:

1. Compare a small real chunk before/after, AA and BB separately; verify masks,
   coefficients, reconstruction and spectral output. For performance-only changes,
   require equal outputs or a justified numerical tolerance.
2. For baseline changes, inject a known line/ripple and vary protected H I while
   holding masks fixed. Check recovery, fit support and held-out channels.
3. For UI changes, test navigation, draft/reset, delayed replies, error recovery,
   narrow windows and detached settings. `QT_QPA_PLATFORM=offscreen` can support
   automated Qt checks, but also inspect rendered screenshots.
4. For cube changes, compare weighted merging with direct gridding on a bounded
   selection, then check coverage, finite support, RFI bounds, WCS and missingness.
5. Confirm existing observation files and completed products remain untouched.

Commit/deploy source and relevant reference data only. Keep temporary previews,
benchmarks and large generated test data outside the software package. Verify the
installed import path and file hashes on each destination; restart old processes.
Do not change code underneath an active science run. Keep a recoverable previous
version and describe changed behavior, validation and remaining limitations.


## Runnable minimal calibration example

See [examples/minimal_calibration/README.md](../examples/minimal_calibration/README.md) for the bundled real Sun17_A_01/M01 extract, readable runner and step-by-step notebook. It calls production calibration functions and checks numerical repeatability against an original-input calculation. The example is included in the source distribution; it is not installed into the Python wheel.


The minimal Sun17 example includes `prepare_viewer.py` for portable sessions and
`run_workflow.py` for the normal four-stage run/wait sequence through resampling.
Its `reference_context_max_rows` limit is enforced by parameter validation and
reference reads. `omit_gridding` blocks only cube launch, omits cube commands from
the Qt command dialog/terminal guide, and labels cube status omitted. Older
`preview_only` sessions retain their stricter batch guard. Ordinary configurations
without these optional flags retain the full workflow. The scheduler detects
macOS available pages using vm_stat because SC_AVPHYS_PAGES is unavailable there;
non-Linux job ownership uses ps command/recipe identity, not PID existence alone.
