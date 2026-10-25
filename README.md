# FAST Pulsar+Spectral pipeline

Source for FAST calibration, pointing, RFI-aware baseline fitting, spectral resampling and imaging.

- `src/fastfish/`: source grouped by function; see [src/README.md](src/README.md).
- `src/fastfish/data/`: shared RFI catalog, published statistics and IERS table.
- `Tcal/`: external noise-diode temperature reference.
- `pyproject.toml`: dependencies and command entry points.
- `bin/`: Mac viewer launcher, where available.

Install in an external environment with `python3 -m pip install -e ".[qt]"`, then start `fish`.

```text
project /bb8/lingrui/my_projects/FAST_IC443
file in org_spec/Dec+2150_06_05/20210204/*_N_0002.fits
list
go where
```

Use `help` for command details. `tsyscal`, `applycal`, `baseline` and `table` open previews; append `run` for batch processing, followed by `wait` before the next stage. `xy_map run --validate` grids and validates the prepared table. Start Qt from a graphical desktop/VNC terminal or a working X11-forwarded session.

Keep inputs, sessions, caches and products outside this package. Data projects contain `org_spec/`, `org_cal/`, and `feed_info/`. Coordinates are cached once in the project's `fish_runs/cache/coordinates/`; project imports write sessions under the data project’s `fish_runs/` directory.

Bundled examples and old documentation have been removed. Use project/file inputs or your own configuration. The Mac launcher uses the external `../fastfish_runtime/qt` environment, or `FISH_QT_PYTHON` when set.

License: GPL; see `LICENSE`.

Batch workers are memory-admitted with a ceiling of 64. One pilot task measures peak memory before scaling up. Scheduling uses MemAvailable and visible cgroup headroom, reserves at least 2 GiB or 10% of total RAM, and limits estimated task memory to 60% of total RAM (or configured memory_budget_gib). New tasks pause under pressure; no-admission timeout defaults to 300 s. No idle heavy worker pool is retained. Gridding runs independent scan/channel blocks with one cygrid thread per task, and combination remains serial/block-bounded. This is conservative admission, not a hard memory limit. Scheduler JSON records active workers and estimates. All caches and generated products belong under the data project’s fish_runs directory; source snapshots are not copied there.

I/O-heavy stages also have an independent concurrency limit: io_workers defaults to 4 for reference, calibration, baseline, resampling, reduction and gridding; stage_io_workers can override individual stages. The global worker ceiling remains 64, and memory admission may lower concurrency further. Scheduler status/logs show both limits. Contiguous SCI row selections use sequential read-ahead; sparse rows and timestamp-only metadata retain random-access behavior. More threads do not bypass storage throughput limits.

Science timestamps: science_timestamp_reference defaults to start. Effective time is UTOBS + each row’s EXPOSURE/2; this same midpoint drives calibration, pointing and Doppler correction. Raw FITS timestamps are unchanged. Set science_timestamp_reference to midpoint only for data already tagged at the center. Scalar metadata records raw and effective bounds, coordinate files record the convention, and FITS cubes include SCITIME/RAWTIME/TIMEOFF provenance. Coordinate cache fingerprints change when the time reader changes; old cached coordinates are not reused for new preparations. Existing completed historical products are not rewritten.

Batch calibration saves fitted coefficients without constructing unused native gain/calibrated arrays. Batch baseline fitting omits preview model/residual stacks, and batches the coarse sine/cosine period search; independent AA/BB fits, masks, clipping and nonlinear refinement remain. Four real Sun17 chunks gave identical coefficients/reconstructed spectra, with measured calibration/baseline computation about 1.7x/1.6x faster and process peak RSS 714 to 380 MiB. These are bounded checks, not whole-run benchmarks.

Cube geometry: coordinate_frame is icrs (default) or galactic in the project configuration. With header=null, go where derives a TAN footprint from valid pointings in that frame, handles longitude wrapping, and includes the actual Gaussian kernel support before rounding to pixel edges. An explicit header retains user-specified geometry. Gridding and its independent validation transform ICRS pointings consistently; Doppler calculations continue to use ICRS. Internal footprint gaps can remain blank. Choose the frame before preparing a new session.

Published RFI references: see [RFI_REFERENCE.md](RFI_REFERENCE.md). The shared catalog now links the verified 5001-row CC0 FAST statistics, selected official site bands, epoch-qualified harmonic patterns and explicit IC443 comparisons. Published statistics remain advisory and include Galactic H I; no automatic masks or historical products are changed.

A downstream comparison can reuse frozen calibration with `python3 -m fastfish.correction --source-run SOURCE --run-root COMPARISON --stage baseline`, then `--stage resample`, then `--stage cube`. It requires the source run's reference-quality audit; this is the reviewed Sun17 comparison workflow, not a general automatic repair. Original products remain immutable. The comparison detects RFI separately in AA/BB, conditionally fits a 1.04–1.14 MHz ripple using withheld frequency blocks, and excludes spectra whose calibration reference has fewer than 90 percent positive finite gain anchors. Rejected rows are recorded separately from RFI and removed consistently from coordinates and table arrays. New cube headers record the ripple/quality policy. Protected H I never anchors the ripple; fit-window sensitivity is not an absolute baseline uncertainty.

The downstream correction keeps the primary ±1500 km/s fit strict. If the alternative ±1000 km/s sensitivity fit has fewer than eight clean anchors on either side, its coefficients are marked unavailable for that row. The final baseline_sensitivity.fits is NaN wherever an unavailable diagnostic contributes; baseline_sensitivity_missing_fraction.fits records that geometric weighted fraction. Science spectra and coverage are retained when the primary fit is valid. A separate resume amendment records compatible diagnostic-only source changes and hashes reused products.


## Configurable baseline components (Qt and FISH stages)

The fixed broad/local controls are replaced by an ordered component table. Use **+ Add baseline component**, choose `sin/cos` or `polynomial`, edit its range, and use the arrows or × to reorder/remove it. Click **Update preview** to preview. The initial two editable entries reproduce the former model choices; neither is compulsory. The same component list is used by `baseline run` and `table run`.

- Range is `full` or `min:max` in LSRK km/s, e.g. `-1500:1500`. Every component fits AND subtracts only in its own range; no extrapolation beyond it. Fit ranges can be asymmetric. All components exclude protected H I, detected RFI and extra baseline-only exclusions.
- Polynomial degree n: 0 constant, 1 linear, 2 quadratic, 3 cubic, up to 20; a stable Chebyshev representation is used internally. At most 16 components. Insufficient support/ill-conditioned fits fail explicitly rather than silently fitting.
- Sine/cosine includes a constant plus sine and cosine. Period is in native-frequency MHz: `auto`, a fixed value such as `1.09`, or a search interval such as `1.04:1.14`. Automatic periods are shared per short chunk/polarization; coefficients are independent at each time and in AA/BB.
- Components fit sequential residuals in list order; this is not a simultaneous joint fit. Each component and the total model are plotted. Preview caching and batch fingerprints include the full component list. Batch records store compact coefficients and periods, not full native model cubes.

CLI equivalent: `set baseline_components [{"function":"sincos","range_kms":null,"period_MHz":null},{"function":"polynomial","degree":2,"range_kms":[-1500,1500]}]`. An empty list applies no baseline correction. Existing configurations without the list are presented with editable sine/cosine + degree-1 entries; old saved science records remain preserved. Restart FISH/viewer after updating.

The old fixed +/-1000 vs +/-1500 sensitivity statistic has no general meaning for arbitrary component lists. New component-list cubes mark that diagnostic unavailable (NaN plus explicit missing-fraction product), rather than reporting false zero sensitivity. The recipe records exact model choices. No existing science cube has been refitted by this software update.

Spectrum navigation: Z selects a horizontal-only zoom rectangle; Shift+Z selects both velocity and intensity ranges. B restores the previous view, F restores full extent, Escape cancels. For exact vertical limits in an already-open window, right-click the spectrum -> Y axis -> Manual, enter minimum/maximum and press Enter. Restart the viewer to load updated keyboard bindings.

Spectrum style: select Centered steps (default) or Lines in the Qt toolbar. This affects observed/calibrated spectra, residuals and resampled/region previews only; fitted components and total models stay smooth, and Tsys time curves remain lines. Steps use channel midpoint edges and exact boundaries of display-averaging groups, including the incomplete last group. Set Display mean points to 1 to inspect native channels. Style changes use cached display arrays without changing calibration, baseline fits, resampling or saved science data. Reopen the viewer after a software update.

Apply calibration preview: upper panel shows the selected raw SCI AA/BB counts; lower panel shows calibrated AA/BB temperature in kelvin (Mean averages each pair in its own units). Both use the same selected timestamp/native channels and LSRK velocity axis. Display styles/averaging affect both consistently; baseline models are hidden on this page. Raw arrays reuse the same bounded preview chunk cache; arrays and Python coefficient records both count toward cache admission. No calibration mathematics or saved products are changed. Reopen the viewer to refresh both local UI and remote worker.

Viewer response and layout (2026-10-23): each page now requests only its own stage. Tsys preview reads the reference without SCI calibration; Apply calibration skips baseline/resampling; Baseline skips product resampling. Local worker replies use lossless uncompressed NPZ; SSH uses fast lossless compression. Raw SCI counts are reused when smoothing changes, and region previews reconstruct only selected rows while retaining the original whole-chunk period fits. Oversized entries no longer flush useful small caches. The worker cache remains 128 MiB and the display cache is limited to 32 MiB/12 entries; these budgets do not include transient fitting/CAL allocations or other library caches.

On the same Sun17 spectrum on .150, three-request median warm end-to-end timings including transfer/decoding improved from 0.77 to 0.19 s for Baseline, 0.72 to 0.06 s for Apply calibration, and 0.79 to 0.20 s for Table. Cached Tsys replies take about 0.003 s. First-time calibration/fitting of a new chunk still takes several seconds. All compared native/calibrated/component/resampled arrays were exactly equal; 48 rendering comparisons and rapid stage-switch tests passed. Plot redraw median improved from 76 to 27 ms in the offscreen test. These are bounded .150 checks, not whole-observation benchmarks.

The viewer uses compact stage-specific settings with a Settings toggle, leaving more space for the side-by-side footprint and spectra. The spectral panels share the velocity/time axis and retain independent vertical ranges. Calibration pages hide baseline/output controls; Baseline shows its component editor, and Table shows output controls. Spectrum styles, display averaging, keyboard zoom and manual vertical ranges are retained. Rapid requests preserve the latest selection and queued parameter edits. Close and reopen the viewer after updating.
The parameter console can also be detached into its own resizable window, with Apply/refit and Attach actions. It reuses the same controls and retains draft values, stage and spectrum selection; detaching/attaching alone does not recompute science. Closing the console safely returns it to the main window.

Footprint scan colors: Scan colors is enabled by default and assigns a stable color to each scan name. Use its arrow menu → Show scan legend for the loaded names and spectrum counts, or turn the button off for a single color. The legend is a separate window so it does not obscure the map. Colors adapt to the plot theme. Changing colors or opening the legend changes display only; selected spectra, displayed point coordinates, calibration and baseline masks stay the same. Reopen the viewer after updating.

Source organization (2026-10-24): functional subpackages replace the flat source directory. Existing commands and former import names remain supported through `_compat/` aliases. See [src/README.md](src/README.md) for the module map. Scientific formulas, defaults and completed products are unchanged; the reorganized source has a new fingerprint.

RFI display key: the spectrum legend labels red bands as Detected RFI (with interval count), purple bands as extra fit exclusions, and gray as protected H I. Detected RFI excludes baseline anchors; channels are retained. No detected intervals is shown as a count of zero, not as proof that all channels are clean. The Tsys time plot hides this spectral key.

Use **RFI catalogue…** to select a package-format JSON catalogue and optional observation profile. In a remote Mac viewer, enter a path on the data server; in a server/local viewer, Browse is available. `packaged` selects the shared default. Leave the profile blank to use published-statistics search priors only; a named profile must cover all loaded scans. JSON frequencies must be topocentric MHz, schema_version 1. Refer to the packaged data/rfi_catalog.json and RFI_REFERENCE.md for the format. If published_statistics is declared, its verified sibling JSON must accompany the input catalogue. Invalid catalogues leave the active session unchanged.

FISH terminal equivalent: `catalog /absolute/path/rfi_catalog.json --profile IC443_2021`, `catalog packaged`, or `catalog` to inspect the choice. Loading validates and freezes the catalogue (and published table) under this session's fish_runs/rfi_catalogues, then selects a new configuration for previews and future stage runs. Original configs/catalogues/products are preserved. Load before stage runs; active jobs block catalogue changes. A different catalogue/config creates different stage fingerprints, so complete required upstream stages with the new settings before running downstream stages. Catalogue entries are search priors, never automatic flags. Detection thresholds and H I protection are unchanged.

Hover over the Baseline title/plot when it says **period at limit** for the component, AA/BB polarization, fitted MHz period and allowed search interval. This reports a constrained automatic search, not an automatic fit failure. Reopen the viewer and FISH shell after updating.


## Guided viewer workflow

Choose Stage, then use **Next step** to move from Footprint through Tsys calibration, Apply calibration, Baseline and Table. The hint below the toolbar explains the current step. The selected scan, beam and raw row are shown above the plots.

Parameter edits are marked **Unapplied edits**. Click **Update preview** to save and preview them, or **Reset edits** to restore the session settings. New edits survive an older preview finishing in the background. Preview activity shows elapsed time; error and period-boundary notices offer **Details…**.

Use **View → Plot ranges…** for exact horizontal and independent upper/lower vertical plot limits. These are display limits. The Table page's **Output velocity** limits instead control the spectral product. **View → Display tools** reveals spectrum styles, display averaging, colors and image controls. **Workflow → RFI catalogue…** opens catalogue selection. Settings can still be detached.

**Run commands…** provides copyable commands for the saved session; it does not launch processing. Apply pending edits first, then run the commands in FISH. In the terminal, `guide` prints the workflow and current session/run folder. Restart FISH and the viewer after installing the update.
