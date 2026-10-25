# Minimal FISH spectral workflow: Sun17_A_01 / M01

Use a compact **real observation** to learn every step from extracted ON/OFF and
science input to the final spectral table. Only spatial gridding is omitted.
No synthetic spectra or known astronomical truth are included.

## Run the complete example

From the package root, in a compatible scientific Python environment:

```bash
python3 -m pip install -e '.[notebook,qt]'
python3 examples/minimal_calibration/run_workflow.py --output "$HOME/fish-suna01-workflow"
```

Use a new empty output directory. The runner uses the normal FISH `run`/`wait`
commands for reference, calibration, baseline and table; it does not maintain a
second implementation of the science. It writes all generated files outside the
released example and uses one worker with two 16-row native chunks.

Open the resulting session in the actual Qt viewer:

```bash
fish-viewer --session-state "$HOME/fish-suna01-workflow/fish_session.json"
```

Choose Tsys calibration → Apply calibration → Baseline → Table in Stage. Click a
real pointing or use Previous/Next, inspect AA/BB, edit a parameter, and Update
preview. Previewing changes does not overwrite previously completed stage products.
To save processing after edits, restore the session in the FISH prompt and run
required stages with the same current settings. Fingerprints distinguish results.

## Inspect first, then save through FISH

Alternatively prepare a session without running the stages:

```bash
python3 examples/minimal_calibration/prepare_viewer.py --output "$HOME/fish-suna01-session"
fish-viewer --session-state "$HOME/fish-suna01-session/fish_session.json"
```

In another terminal in the same environment:

```bash
fish --viewer qt
```

Then in FISH (replace the path with your actual home directory):

```text
restore /your/home/fish-suna01-session/fish_session.json
tsyscal
applycal
baseline
table
```

After inspecting settings, save the four stages:

```text
tsyscal run
wait
applycal run
wait
baseline run
wait
table run
wait
status
```

The compact session permits all these stages. Only `xy_map run` is disabled,
because the short one-beam track is unsuitable for a meaningful spatial cube.
The viewer's Run commands dialog copies exactly the four-stage sequence.

## Input data and exact volume

| Input | Selection | Bytes |
| --- | --- | ---: |
| `data/science.fits` | Real M01 N_0001, zero-based rows1040:1072, AA/BB, 65,536 native channels | 16,796,160 |
| `data/on_off.fits` | M01 CAL_0017, 64 times, OFF/ON, AA/BB, original channels3314:3575 | 279,360 |
| `data/tcal_reference.fits` | Original CAL.20201014.low.W.fits, M01 channels53006:57203 | 106,560 |
| `data/pointings.npz` | Actual 32 RA/Dec positions, midpoint times and LSRK frequency factors | 2,013 |
| **Data total** | Excluding small JSON metadata | **17,184,093 (16.39 MiB)** |

Observation date is 2021-10-10. SCI exposure totals 6.442450944 s (32 integrations
of 0.201326592 s). Native frequency coverage is 1404.784679–1436.034203 MHz. No
native spectral averaging was used when extracting the input; cross-polarizations
were omitted. CAL retains context for both default 10-timestamp filters; Tcal
retains context for its 16-frequency-channel filter. The example rejects larger
time windows, which require a longer extract. Actual telemetry-derived pointings
and per-spectrum LSRK factors were saved from the verified Sun17 reduction; they
are not fabricated from the SCI header's placeholder coordinates. Viewer raw rows
0..31 correspond to original rows1040..1071.

`data/provenance.json` records basenames, slices, checksums and extraction checks.
The original data remain unchanged on the server. For comparison, the whole A01
scan has 82,037,318,400 bytes of narrow-band SCI (57 files) and 10,288,290,240 bytes
of extracted CAL (1,235 files) across 19 beams; these are unnecessary here.

## Processing stages and scientific conventions

1. **Load data and positions.** SCI is `(time, channel, AA/BB)`; reference arrays
   are `(time, AA/BB, CAL channel)`. Convert original SCI UTOBS start times to
   midpoints using EXPOSURE/2; CAL timestamps already denote midpoints.
2. **Temperature reference.** Production `rebuild_reference` normalizes ON−OFF,
   filters it with `n_t_cal`, smooths/interpolates Tcal, derives Tpsr and retains
   float32 reference rounding. `n_t_apply` independently filters Tpsr before
   interpolation to SCI times. Neither control time-averages SCI spectra.
3. **Calibration and RFI.** Production `calibrate` searches observed spectra using
   catalogue priors, saves detected intervals, and fits positive frequency-dependent
   gains outside the gain-anchor H I exclusion (±700 km/s) and detected RFI. SCI
   channels are retained. No detections in this chunk is not proof of clean data.
4. **Baseline.** Ordered sin/cos and degree-1 polynomial components are fitted
   independently to AA/BB. Defaults protect ±250 km/s; the first component uses the
   full native band, the local component fits/applies only within ±1500 km/s.
   Detected RFI and extra exclusions do not anchor the fit. Models are subtracted
   inside protected H I without using that H I as fitting data. Inspect the report:
   the default run reports one chunk with a period at its allowed boundary.
5. **Doppler alignment and table.** Saved real LSRK factors align each residual
   spectrum; overlap rebinning uses n_b=2 and the output range ±250 km/s. The output
   has 32×2485 values, not a spatial cube. RFI fractions and explicitly unavailable
   sensitivity diagnostics are retained in the production table records.

Noise-diode timing: the official [FAST backend manual, Section 4, p.4](https://fast.bao.ac.cn/static/uploadfiles/%E5%B8%AE%E5%8A%A9%E6%96%87%E6%A1%A3/FAST_Digital_Backends_and_Science_Data.pdf#page=4)
defines ON/OFF settings in 4 ns units. The original code uses ON=20480 and OFF=28672:
81.92 μs ON, 114.688 μs OFF, and a 196.608 μs cycle. The CAL header confirms a
98.304 μs PSR sample. The extractor averages alternate equal-length samples.
With the assumed alignment the ON-labelled sample receives injection for 5/6 of
its duration, hence D=(ON−OFF)/20480×(20480+28672)/2. This is a 1.2 correction,
not the inverse of the whole-cycle ON fraction. Different timing/phase needs a
new calculation. See also the [CRAFTS paper, Section II.1](https://arxiv.org/html/2412.08173v2#S2.SS1).

## Baseline inspection and parameter exercises

- Inspect the mask bands and AA/BB separately. The green/local model in
  `baseline.png` is plotted alongside the residual; the H I peak is deliberately
  clipped so small offsets are visible. The display averages 32 channels; fitting
  and stored native spectra are unchanged. Local model steps at ±1500 km/s follow
  the explicit no-extrapolation policy.
- In the Qt component editor, change a function/range/degree/period, then Update
  preview. Components fit sequentially. A period at its limit is a diagnostic,
  not proof of success or failure. Open Details and compare alternatives.
- Adjust protected H I or extra exclusions only with a scientific reason. Do not
  label a sky line as RFI just because it is bright. Catalogue entries are priors.
- In FISH, for example: `set baseline_exclusions_kms [[800,850]]` adds a baseline-only
  exclusion, not an RFI identification. This is syntax guidance, not a claim that
  an RFI exists there in this sample.
- Try rebin n_b=2 versus another value or Gaussian n_g. n_g is target Gaussian
  FWHM in native-channel units; n_b is the rebin factor. The method changes the
  spectral response. Display averaging is separate.

## Saved outputs and expected results

- `table_manifest.json`: output frequencies, settings, and disk-backed table paths.
- `stages/reference/...`: temperature-reference arrays and provenance.
- `stages/calibration/...`: gain coefficients and per-chunk RFI intervals.
- `stages/baseline/...`: compact per-polarization component fits and fit reports.
- `stages/resample/...`: residual spectra, RFI fractions, coordinates and diagnostic fields.
- `workflow_summary.json`, `baseline.png`, `table-spectrum.png`: compact checks and
  figures made from saved fits, not a separate re-fit for presentation.

The default table is finite with 32×2485 values. Saved-table consistency and
known input hashes are numerical checks, not an absolute calibration or baseline
truth test. Inside protected ±250 km/s the baseline is constrained by outside
anchors and model assumptions; its true error cannot be measured by the fit alone.
Do not regard the default model as uniquely optimal. Native-channel residual
scatter includes radiometer noise; the displayed channel means have smaller scatter.

The compact CAL/Tcal produces exactly the same default effective reference and
calibrated spectra as surrounding original CAL files/full Tcal on the server
(maximum difference 0 K). Some nonpositive Tpsr values occur in the excluded H I
channel at 1420.28808594 MHz; disabling filters also exposes nonfinite reference
values at 1425.415 and 1426.514 MHz. Production fitting rejects unusable anchors;
fractions/frequencies are reported by the calibration-only script. This is not a
blanket assertion that invalid references are harmless.

`calibrate.py` remains a calibration-only sub-example; its output is topocentric
(scale=1) and does not subtract a baseline. `run_workflow.py` and the viewer use
real LSRK factors and complete the spectral workflow. The notebook covers both.
The optional `--n-t-cal`, `--n-t-apply`, and `--components-file` flags customize a
fresh unattended run. Preserve released expected results while experimenting.
