# Package-local FISH example

Open `FISH_example.ipynb` from this folder or the repository root. The runnable example starts from extracted pulsar CAL products; the original notebook is preserved in its documentation appendix.

## Downloaded subset

Sun15_56_ra, beam M01: one raw noise-off spectral file (`N_0002`), 17 extracted CAL files (`0028`–`0044`), the scan pointing table, and the original Tcal file. Total: 936,391,597 bytes. No raw pulsar files or already calibrated spectral files are needed. `input_manifest.json` records server paths, local paths, sizes, and remote SHA-256 checksums.

```text
data/cal/Sun15_56_ra/M01/    # extracted CAL inputs
data/spec/Sun15_56_ra/      # raw spectral input
data/pointing/Sun15_56_ra.xlsx
../Tcal/CAL.20201014.low.W.fits  # shared instrument calibration reference
output/                     # generated products
```

The CAL subset covers the spectral timestamps and includes 60 seconds of preceding smoothing context. Its earlier boundary differs from the full scan, so calibration need not be numerically identical to the full-scan product. Input existence checks do not validate scientific accuracy. No observation data are tracked by Git.

Dependencies: Python, Jupyter, NumPy, SciPy, Matplotlib, Astropy, pandas, openpyxl, and cygrid. cygrid may require platform-specific installation.

The grid-list reader now accepts a single spectral file. Calibration algorithms are unchanged. Known limitations include unchecked interpolation extrapolation, hard-coded diode correction, the gridder spectral WCS offset, and per-file Doppler correction. RFI flagging and baseline subtraction remain unimplemented. Review these before scientific use. The legacy raw PSR extraction route remains documented in the appendix.

The external Tcal reference is shared at package level in `Tcal/`, alongside `FISH_scripts/`. The example manifest and download script use this location.
