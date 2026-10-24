# External noise-diode calibration reference

This shared reference sits alongside `src/` and `Example/`. It contains frequency-dependent diode temperatures for 19 beams and two polarizations. It is an instrument input, not a test output.

`CAL.20201014.low.W.fits` is retained from the original example. Its internal FROM field points to `lowcal_20201014_W_tny.npz`; the original provider and applicability to other observations have not been established here. The old server source path is unavailable, so this local reference is preserved. Its size and SHA256 are recorded in `docs/cleanup_20261018.json`.

Existing extracted temperature references already incorporate their diode calibration. Use this file explicitly only when regenerating a reference from extracted CAL, as shown in the notebooks.
