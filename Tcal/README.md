# External noise-diode calibration references

This directory sits alongside `FISH_scripts/` and holds shared instrument calibration inputs, rather than observation products.

`CAL.20201014.low.W.fits` contains frequency-dependent noise-diode temperatures for 19 beams and two polarizations. It was copied from the reference path used by the original example notebook. Its internal FROM field points to `lowcal_20201014_W_tny.npz`; the original provider and applicability to other observations have not been established here.

FITS files in this directory are ignored by Git. The example download manifest records the remote source, byte count, and SHA-256 checksum. Run `python3 Example/download_inputs.py` from the repository root to retrieve or verify the example inputs.
