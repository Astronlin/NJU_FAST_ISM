#!/usr/bin/env python3
"""Create a portable FISH session for the bundled real-data viewer example.

Writes only a new external workspace; observations are never changed.
"""
from pathlib import Path
import argparse
import json
import shutil
import numpy as np
from astropy.io import fits
from calibrate import HERE, load_inputs, REST_GHZ
from fastfish.commands.shell import defaults, validate_params
from fastfish.io.fits_stream import spectral_metadata
from fastfish.processing.session import write_session


def prepare(output):
    root = Path(output).expanduser().resolve()
    if root == HERE or HERE in root.parents:
        raise ValueError('Choose a workspace outside this example')
    if root.exists() and any(root.iterdir()):
        raise ValueError('Use a new empty workspace; existing files are preserved')
    raw, frequency, times, provenance = load_inputs()
    data = HERE / 'data'
    with np.load(data / 'pointings.npz') as c:
        coords = {k:c[k].copy() for k in c.files}
    np.testing.assert_allclose(coords['mjd'], times, atol=1e-10, rtol=0)
    np.testing.assert_array_equal(coords['ids'], np.arange(len(times)))
    assert np.isfinite(coords['ra']).all() and np.isfinite(coords['dec']).all()
    assert np.isfinite(coords['frequency_scale']).all()
    root.mkdir(parents=True, exist_ok=True)
    calfolder = root/'org_cal/Sun17_A_01/20211010'
    calfolder.mkdir(parents=True, exist_ok=True)
    # Original basename is needed by CAL time-index discovery, not invented data.
    shutil.copy2(data/'on_off.fits', calfolder/provenance['CAL_original_basename'])
    coordinate_file = root/'pointings.npz'
    shutil.copy2(data/'pointings.npz', coordinate_file)
    cfg = dict(omit_gridding=True,scans=['Sun17_A_01'], beams=['M01'], rest_freq_ghz=REST_GHZ,
        chunk_rows=16, gain_degree=3, n_t_cal=10, n_t_apply=10,
        reference_context_max_rows=10, baseline=True, protect_kms=250.,
        local_fit_kms=1500., alternative_fit_kms=1000., baseline_exclusions_kms=[],
        baseline_components=[dict(function='sincos',range_kms=None,period_MHz=None),
            dict(function='polynomial',degree=1,range_kms=[-1500.,1500.])],
        baseline_fit_bin=128, spectral_mode='rebin',n_b=2,n_g=2.,
        native_fwhm_chan=1.,native_df_hz=float(np.median(np.diff(frequency)))*1e9,
        velocity_range_kms=[-250.,250.],rfi_catalog='packaged',rfi_profile=None,
        extracted_cal_root=str(root/'org_cal'), tcal_path=str(data/'tcal_reference.fits'),
        workers=1,pointing_workers=1,io_workers=1,
        science_timestamp_reference='start')
    config = root/'config.json';config.write_text(json.dumps(cfg,indent=2)+'\n')
    entry = spectral_metadata(data/'science.fits')
    selection = root/'selection.json'
    selection.write_text(json.dumps(dict(items=[dict(scan='Sun17_A_01',beam='M01',
        reference=str(root/'reference_from_extracted_cal.fits'),
        files=[dict(entry=entry,coordinates=str(coordinate_file))])]),indent=2)+'\n')
    session = root/'fish_session.json'
    write_session(session,dict(config=str(config),selection=str(selection),root=str(root),
        params=validate_params(cfg,defaults(cfg)),index=0,omit_gridding=True,
        example_original_row_offset=1040,
        example_note='32 real M01 samples. Viewer raw rows 0..31 correspond to original rows 1040..1071. Cached real pointings and LSRK factors are included.'))
    print('Created session:',session)
    print('Viewer raw rows 0..31 correspond to original SCI rows 1040..1071.')
    print('Open with: fish-viewer --session-state '+str(session))
    print('Initial stage: Footprint. Choose Tsys calibration, then Apply calibration.')
    return session


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    prepare(parser.parse_args().output)
