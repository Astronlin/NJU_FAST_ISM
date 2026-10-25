"""Shared editable RFI reference, packaged with fastfish and frozen per run."""
from fastfish._paths import DATA_ROOT
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import math

DEFAULT_CATALOG = DATA_ROOT / 'rfi_catalog.json'


def catalog_path(path='packaged', base_dir=None):
    """Resolve custom paths relative to the example config, not the working dir."""
    if path is None or str(path) == 'packaged':
        return DEFAULT_CATALOG.resolve()
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = Path(base_dir if base_dir is not None else Path.cwd()) / p
    return p.resolve()


def _positive(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(label + ' must be a finite positive number')


def validate_catalog(catalog):
    if not isinstance(catalog, dict) or catalog.get('schema_version') != 1:
        raise ValueError('Unsupported RFI catalog schema; expected schema_version=1')
    if not isinstance(catalog.get('catalog_version'), str) or not catalog['catalog_version']:
        raise ValueError('RFI catalog needs a nonempty catalog_version')
    if catalog.get('frequency_unit') != 'MHz' or catalog.get('frequency_frame') != 'topocentric':
        raise ValueError('RFI catalog frequencies must be topocentric MHz')
    try:
        pattern = catalog['patterns']['repeated_broad_features']
        spacings = pattern['spacing_MHz']
        tolerance = pattern['match_tolerance_MHz']
    except (KeyError, TypeError):
        raise ValueError('RFI catalog is missing the repeated_broad_features pattern') from None
    if not isinstance(spacings, list) or not spacings:
        raise ValueError('RFI pattern needs at least one frequency spacing')
    for value in spacings:
        _positive(value, 'RFI frequency spacing')
    _positive(tolerance, 'RFI spacing tolerance')
    profiles = catalog.get('profiles')
    if not isinstance(profiles, dict):
        raise ValueError('RFI catalog profiles must be an object')
    for name, profile in profiles.items():
        if not isinstance(profile, dict) or not isinstance(profile.get('scan_priors'), dict):
            raise ValueError('RFI profile needs scan_priors: ' + name)
        for scan, entries in profile['scan_priors'].items():
            if not isinstance(entries, list):
                raise ValueError('RFI scan priors must be a list: ' + scan)
            for entry in entries:
                if not isinstance(entry, dict):
                    raise ValueError('RFI feature must be an object: ' + scan)
                _positive(entry.get('centre_MHz'), 'RFI centre_MHz')
                partners = entry.get('partners_MHz', [])
                if not isinstance(partners, list):
                    raise ValueError('RFI partners_MHz must be a list')
                for value in partners:
                    _positive(value, 'RFI partner frequency')
    published=catalog.get('published_statistics')
    if published is not None:
        if not isinstance(published,dict) or published.get('mask_policy')!='reference_only':
            raise ValueError('Published statistics must remain reference-only')
        path=published.get('path')
        if not isinstance(path,str) or Path(path).name!=path or not path.endswith('.json'):
            raise ValueError('Published statistics must name a sibling JSON file')
        digest=published.get('sha256','')
        if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Published statistics require a SHA256 fingerprint')
    return catalog


def load_catalog(path='packaged', base_dir=None):
    """Load one reference; no observation-specific profile is applied implicitly."""
    return validate_catalog(json.loads(catalog_path(path, base_dir).read_text()))


def catalog_fingerprint(path='packaged', base_dir=None):
    p = catalog_path(path, base_dir)
    return hashlib.sha256(p.read_bytes()).hexdigest()


def profile_priors(catalog, profile_name, scans):
    """Select only a named observation profile; never reuse IC443 for other scans."""
    validate_catalog(catalog)
    if profile_name not in catalog['profiles']:
        raise ValueError('Unknown RFI profile: ' + str(profile_name))
    priors = catalog['profiles'][profile_name]['scan_priors']
    missing = [scan for scan in scans if scan not in priors]
    if missing:
        raise ValueError('RFI profile '+str(profile_name)+' has no priors for '+', '.join(missing))
    return {scan: deepcopy(priors[scan]) for scan in scans}


def load_published_statistics(catalog=None, base_dir=None):
    """Load verified historical statistics for inspection, never automatic masking.

    For a custom catalog dictionary, base_dir is its catalog-file directory.
    The source table contains astronomical H I; high occupancy is not proof of RFI.
    """
    catalog=load_catalog() if catalog is None else validate_catalog(catalog)
    record=catalog.get('published_statistics')
    if record is None:raise ValueError('This catalog has no published statistics')
    folder=DEFAULT_CATALOG.parent if base_dir is None else Path(base_dir)
    payload=(folder/record['path']).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=record['sha256']:
        raise ValueError('Published RFI statistics checksum mismatch')
    table=json.loads(payload);rows=table.get('rows',[])
    if len(rows)!=5001:raise ValueError('Expected the complete 5001-row published table')
    for index,row in enumerate(rows):
        if len(row)!=5 or not all(isinstance(v,(float,int)) and math.isfinite(v) for v in row):
            raise ValueError('Invalid published RFI statistics row')
        if abs(row[0]-(1000+.1*index))>1e-6 or not all(0<=row[k]<=100 for k in (1,3)) or min(row[2],row[4])<0:
            raise ValueError('Invalid published RFI statistics units/grid')
    return table


def lookup_published_statistics(frequency_MHz, catalog=None, base_dir=None):
    """Return the nearest historical sample and cautions; not a detection result."""
    frequency=float(frequency_MHz)
    if not math.isfinite(frequency) or not 1000<=frequency<=1500:
        raise ValueError('Published table covers 1000 through 1500 MHz')
    table=load_published_statistics(catalog,base_dir)
    row=table['rows'][int(math.floor((frequency-1000)/.1+.5))]
    return dict(zip(table['columns'],row),query_frequency_MHz=frequency,
        offset_MHz=frequency-row[0],mask_policy='reference_only',
        caution='Different epochs; Galactic H I is present near 1420.5 MHz. Also protect the target H I range. No RFI identification is implied.')
