"""Load the user's cached leap seconds without downloading or changing EOP data."""
from fastfish._paths import DATA_ROOT
from pathlib import Path
import warnings


def initialize_leap_seconds():
    """Initialize each new process before pipeline time conversions.

    No cached file leaves Astropy's normal behavior intact. Astropy is optional
    for the remote Qt client, which only needs Qt/NumPy on the local computer.
    """
    path=Path.home()/'.astropy'/'Leap_Second.dat'
    if not path.is_file():return False
    try:
        from astropy.utils import iers
    except ImportError:
        return False
    previous=iers.conf.system_leap_second_file
    try:
        iers.conf.system_leap_second_file=str(path)
        table=iers.LeapSeconds.open(str(path))
        table.update_erfa_leap_seconds()
    except Exception as exc:
        iers.conf.system_leap_second_file=previous
        warnings.warn('Could not initialize cached leap seconds from {}: {}'.format(path,exc),RuntimeWarning)
        return False
    return True


from functools import lru_cache

@lru_cache(maxsize=1)
def initialize_earth_orientation():
    """Pin the same packaged EOP table when deferred Doppler starts in a worker."""
    from astropy.table import Table
    from astropy.utils import iers
    iers.conf.auto_download=False
    iers.conf.auto_max_age=None
    iers.earth_orientation_table.set(iers.IERS_B(Table.read(DATA_ROOT / 'iers_b.fits')))
