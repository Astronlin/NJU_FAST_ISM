"""FAST Pulsar+Spectral calibration and imaging, step by step."""
from pathlib import Path as _Path

__version__ = "0.1.0"

# Resolve former flat imports and python -m commands through small aliases.
# Functional modules import the organized subpackages directly.
__path__.append(str(_Path(__file__).with_name('_compat')))

# Runs before any pipeline submodule can convert UTC/TAI timestamps.
from .pointing.time_tables import initialize_leap_seconds
initialize_leap_seconds()
