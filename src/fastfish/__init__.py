"""FAST Pulsar+Spectral calibration and imaging, step by step."""
__version__ = "0.1.0"

# Runs before any pipeline submodule can convert UTC/TAI timestamps.
from .time_tables import initialize_leap_seconds
initialize_leap_seconds()
