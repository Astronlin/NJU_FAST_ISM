"""Compatibility alias; implementation: fastfish.baseline.standing_wave."""
from importlib import import_module
import sys

if __name__ == "__main__":
    from runpy import run_module
    run_module('fastfish.baseline.standing_wave', run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = import_module('fastfish.baseline.standing_wave')
