"""Compatibility alias; implementation: fastfish.processing.scheduler."""
from importlib import import_module
import sys

if __name__ == "__main__":
    from runpy import run_module
    run_module('fastfish.processing.scheduler', run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = import_module('fastfish.processing.scheduler')
