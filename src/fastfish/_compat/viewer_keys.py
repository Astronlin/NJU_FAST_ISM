"""Compatibility alias; implementation: fastfish.viewer.keyboard."""
from importlib import import_module
import sys

if __name__ == "__main__":
    from runpy import run_module
    run_module('fastfish.viewer.keyboard', run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = import_module('fastfish.viewer.keyboard')
