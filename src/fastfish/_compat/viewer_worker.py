"""Compatibility alias; implementation: fastfish.viewer.worker."""
from importlib import import_module
import sys

if __name__ == "__main__":
    _module = import_module('fastfish.viewer.worker')
    raise SystemExit(_module.main())
else:
    sys.modules[__name__] = import_module('fastfish.viewer.worker')
