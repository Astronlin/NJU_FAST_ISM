"""Compatibility alias; implementation: fastfish.commands.shell."""
from importlib import import_module
import sys

if __name__ == "__main__":
    _module = import_module('fastfish.commands.shell')
    raise SystemExit(_module.main())
else:
    sys.modules[__name__] = import_module('fastfish.commands.shell')
