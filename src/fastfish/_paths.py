"""Shared package locations, independent of individual module depth."""
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent.parent
DATA_ROOT = PACKAGE_ROOT / 'data'
