"""Shared fixtures and constants for per-sec extract tests."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Add scripts/ to import path so tests can `import extract_per_sec_day`.
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "analysis"))

# Known local FITS files used as integration fixtures.
DATA_DIR = REPO_ROOT / "data"

HE_ENG_2017_BOXA = DATA_DIR / "1B/2017/20171001/0766/HXMT_1B_0766_20171001T000000_G002572_000_003.fits"
HE_EVT_20260410_HR07 = DATA_DIR / "1K/Y202604/20260410-3222/HXMT_20260410T07_HE-Evt_FFFFFF_V1_1K.FITS"
ORBIT_20260410_HR07  = DATA_DIR / "1K/Y202604/20260410-3222/HXMT_20260410T07_Orbit_FFFFFF_V1_1K.FITS"
ATT_20260410_HR07    = DATA_DIR / "1K/Y202604/20260410-3222/HXMT_20260410T07_Att_FFFFFF_V1_1K.FITS"


@pytest.fixture
def require_file():
    """Returns a function that skips if a path doesn't exist."""
    def _check(path: Path):
        if not path.exists():
            pytest.skip(f"Test fixture missing: {path}")
    return _check
