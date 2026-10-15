"""Shared test configuration.

CASA is not installable everywhere (wheels exist only for some platforms and Python versions), so
CASA-free tests stub ``casatasks``/``casatools`` when they are missing, and the ``casa`` fixture
replaces every CASA task and tool that snow uses with mocks. Tests that need real CASA must be
marked ``@pytest.mark.casa`` and live under ``tests/integration``.
"""
import importlib.util
import itertools
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

for _module in ("casatasks", "casatools"):
    if importlib.util.find_spec(_module) is None:
        sys.modules[_module] = mock.MagicMock(name=_module)

from snow.imaging import tclean  # noqa: E402
from snow.selfcalibration import ampcal, apcal, phasecal, selfcal  # noqa: E402
from snow import utils  # noqa: E402
from snow.utils import image_utils, selfcal_utils  # noqa: E402

CASA_NAMES = (
    "gaincal", "applycal", "rmtables", "flagmanager", "delmod", "clearcal", "flagdata", "split",
    "statwt", "uvsub", "tb", "tclean", "imstat"
)
CASA_MODULES = (selfcal, phasecal, ampcal, apcal, tclean, image_utils, selfcal_utils)


@pytest.fixture
def casa(monkeypatch) -> SimpleNamespace:
    """
    Replace every CASA task and tool used by snow with mocks, so whole self-calibration runs
    execute without CASA. Set ``spw_ids``/``nspw`` to describe the measurement set. Image
    statistics report a peak that grows on every call, so each iteration improves the PSNR.
    """
    mocks = SimpleNamespace(spw_ids=[0, 1, 2, 3], nspw=4)
    for name in CASA_NAMES:
        setattr(mocks, name, mock.MagicMock(name=name))
    # Like gaincal, create the output table so later steps can check that it exists
    mocks.gaincal.side_effect = lambda **kwargs: Path(kwargs["caltable"]).mkdir(exist_ok=True)
    # A calibrated measurement set without a model yet
    mocks.tb.colnames.return_value = ["DATA", "CORRECTED_DATA"]
    peaks = itertools.count(1.0)
    mocks.imstat.side_effect = lambda *args, **kwargs: {"max": [next(peaks)], "rms": [0.1]}
    for module in CASA_MODULES:
        for name in CASA_NAMES:
            if hasattr(module, name):
                monkeypatch.setattr(module, name, getattr(mocks, name))
    for module in (selfcal, utils):
        monkeypatch.setattr(module, "get_spw_ids", lambda ms_name, spw="": mocks.spw_ids)
        monkeypatch.setattr(module, "count_spectral_windows", lambda ms_name: mocks.nspw)
    return mocks
