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

import numpy as np
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
    "statwt", "uvsub", "tclean", "imstat"
)
CASA_MODULES = (selfcal, phasecal, ampcal, apcal, tclean, image_utils, selfcal_utils)


@pytest.fixture
def casa(monkeypatch) -> SimpleNamespace:
    """
    Replace every CASA task and tool used by snow with mocks, so whole self-calibration runs
    execute without CASA. Set ``spw_ids``/``nspw`` to describe the measurement set and
    ``columns`` for its data columns; ``write_model()`` adds MODEL_DATA like an imager does, and
    ``getcell`` returns the same values for every column. Image statistics report a peak that grows
    on every call, so each iteration improves the PSNR.
    """
    mocks = SimpleNamespace(spw_ids=[0, 1, 2, 3], nspw=4, tb=mock.MagicMock(name="table"))
    for name in CASA_NAMES:
        setattr(mocks, name, mock.MagicMock(name=name))
    # Like gaincal, create the output table so later steps can check that it exists
    mocks.gaincal.side_effect = lambda **kwargs: Path(kwargs["caltable"]).mkdir(exist_ok=True)
    # A calibrated measurement set without a model yet; tclean writes it, delmod deletes it
    mocks.columns = ["DATA", "CORRECTED_DATA"]
    mocks.tb.colnames.side_effect = lambda: list(mocks.columns)
    mocks.tb.nrows.return_value = 10
    mocks.tb.getcell.side_effect = lambda column, row: np.ones((2, 4), complex)

    def write_model(**kwargs):
        if "MODEL_DATA" not in mocks.columns:
            mocks.columns.append("MODEL_DATA")

    def delete_model(**kwargs):
        if "MODEL_DATA" in mocks.columns:
            mocks.columns.remove("MODEL_DATA")

    mocks.write_model = write_model
    mocks.tclean.side_effect = write_model
    mocks.delmod.side_effect = delete_model
    peaks = itertools.count(1.0)
    mocks.imstat.side_effect = lambda *args, **kwargs: {"max": [next(peaks)], "rms": [0.1]}
    for module in CASA_MODULES:
        for name in CASA_NAMES:
            if hasattr(module, name):
                monkeypatch.setattr(module, name, getattr(mocks, name))
    # open_table() creates its table tool through this factory
    monkeypatch.setattr(selfcal_utils, "table", lambda: mocks.tb)
    for module in (selfcal, utils):
        monkeypatch.setattr(module, "get_spw_ids", lambda ms_name, spw="": mocks.spw_ids)
        monkeypatch.setattr(module, "count_spectral_windows", lambda ms_name: mocks.nspw)
    return mocks
