from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from snow.selfcalibration import ampcal, apcal, phasecal, selfcal

from .helpers import FakeImager

CASA_NAMES = (
    "gaincal", "applycal", "rmtables", "flagmanager", "delmod", "clearcal", "flagdata", "split",
    "statwt", "uvsub", "tb"
)


@pytest.fixture
def ms(tmp_path) -> str:
    path = tmp_path / "target.ms"
    path.mkdir()
    (path / "table.dat").write_text("visibilities")
    return str(path)


@pytest.fixture
def imager(tmp_path) -> FakeImager:
    return FakeImager(output=str(tmp_path / "out_"))


@pytest.fixture
def casa(monkeypatch) -> SimpleNamespace:
    """
    Replace every CASA task and tool used by the self-calibration modules with mocks, so whole
    run() loops execute without CASA. Set ``spw_ids``/``nspw`` to describe the measurement set.
    """
    mocks = SimpleNamespace(spw_ids=[0, 1, 2, 3], nspw=4)
    for name in CASA_NAMES:
        setattr(mocks, name, mock.MagicMock(name=name))
    # Like gaincal, create the output table so later steps can check that it exists
    mocks.gaincal.side_effect = lambda **kwargs: Path(kwargs["caltable"]).mkdir(exist_ok=True)
    for module in (selfcal, phasecal, ampcal, apcal):
        for name in CASA_NAMES:
            if hasattr(module, name):
                monkeypatch.setattr(module, name, getattr(mocks, name))
    monkeypatch.setattr(selfcal, "get_spw_ids", lambda ms_name, spw: mocks.spw_ids)
    monkeypatch.setattr(selfcal, "count_spectral_windows", lambda ms_name: mocks.nspw)
    return mocks
