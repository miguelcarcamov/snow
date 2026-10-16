"""The new self-calibration loop (#22c) end to end with real CASA."""
import os

import numpy as np
import pytest
from casatasks import flagmanager

from snow.config import ImageGeometry, SolveOptions
from snow.imaging.backends import Tclean
from snow.selfcalibration.acceptance import AllOf, PSNRImproves
from snow.selfcalibration.loop import Selfcal
from snow.selfcalibration.result import SelfcalResult
from snow.selfcalibration.steps import Phasecal
from snow.testing import simulate_point_source_ms
from snow.utils.selfcal_utils import model_column_problem, open_table

IMAGER = Tclean(ImageGeometry(imsize=128, cell="2arcsec"), niter=200, nsigma=3.0)
SOLVE = SolveOptions(refant="A0", minblperant=3)


def snow_flag_versions(ms: str) -> list:
    versions = flagmanager(vis=ms, mode="list")
    return [v["name"] for k, v in versions.items() if k != "MS" and v["name"].startswith("snow_")]


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    ms = simulate_point_source_ms(str(tmp_path_factory.mktemp("loop") / "target.ms"))
    result = Selfcal(ms, IMAGER, [Phasecal("inf"), Phasecal("10min")], solve=SOLVE).run()
    return result


def test_every_round_improves_and_is_kept(accepted):
    psnr = [record.psnr for record in accepted.rounds]

    assert all(record.accepted for record in accepted.rounds)
    assert psnr[0] < psnr[1] < psnr[2]
    assert psnr[2] > 2 * psnr[0]
    assert len(accepted.stack) == 2


def test_the_measurement_set_holds_the_accepted_state(accepted):
    assert model_column_problem(accepted.ms) == ""
    with open_table(accepted.ms) as ms_table:
        assert not np.allclose(ms_table.getcell("DATA", 0), ms_table.getcell("CORRECTED_DATA", 0))
    assert snow_flag_versions(accepted.ms) == []


def test_the_result_record_round_trips(accepted):
    path = os.path.join(os.path.dirname(accepted.final.image), "selfcal.toml")

    assert SelfcalResult.from_toml(path) == accepted


def test_a_rejected_round_is_rolled_back(tmp_path):
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"))
    impossible = AllOf(PSNRImproves(min_ratio=100.0))

    result = Selfcal(ms, IMAGER, [Phasecal("inf", acceptance=impossible)], solve=SOLVE).run()

    assert [record.accepted for record in result.rounds] == [True, False]
    assert len(result.stack) == 0
    with open_table(ms) as ms_table:
        # No calibration accepted: CORRECTED_DATA is DATA again, and the round-0 model is back
        assert np.allclose(ms_table.getcell("DATA", 0), ms_table.getcell("CORRECTED_DATA", 0))
    assert model_column_problem(ms) == ""
    assert snow_flag_versions(ms) == []
