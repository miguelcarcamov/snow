"""
End-to-end self-calibration with real CASA: simulate a small VLA observation of a point source
with seeded phase errors and noise, then run phase self-calibration imaged with tclean.
"""
import os

import pytest
from casatools import table

from snow.imaging import Tclean
from snow.selfcalibration import Phasecal
from snow.testing import simulate_point_source_ms


def columns(path: str) -> list:
    tb = table()
    tb.open(path)
    try:
        return tb.colnames()
    finally:
        tb.close()


@pytest.fixture(scope="module")
def phasecal(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("selfcal")
    # CASA writes scratch files to the working directory
    previous_cwd = os.getcwd()
    os.chdir(workdir)
    ms = simulate_point_source_ms(str(workdir / "target.ms"))
    imager = Tclean(
        output=str(workdir / "target"),
        cell="2arcsec",
        M=128,
        N=128,
        niter=200,
        specmode="mfs",
        use_mask="user",  # no mask file: clean the whole image
        nsigma=3.0
    )
    phscal = Phasecal(
        visfile=ms,
        imager=imager,
        solint=["inf", "10min"],
        refant="A0",
        minblperant=3,
        want_plot=False
    )
    try:
        phscal.run()
    finally:
        os.chdir(previous_cwd)
    return phscal


def test_every_iteration_solves_and_images(phasecal):
    assert len(phasecal._caltables) == 2
    assert len(phasecal._psnr_history) == 3


def test_calibrated_data_and_model_are_written(phasecal):
    assert {"MODEL_DATA", "CORRECTED_DATA"} <= set(columns(phasecal.visfile))


def test_self_calibration_improves_the_image(phasecal):
    # The corruption is seeded, so this is deterministic: about 33 -> 55 -> 99
    psnr = phasecal._psnr_history
    assert psnr[0] < psnr[1] < psnr[2]
    assert psnr[2] > 2 * psnr[0]
