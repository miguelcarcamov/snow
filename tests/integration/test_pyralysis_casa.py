"""
The Pyralysis imager (#24) with real CASA and Pyralysis: the imager contract, and self-calibration
end to end. Skipped where Pyralysis isn't installed. Pyralysis is never imported here: it runs in
its own process, since it can't share one with casatools.
"""
import importlib.util
import warnings

import numpy as np
import pytest
from astropy.io import fits

from snow.config import ImageGeometry, SolveOptions
from snow.imaging import Pyralysis
from snow.selfcalibration import Phasecal, Selfcal
from snow.testing import ImagerContractSuite, simulate_point_source_ms
from snow.utils.selfcal_utils import open_table

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pyralysis") is None, reason="Pyralysis isn't installed"
)

GEOMETRY = ImageGeometry(imsize=128, cell="2arcsec")


def pyralysis_imager(**options) -> Pyralysis:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # no mask: the contract doesn't need one
        return Pyralysis(
            GEOMETRY, terms=["chi_squared"], optimizer="lbfgsb", max_iter=200, **options
        )


class TestPyralysisContract(ImagerContractSuite):
    imager = pyralysis_imager() if importlib.util.find_spec("pyralysis") else None


def test_options_are_checked_by_pyralysis():
    with pytest.raises(TypeError, match="Invalid Pyralysis options"):
        pyralysis_imager(linesearch="no-such-linesearch")
    with pytest.raises(TypeError, match="not a Pyralysis option"):
        pyralysis_imager(maxiter=10)


def test_settings_snow_controls_are_rejected():
    with pytest.raises(TypeError, match="ImageGeometry"):
        pyralysis_imager(imsize=256)


def test_self_calibration_with_pyralysis(tmp_path):
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"))
    # A box around the source, as for CLEAN: a model free in every pixel absorbs the gain errors
    mask = np.zeros((128, 128), dtype=np.int16)
    mask[60:69, 60:69] = 1
    fits.PrimaryHDU(mask).writeto(tmp_path / "mask.fits")
    imager = pyralysis_imager(mask_mode="file", mask_path=str(tmp_path / "mask.fits"))

    result = Selfcal(
        ms,
        imager, [Phasecal("inf"), Phasecal("10min")],
        solve=SolveOptions(refant="A0", minblperant=3)
    ).run()

    psnr = [record.psnr for record in result.rounds]
    assert all(record.accepted for record in result.rounds)
    assert psnr[0] < psnr[1] < psnr[2] and psnr[2] > 2 * psnr[0]
    # gaincal found the phase errors: the model didn't absorb them
    with open_table(result.stack.tables[0].path) as cal_table:
        phases = np.degrees(np.angle(cal_table.getcol("CPARAM")))
    assert np.max(np.abs(phases)) > 10.0
