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

from wenuray.config import ImageGeometry, SolveOptions
from wenuray.imaging import Pyralysis
from wenuray.imaging.contract import ImagingRequest
from wenuray.selfcalibration import Phasecal, Selfcal
from wenuray.testing import ImagerContractSuite, simulate_point_source_ms
from wenuray.utils.selfcal_utils import open_table

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pyralysis") is None, reason="Pyralysis isn't installed"
)

# Nyquist sampling: the longest baselines (~59 klambda) need cells of 1.74 arcsec or less. With
# coarser cells they fall outside the uv grid, which Pyralysis's fit doesn't represent
GEOMETRY = ImageGeometry(imsize=128, cell="1arcsec")


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


def test_settings_wenuray_controls_are_rejected():
    with pytest.raises(TypeError, match="ImageGeometry"):
        pyralysis_imager(imsize=256)


def box_mask(path) -> str:
    """A box around the source, as for CLEAN: a model free in every pixel absorbs the gain errors"""
    mask = np.zeros((128, 128), dtype=np.int16)
    mask[60:69, 60:69] = 1
    fits.PrimaryHDU(mask).writeto(path)
    return str(path)


def test_the_masked_model_holds_the_source_flux(tmp_path):
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"), phase_errors_in=("DATA", ))
    imager = pyralysis_imager(mask_mode="file", mask_path=box_mask(tmp_path / "mask.fits"))

    statistics = imager.image(ImagingRequest(ms, str(tmp_path / "masked"))).statistics()

    assert statistics.model_flux == pytest.approx(1.0, rel=0.05)


def test_a_measurement_set_read_in_several_chunks(tmp_path):
    """Dask runs the chunks in parallel threads, which Numba allows with OpenMP or TBB only"""
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"), phase_errors_in=("DATA", ))
    mask = box_mask(tmp_path / "mask.fits")

    whole, chunked = (
        pyralysis_imager(mask_mode="file", mask_path=mask,
                         **options).image(ImagingRequest(ms, str(tmp_path / name))).statistics()
        for name, options in (("whole", {}), ("chunked", {
            "row_chunks": 37
        }))
    )

    assert chunked.peak == pytest.approx(whole.peak, rel=0.01)
    assert chunked.model_flux == pytest.approx(whole.model_flux, rel=0.01)


def test_self_calibration_with_pyralysis(tmp_path):
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"))
    imager = pyralysis_imager(mask_mode="file", mask_path=box_mask(tmp_path / "mask.fits"))

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
