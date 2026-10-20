"""The built-in imagers pass the public imager contract suite (#23) with real CASA."""
import pytest

from snow.config import ImageGeometry
from snow.imaging.backends import Tclean
from snow.imaging.contract import CallableImager, ImagerContractError, ImagingRequest, ImagingResult
from snow.testing import ImagerContractSuite, simulate_point_source_ms


class TestTcleanContract(ImagerContractSuite):
    imager = Tclean(ImageGeometry(imsize=128, cell="1arcsec"), niter=200, nsigma=3.0)


def test_an_imager_that_writes_no_model_is_caught(tmp_path):
    """With real CASA: the model tclean wrote earlier is removed, so it can't hide the omission."""
    ms = simulate_point_source_ms(str(tmp_path / "target.ms"))
    tclean_imager = Tclean(ImageGeometry(imsize=64, cell="4arcsec"), niter=10)
    tclean_imager.image(ImagingRequest(ms, str(tmp_path / "first")))
    lazy = CallableImager(lambda request: ImagingResult("x", "x", "x"))

    with pytest.raises(ImagerContractError, match="MODEL_DATA does not exist"):
        lazy.image(ImagingRequest(ms, str(tmp_path / "second")))
