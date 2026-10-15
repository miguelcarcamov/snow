"""Imager statistics: an explicit stdv_pixels must work for CASA and FITS images (#29)."""
import pytest

from snow.imaging import imager as imager_module
from snow.imaging.imager import Imager


class StatisticsImager(Imager):

    def run(self, imagename=""):
        pass


@pytest.mark.parametrize(
    "method, calculator", [
        ("_calculate_statistics_msimage", "calculate_psnr_ms"),
        ("_calculate_statistics_fits", "calculate_psnr_fits")
    ]
)
@pytest.mark.parametrize("stdv_pixels, expected_pixels", [(None, 64), (32, 32)])
def test_statistics_use_explicit_or_default_noise_pixels(
    monkeypatch, method, calculator, stdv_pixels, expected_pixels
):
    received = {}

    def fake_calculator(signal, residual, pixels):
        received["pixels"] = pixels
        return 20.0, 2.0, 0.1

    monkeypatch.setattr(imager_module, calculator, fake_calculator)
    image = StatisticsImager(noise_pixels=64)

    getattr(image, method)("image", "residual", stdv_pixels=stdv_pixels)

    assert received["pixels"] == expected_pixels
    assert (image.psnr, image.peak, image.stdv) == (20.0, 2.0, 0.1)
