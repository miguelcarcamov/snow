"""Test doubles and helpers shared by the self-calibration tests."""
from dataclasses import dataclass, field
from typing import Callable, Optional

from snow.imaging.imager import Imager
from snow.selfcalibration.caltable import CalTable


@dataclass
class FakeImager(Imager):
    """
    Imager that reports a scripted sequence of PSNR values instead of imaging. ``write_model`` is
    called on every run, like a real imager saving its model to MODEL_DATA.
    """
    psnrs: list = field(default_factory=list)
    write_model: Optional[Callable[[], None]] = None

    def run(self, imagename=""):
        if self.write_model is not None:
            self.write_model()
        self.psnr = self.psnrs.pop(0) if self.psnrs else 1.0
        self.peak = self.psnr
        self.stdv = 1.0


def make_caltable(tmp_path, name) -> str:
    path = tmp_path / name
    path.mkdir()
    return str(path)


def mark_as_run(selfcal_object, caltable) -> None:
    """Simulate a finished run() that accepted one calibration table."""
    selfcal_object._caltables.append(caltable)
    selfcal_object._accept_caltable(CalTable(caltable, selfcal_object._calmode))
