"""Regression tests for #25: chaining self-calibrations must not silently discard solutions."""
from pathlib import Path

import pytest

from snow.imaging.imager import Imager
from snow.selfcalibration import AmpPhasecal, Ampcal, Phasecal


class FakeImager(Imager):

    def run(self, imagename=""):
        pass


@pytest.fixture
def ms(tmp_path) -> str:
    path = tmp_path / "target.ms"
    path.mkdir()
    (path / "table.dat").write_text("visibilities")
    return str(path)


@pytest.fixture
def imager(tmp_path) -> FakeImager:
    return FakeImager(output=str(tmp_path / "out"))


def make_caltable(tmp_path, name) -> str:
    path = tmp_path / name
    path.mkdir()
    return str(path)


def mark_as_run(selfcal, caltable) -> None:
    """Simulate a finished run() that accepted one calibration table."""
    selfcal._caltables.append(caltable)


@pytest.fixture
def phasecal(ms, imager, tmp_path) -> Phasecal:
    phs = Phasecal(visfile=ms, imager=imager, solint=["inf"])
    mark_as_run(phs, make_caltable(tmp_path, "pcal0"))
    return phs


def copies_of(ms) -> list:
    return sorted(p.name for p in Path(ms).parent.glob(Path(ms).stem + "_*"))


class TestVisfileFromPreviousSelfcal:

    def test_empty_visfile_continues_from_previous_measurement_set(self, phasecal, imager):
        apcal = AmpPhasecal(previous_selfcal=phasecal, imager=imager, solint=["inf"])

        assert Path(apcal.visfile).name == Path(phasecal.visfile).stem + "_ap0.ms"
        assert apcal.input_caltable == phasecal._caltables[-1]

    def test_explicit_previous_measurement_set_is_accepted(self, phasecal, imager):
        apcal = AmpPhasecal(
            visfile=phasecal.visfile, previous_selfcal=phasecal, imager=imager, solint=["inf"]
        )

        assert apcal.input_caltable == phasecal._caltables[-1]

    def test_original_measurement_set_raises_before_copying(self, ms, phasecal, imager):
        before = copies_of(ms)

        with pytest.raises(ValueError, match="not the measurement set of previous_selfcal"):
            AmpPhasecal(visfile=ms, previous_selfcal=phasecal, imager=imager, solint=["inf"])

        assert copies_of(ms) == before


class TestAmplitudeNeedsPhaseSolutions:

    @pytest.mark.parametrize("cls", [Ampcal, AmpPhasecal])
    def test_without_phase_solutions_raises_before_copying(self, cls, ms, imager):
        with pytest.raises(ValueError, match="needs phase solutions"):
            cls(visfile=ms, imager=imager, solint=["inf"])

        assert copies_of(ms) == []

    @pytest.mark.parametrize("cls", [Ampcal, AmpPhasecal])
    def test_previous_selfcal_that_has_not_run_raises(self, cls, ms, imager):
        phs = Phasecal(visfile=ms, imager=imager, solint=["inf"])

        with pytest.raises(ValueError, match="needs phase solutions"):
            cls(previous_selfcal=phs, imager=imager, solint=["inf"])

    @pytest.mark.parametrize("cls", [Ampcal, AmpPhasecal])
    def test_input_caltable_provides_phase_solutions(self, cls, ms, imager, tmp_path):
        caltable = make_caltable(tmp_path, "external_pcal")

        selfcal = cls(visfile=ms, imager=imager, solint=["inf"], input_caltable=caltable)

        assert selfcal.input_caltable == caltable

    def test_phase_solutions_found_through_the_chain(self, phasecal, imager, tmp_path):
        ampcal = Ampcal(previous_selfcal=phasecal, imager=imager, solint=["inf"])
        mark_as_run(ampcal, make_caltable(tmp_path, "acal0"))

        apcal = AmpPhasecal(previous_selfcal=ampcal, imager=imager, solint=["inf"])

        assert apcal.input_caltable == ampcal._caltables[-1]

    def test_non_incremental_ampphasecal_raises(self, phasecal, imager):
        with pytest.raises(ValueError, match="incremental=False"):
            AmpPhasecal(previous_selfcal=phasecal, imager=imager, solint=["inf"], incremental=False)


class TestPhasecalDoesNotDiscardSolutions:

    def test_after_amplitude_selfcal_raises(self, phasecal, imager, tmp_path):
        apcal = AmpPhasecal(previous_selfcal=phasecal, imager=imager, solint=["inf"])
        mark_as_run(apcal, make_caltable(tmp_path, "apcal0"))

        with pytest.raises(ValueError, match="would discard its amplitude solutions"):
            Phasecal(previous_selfcal=apcal, imager=imager, solint=["30s"])

    def test_after_phasecal_is_allowed(self, phasecal, imager):
        refined = Phasecal(previous_selfcal=phasecal, imager=imager, solint=["30s"])

        assert Path(refined.visfile).name == Path(phasecal.visfile).stem + "_p0.ms"

    def test_with_input_caltable_raises(self, ms, imager, tmp_path):
        with pytest.raises(ValueError, match="does not pre-apply input_caltable"):
            Phasecal(
                visfile=ms,
                imager=imager,
                solint=["inf"],
                input_caltable=make_caltable(tmp_path, "pcal_ext")
            )
