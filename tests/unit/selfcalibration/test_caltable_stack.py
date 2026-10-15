"""Regression tests for #26: every accepted table stays applied, each with its own spwmap."""
from pathlib import Path

import pytest

from snow.selfcalibration import AmpPhasecal, Ampcal, Phasecal

from .helpers import make_caltable


def names(paths) -> list:
    return [Path(path).name for path in paths]


def gaintables(task) -> list:
    return [names(call.kwargs["gaintable"]) for call in task.call_args_list]


def spwmaps(task) -> list:
    return [call.kwargs["spwmap"] for call in task.call_args_list]


@pytest.fixture
def phase_then_ampphase(casa, ms, imager):
    """Phase self-cal combining spws 1~3 of a 4-spw dataset, then amplitude-phase per spw."""
    casa.spw_ids = [1, 2, 3]
    imager.spw = "1~3"
    phscal = Phasecal(visfile=ms, imager=imager, solint=["inf", "60s"], combine="spw")
    phscal.run()
    casa.gaincal.reset_mock()
    casa.applycal.reset_mock()
    apcal = AmpPhasecal(previous_selfcal=phscal, imager=imager, solint=["inf", "30min"])
    apcal.run()
    return phscal, apcal


class TestPhasecal:

    def test_solves_fresh_and_replaces_previous_phase_table(self, casa, ms, imager):
        phscal = Phasecal(visfile=ms, imager=imager, solint=["inf", "60s"], combine="spw")
        phscal.run()

        assert gaintables(casa.gaincal) == [[], []]
        assert gaintables(casa.applycal) == [["out_pcal0"], ["out_pcal1"]]
        assert names(t.path for t in phscal._caltable_stack) == ["out_pcal1"]


class TestIncrementalStack:

    def test_ampphasecal_keeps_phase_table_in_every_iteration(self, casa, phase_then_ampphase):
        assert gaintables(casa.gaincal) == [["out_pcal1"], ["out_pcal1", "out_apcal_0"]]
        assert gaintables(casa.applycal) == [
            ["out_pcal1", "out_apcal_0"],
            ["out_pcal1", "out_apcal_0", "out_apcal_1"],
        ]

    def test_each_table_gets_its_own_spwmap(self, casa, phase_then_ampphase):
        phase_map = [1, 1, 1, 1]  # combined over spws 1~3 -> stored under spw 1
        assert spwmaps(casa.applycal) == [[phase_map, []], [phase_map, [], []]]

    def test_next_selfcal_preapplies_the_whole_stack(self, casa, imager, phase_then_ampphase):
        _, apcal = phase_then_ampphase
        casa.gaincal.reset_mock()

        Ampcal(previous_selfcal=apcal, imager=imager, solint=["1h"]).run()

        assert gaintables(casa.gaincal) == [["out_pcal1", "out_apcal_0", "out_apcal_1"]]


class TestRejectedIterations:

    def test_rejected_table_is_not_kept(self, casa, ms, imager):
        imager.psnrs = [10.0, 20.0, 15.0]  # initial, accepted, rejected
        phscal = Phasecal(
            visfile=ms, imager=imager, solint=["inf", "60s"], combine="spw", restore_psnr=True
        )
        phscal.run()
        casa.gaincal.reset_mock()

        AmpPhasecal(previous_selfcal=phscal, imager=imager, solint=["inf"]).run()

        assert names(t.path for t in phscal._caltable_stack) == ["out_pcal0"]
        assert gaintables(casa.gaincal) == [["out_pcal0"]]


class TestExplicitSpwmap:

    def test_explicit_spwmap_is_kept_but_warns_when_it_disagrees(self, casa, ms, imager):
        casa.spw_ids = [1, 2, 3]
        phscal = Phasecal(
            visfile=ms, imager=imager, solint=["inf"], combine="spw", spwmap=[0, 0, 0, 0]
        )

        with pytest.warns(UserWarning, match="lowest selected spectral window"):
            phscal.run()

        assert spwmaps(casa.applycal) == [[[0, 0, 0, 0]]]

    def test_shared_spwmap_is_not_applied_to_per_spw_solutions(self, casa, ms, imager, tmp_path):
        """The README passes spwmap=[0,0,0,0] meant for the combined phase table."""
        pcal = make_caltable(tmp_path, "pcal_ext")
        apcal = AmpPhasecal(
            visfile=ms, imager=imager, solint=["inf"], spwmap=[0, 0, 0, 0], input_caltable=pcal
        )

        apcal.run()

        assert gaintables(casa.applycal) == [["pcal_ext", "out_apcal_0"]]
        assert spwmaps(casa.applycal) == [[[0, 0, 0, 0], []]]
