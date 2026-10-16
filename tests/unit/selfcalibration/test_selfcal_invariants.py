"""Self-calibration invariants (#27): one selection, a saved model, calibrated DATA, clean flagging."""
import warnings

import numpy as np
import pytest

from snow.imaging import Tclean
from snow.selfcalibration import AmpPhasecal, Phasecal


def calls(task) -> list:
    """The calls of a task, without the flagdata summaries that report the flagged fraction."""
    return [call for call in task.call_args_list if call.kwargs.get("mode") != "summary"]


def selections(task) -> list:
    return [(call.kwargs["field"], call.kwargs["spw"]) for call in calls(task)]


class TestOneSelection:

    def test_imager_selection_is_used_to_solve_apply_and_flag(self, casa, ms, imager):
        imager.field, imager.spw = "target", "1~3"

        Phasecal(visfile=ms, imager=imager, solint=["inf"], flag_dataset=True).run()

        for task in (casa.gaincal, casa.applycal, casa.flagdata):
            assert selections(task) == [("target", "1~3")]

    def test_selfcal_selection_is_given_to_the_imager(self, ms, imager):
        phscal = Phasecal(visfile=ms, imager=imager, solint=["inf"], field="target", spw="0")

        assert (imager.field, imager.spw) == ("target", "0")
        assert (phscal.field, phscal.spw) == ("target", "0")

    @pytest.mark.parametrize("name", ["field", "spw"])
    def test_different_selections_raise(self, ms, imager, name):
        setattr(imager, name, "1")

        with pytest.raises(ValueError, match=name):
            Phasecal(visfile=ms, imager=imager, solint=["inf"], **{name: "2"})

    @pytest.mark.parametrize("varchange", ["varchange_selfcal", "varchange_imager"])
    def test_selection_changed_per_iteration_applies_everywhere(self, casa, ms, imager, varchange):
        Phasecal(
            visfile=ms, imager=imager, solint=["inf", "60s"], **{
                varchange: {
                    "spw": ["0", "1"]
                }
            }
        ).run()

        assert selections(casa.gaincal) == selections(casa.applycal) == [("", "0"), ("", "1")]
        assert imager.spw == "1"


class TestModel:

    def test_save_model_is_forced_on(self, ms, imager):
        imager.save_model = False

        with pytest.warns(UserWarning, match="save_model=True"):
            Phasecal(visfile=ms, imager=imager, solint=["inf"])

        assert imager.save_model

    def test_tclean_saves_the_model_column_when_forced_after_creation(self, casa, ms, tmp_path):
        clean = Tclean(output=str(tmp_path / "out"), cell="0.1arcsec", save_model=False)

        with pytest.warns(UserWarning, match="save_model"):
            Phasecal(visfile=ms, imager=clean, solint=["inf"]).run()

        assert {call.kwargs["savemodel"] for call in casa.tclean.call_args_list} == {"modelcolumn"}

    def test_fresh_run_removes_virtual_models_before_imaging(self, casa, ms, imager):
        casa.delmod.side_effect = lambda **kwargs: casa.events.append("delmod")
        casa.events = []
        imager.write_model = lambda: (casa.write_model(), casa.events.append("image"))

        Phasecal(visfile=ms, imager=imager, solint=["inf"]).run()

        assert casa.events[:2] == ["delmod", "image"]
        assert casa.delmod.call_args_list[0].kwargs == {
            "vis": imager.inputvis,
            "otf": True,
            "scr": False
        }

    def test_imager_that_writes_no_model_raises_before_solving(self, casa, ms, imager):
        imager.write_model = None

        with pytest.raises(RuntimeError, match="MODEL_DATA"):
            Phasecal(visfile=ms, imager=imager, solint=["inf"]).run()

        casa.gaincal.assert_not_called()


class TestCalibratedData:

    @pytest.fixture
    def corrected_differs(self, casa):
        casa.tb.getcell.side_effect = lambda column, row: np.full(
            (2, 4), 2.0 if column == "CORRECTED_DATA" else 1.0
        )

    def test_warns_when_corrected_data_holds_a_calibration(self, ms, imager, corrected_differs):
        with pytest.warns(UserWarning, match="split"):
            Phasecal(visfile=ms, imager=imager, solint=["inf"]).run()

    def test_no_warning_when_corrected_data_equals_data(self, ms, imager):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Phasecal(visfile=ms, imager=imager, solint=["inf"]).run()

    def test_continuing_self_calibration_is_not_checked(self, ms, imager):
        phscal = Phasecal(visfile=ms, imager=imager, solint=["inf"])
        phscal.run()

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            AmpPhasecal(previous_selfcal=phscal, imager=imager, solint=["inf"]).run()


class TestResidualFlagging:

    @pytest.mark.parametrize(
        "mode, used, unused", [
            ("rflag", "timedevscale", "timecutoff"),
            ("tfcrop", "timecutoff", "timedevscale"),
        ]
    )
    def test_only_parameters_of_the_mode_are_passed(self, casa, ms, imager, mode, used, unused):
        Phasecal(visfile=ms, imager=imager, solint=["inf"], flag_dataset=True, flag_mode=mode).run()

        kwargs = calls(casa.flagdata)[-1].kwargs
        assert kwargs["mode"] == mode
        assert used in kwargs and unused not in kwargs
        assert kwargs["flagbackup"] is False

    def test_mode_that_cannot_flag_residuals_raises_at_creation(self, casa, ms, imager):
        with pytest.raises(ValueError, match="flag_mode"):
            Phasecal(visfile=ms, imager=imager, solint=["inf"], flag_dataset=True, flag_mode="clip")


def test_flagged_fraction_is_reported_before_and_after_each_iteration(casa, ms, imager, capsys):
    phscal = Phasecal(visfile=ms, imager=imager, solint=["inf", "60s"], spw="0")
    phscal.run()

    assert phscal._flagged_history == pytest.approx([0.0, 0.1, 0.2])
    summaries = [c.kwargs for c in casa.flagdata.call_args_list if c.kwargs["mode"] == "summary"]
    assert {(kwargs["field"], kwargs["spw"]) for kwargs in summaries} == {("", "0")}
    assert "Flagged: 20.00% (+10.00% in this iteration)" in capsys.readouterr().out
