"""The README and main_files examples must run against the current API (#31)."""
import re
import runpy
import sys
from pathlib import Path
from unittest import mock

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "main_files"


def readme_example() -> str:
    readme = (ROOT / "README.md").read_text()
    section = readme[readme.index("### Basic Usage Example"):]
    return re.search(r"```python\n(.*?)```", section, re.DOTALL).group(1)


def run_example(script: Path, *args: str) -> None:
    with mock.patch.object(sys, "argv", [script.name, *args]):
        runpy.run_path(str(script), run_name="__main__")


@pytest.fixture
def workdir(tmp_path, monkeypatch) -> Path:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "target.ms").mkdir()
    return tmp_path


@pytest.fixture
def casatasks_module(monkeypatch):
    """casatasks as imported by the example scripts themselves (e.g. flagdata, mstransform)."""
    module = mock.MagicMock(name="casatasks")
    module.mstransform.side_effect = lambda **kwargs: Path(kwargs["outputvis"]).mkdir()
    monkeypatch.setitem(sys.modules, "casatasks", module)
    return module


def test_readme_example_is_main_alma():
    assert readme_example() == (EXAMPLES / "main_ALMA.py").read_text()


@pytest.mark.parametrize(
    "script, last_phase_table",
    [("main_ALMA.py", "target_scpcal1"), ("main_VLA.py", "target_scpcal0")]
)
def test_phase_then_ampphase_examples_run(script, last_phase_table, casa, workdir):
    run_example(EXAMPLES / script, "target.ms", "target_sc")

    # The output comes from the amplitude-phase copy of the phase self-calibrated dataset, and the
    # phase solutions are still applied at the end
    assert Path(casa.split.call_args.kwargs["vis"]).name == "target_p0_ap0.ms"
    assert [Path(t).name for t in casa.applycal.call_args.kwargs["gaintable"]
            ] == [last_phase_table, "target_scapcal_0"]


@pytest.mark.parametrize(
    "flag_args, flagged_spws", [
        ([], []), (["--flag"], [""]),
        (["--flag", "--rflag-spw", "0", "--rflag-spw", "1~7"], ["0", "1~7"])
    ]
)
def test_emerlin_example_flags_before_selfcal(
    flag_args, flagged_spws, casa, casatasks_module, workdir
):
    casa.nspw, casa.spw_ids = 8, list(range(8))
    (workdir / "pcal").mkdir()

    run_example(EXAMPLES / "main_eMERLIN.py", "target.ms", "target_sc", "1252+5634", *flag_args)

    rflag_calls = [
        call for call in casatasks_module.flagdata.call_args_list if call.kwargs["mode"] == "rflag"
    ]
    assert [call.kwargs["spw"] for call in rflag_calls] == flagged_spws


@pytest.mark.parametrize(
    "phase_args, phase_spwmap", [([], [0] * 8), (["--phase-caltable-per-spw"], [])]
)
def test_emerlin_example_preapplies_phase_table(
    phase_args, phase_spwmap, casa, casatasks_module, workdir
):
    casa.nspw, casa.spw_ids = 8, list(range(8))
    (workdir / "pcal").mkdir()

    run_example(EXAMPLES / "main_eMERLIN.py", "target.ms", "target_sc", "1252+5634", *phase_args)

    applied = casa.applycal.call_args.kwargs
    assert applied["gaintable"][0] == "pcal"
    assert applied["spwmap"][0] == phase_spwmap
