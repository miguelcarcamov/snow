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


@pytest.mark.parametrize("flag", [False, True])
def test_emerlin_example_runs(flag, casa, casatasks_module, workdir):
    casa.nspw, casa.spw_ids = 8, list(range(8))
    (workdir / "pcal1").mkdir()
    args = ["target.ms", "target_sc", "1252+5634"] + (["--flag"] if flag else [])

    run_example(EXAMPLES / "main_eMERLIN.py", *args)

    assert casatasks_module.flagdata.called == flag
    assert casa.applycal.call_args.kwargs["gaintable"][0] == "pcal1"
