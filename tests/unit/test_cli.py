"""The snow command (#36): init, validate, show and run over TOML recipes."""
import os
import textwrap
import tomllib
from types import SimpleNamespace

import pytest

from snow import cli
from snow.config import SolveOptions
from snow.selfcalibration import Selfcal

RECIPE = """
    ms = "target.ms"

    [solve]
    refant = "DA51"

    [imager]
    type = "tclean"
    niter = 100
    geometry = { imsize = 256, cell = "0.1arcsec" }

    [[steps]]
    type = "phase"
    solint = "inf"

    [[steps]]
    type = "ampphase"
    solint = "inf"
"""


@pytest.fixture
def recipe(tmp_path):
    path = tmp_path / "selfcal.toml"
    path.write_text(textwrap.dedent(RECIPE))
    return str(path)


def snow(capsys, *argv):
    status = cli.main(list(argv))
    output = capsys.readouterr()
    return SimpleNamespace(status=status, out=output.out, err=output.err)


class TestInit:

    def test_writes_a_recipe_that_validates(self, capsys, tmp_path):
        result = snow(capsys, "init")
        (tmp_path / "selfcal.toml").write_text(result.out)

        assert result.status == 0
        assert snow(capsys, "validate", str(tmp_path / "selfcal.toml")).status == 0

    def test_to_a_file_that_it_does_not_overwrite(self, capsys, tmp_path):
        path = str(tmp_path / "selfcal.toml")

        assert snow(capsys, "init", "-o", path).status == 0
        result = snow(capsys, "init", "-o", path)
        assert result.status == 1
        assert "already exists" in result.err

    def test_ships_no_instrument_values(self, capsys):
        """snow's defaults come from theory or the data (#36)"""
        assert "preset = " not in snow(capsys, "init").out.replace("# preset = ", "")


class TestValidate:

    def test_lists_the_run(self, capsys, recipe, tmp_path):
        result = snow(capsys, "validate", recipe)

        assert result.status == 0
        assert "valid" in result.out
        assert str(tmp_path / "target.ms") in result.out
        assert "step 2:  AmpPhasecal('inf')" in result.out
        assert "does not exist" in result.out  # a note, not an error

    def test_a_typo_is_the_constructor_error(self, capsys, recipe):
        with open(recipe) as file:
            text = file.read()
        with open(recipe, "w") as file:
            file.write(text.replace("refant", "refnat"))

        result = snow(capsys, "validate", recipe)

        assert result.status == 1
        assert result.err.startswith("snow: error:") and "refnat" in result.err

    def test_an_invalid_schedule(self, capsys, recipe):
        with open(recipe) as file:
            text = file.read()
        with open(recipe, "w") as file:
            file.write(text.replace('type = "phase"', 'type = "amp"'))

        result = snow(capsys, "validate", recipe)

        assert result.status == 1
        assert "phase" in result.err

    def test_traceback_reraises(self, recipe):
        with pytest.raises(FileNotFoundError):
            cli.main(["--traceback", "validate", recipe + ".missing"])


class TestShow:

    def test_prints_the_recipe_with_its_preset(self, capsys, recipe, tmp_path):
        (tmp_path / "group.toml").write_text('[solve]\nminsnr = 5.0\n')
        with open(recipe) as file:
            text = file.read()
        with open(recipe, "w") as file:
            file.write('preset = "group.toml"\n' + text)

        shown = tomllib.loads(snow(capsys, "show", recipe).out)

        assert shown["solve"] == {"refant": "DA51", "minsnr": 5.0}
        assert shown["ms"] == "target.ms"

    def test_defaults_prints_every_value(self, capsys, recipe):
        shown = tomllib.loads(snow(capsys, "show", "--defaults", recipe).out)

        assert shown["solve"]["minblperant"] == 4
        assert shown["output"] == "target_selfcal"
        assert shown["imager"]["weighting"] == {
            "weighting": "natural",
            "robust": 0.5,
            "uvtaper": []
        }
        assert shown["steps"][0]["flagging"] == {"type": "none"}

    def test_the_resolved_recipe_builds_the_same_run(self, capsys, recipe, tmp_path):
        resolved = tmp_path / "resolved.toml"
        resolved.write_text(snow(capsys, "show", "--defaults", recipe).out)

        assert Selfcal.from_toml(resolved) == Selfcal.from_toml(recipe)

    def test_ms_and_output_overrides(self, capsys, recipe, tmp_path):
        other = str(tmp_path / "data" / "other.ms")

        shown = tomllib.loads(snow(capsys, "show", "--defaults", recipe, "--ms", other).out)

        # Without --output, the output follows the new measurement set
        assert shown["ms"] == os.path.join("data", "other.ms")
        assert shown["output"] == os.path.join("data", "other_selfcal")


def test_run_runs_the_recipe(capsys, recipe, monkeypatch):
    runs = []

    def run(self, observers):
        runs.append((self, observers))
        record = SimpleNamespace(image="round2.image")
        return SimpleNamespace(rounds=[record] * 3, accepted_rounds=[record] * 2, final=record)

    monkeypatch.setattr(Selfcal, "run", run)

    result = snow(capsys, "run", "--quiet", "--plots", recipe)

    assert result.status == 0
    assert "2 of 3 rounds accepted; final image round2.image" in result.out
    (selfcal, observers), = runs
    assert selfcal.solve == SolveOptions(refant="DA51")
    assert [type(observer).__name__ for observer in observers] == \
        ["LoggingObserver", "GainPlotObserver"]
