"""
Runs one Pyralysis reconstruction for snow's Pyralysis imager (#24), in a process of its own.

casatools and python-casacore (which Pyralysis reads measurement sets with) each bundle a casacore
build, and the two crash in one process. snow therefore runs this file as a script, so it imports
neither snow nor casatools. It reads a JSON request from stdin and prints a JSON reply as the last
line of its output (Pyralysis prints its progress before it):

    {"check": true, "options": {...}}                    -> {"error": ""} (validates the options)
    {"ms": ..., "taql": ..., "options": {...}, "outputs": {"model": ..., "residual": ..., "restored": ...}}
                                                         -> {"error": ""} after writing MODEL_DATA
"""
import argparse
import contextlib
import io
import json
import os
import sys

# snow's pyralysis.py lives next to this file and would shadow the pyralysis package
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:] = [path for path in sys.path if os.path.abspath(path or ".") != _HERE]
from dataclasses import dataclass
from pathlib import Path


def parser() -> argparse.ArgumentParser:
    from pyralysis.cli.optimization import add_optimization_cli_args

    result = argparse.ArgumentParser(prog="pyralysis", exit_on_error=False)
    add_optimization_cli_args(result)
    return result


def to_argv(cli: argparse.ArgumentParser, options: dict) -> list:
    """
    Pyralysis CLI arguments from options named like the CLI's destinations (max_iter=100, ...)
    """
    actions = {action.dest: action for action in cli._actions if action.option_strings}
    argv = []
    for name, value in options.items():
        if name not in actions:
            raise ValueError("{0!r} is not a Pyralysis option".format(name))
        action = actions[name]
        flags = action.option_strings
        if isinstance(action, argparse.BooleanOptionalAction):
            argv.append(flags[0] if value else "--no-" + flags[0][2:])
        elif isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            if bool(value) == isinstance(action, argparse._StoreTrueAction):
                argv.append(flags[0])
        elif isinstance(value, (list, tuple)):
            if action.nargs in (None, "?", 1):  # a repeated option, e.g. --term a --term b
                for item in value:
                    argv += [flags[0], str(item)]
            else:
                argv += [flags[0]] + [str(item) for item in value]
        else:
            argv += [flags[0], str(value)]
    return argv


def parse(options: dict, extra: list) -> argparse.Namespace:
    """
    Parses the options with Pyralysis's own CLI parser, so they are checked by Pyralysis
    """
    cli = parser()
    messages = io.StringIO()
    with contextlib.redirect_stderr(messages):
        try:
            return cli.parse_args(to_argv(cli, options) + extra)
        except argparse.ArgumentError as error:
            raise ValueError(str(error)) from None
        except SystemExit:  # argparse reports some errors by exiting
            lines = messages.getvalue().strip().splitlines()
            raise ValueError(lines[-1] if lines else "invalid Pyralysis options") from None


def run(request: dict) -> None:
    from pyralysis.cli.optimization import (
        build_optimizer, build_pipeline_imager_config, checkpoint_final_model_residuals,
        normalize_dask_runtime_for_gpu, optimize_kwargs_from_args, parse_dask_runtime_options,
        parse_terms
    )
    from pyralysis.io import DaskMS
    from pyralysis.optimization import ObjectiveFunction
    from pyralysis.pipelines.base import Step
    from pyralysis.pipelines.imager import (
        ImagerContext, ImagerPipeline, SetupDaskCluster, TeardownDaskCluster
    )
    from pyralysis.pipelines.imager.steps import (
        ApplyVisibilityWeighting, ExportOptimizedModel, ExportResults, RestoreImage, SetupImageGrid,
        SetupMeasurementOperator
    )

    ms, outputs = request["ms"], request["outputs"]
    args = parse(
        request["options"], [
            "--input-path", ms, "--output-model", outputs["model"], "--output-residual",
            outputs["residual"], "--output-restored", outputs["restored"], "--restore",
            "--no-hermitian-symmetry"
        ]
    )
    for name in ("input_path", "output_model", "output_residual", "output_restored"):
        setattr(args, name, Path(getattr(args, name)))
    dask_options = normalize_dask_runtime_for_gpu(parse_dask_runtime_options(args))
    config = build_pipeline_imager_config(args, dask_options.to_pipeline_config())

    # LoadData can't load CORRECTED_DATA or select data yet (pyralysis#358)
    @dataclass
    class LoadCorrectedData(Step[ImagerContext]):

        def execute(self, context):
            context.dataset = DaskMS(input_name=ms).read(
                columns=["CORRECTED_DATA"],
                taql_query=request.get("taql") or None,
                calculate_psf=True
            )
            return context

    @dataclass
    class Reconstruct(Step[ImagerContext]):

        def execute(self, context):
            terms = parse_terms(args.term, context.measurement_operator, dataset=context.dataset)
            context.objective = ObjectiveFunction(
                terms=terms, parameter=context.image, persist_gradient=True
            )
            optimizer = build_optimizer(
                args.optimizer,
                context.objective,
                context.image,
                args,
                mask=context.mask,
                dataset=context.dataset
            )
            context.image = optimizer.optimize(**optimize_kwargs_from_args(args, optimizer))
            context.optimizer = optimizer
            return context

    # The public API for this is pyralysis#359
    @dataclass
    class WriteModelData(Step[ImagerContext]):

        def execute(self, context):
            checkpoint_final_model_residuals(
                dataset=context.dataset,
                result=context.image,
                model_vis=context.measurement_operator,
                input_path=Path(ms),
                transient_resources=(
                    context.objective, context.optimizer, context.measurement_operator
                ),
            )
            context.config["residuals_ready"] = True
            return context

    steps = [
        SetupDaskCluster(),
        LoadCorrectedData(),
        ApplyVisibilityWeighting(),
        SetupImageGrid(),
        SetupMeasurementOperator(),
        Reconstruct(),
        ExportOptimizedModel(),
        WriteModelData(),
        RestoreImage(),
        ExportResults(),
        TeardownDaskCluster(),
    ]
    context = ImagerContext(config=config)
    try:
        ImagerPipeline(steps=steps).run(context)
    finally:
        TeardownDaskCluster().execute(context)


def main() -> int:
    request = json.load(sys.stdin)
    try:
        if request.get("check"):
            parse(request["options"], ["--input-path", "check.ms"])
        else:
            run(request)
    except Exception as error:  # reported to snow, which raises it
        print(json.dumps({"error": "{0}: {1}".format(type(error).__name__, error)}))
        return 1
    print(json.dumps({"error": ""}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
