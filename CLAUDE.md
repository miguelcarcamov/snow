# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this project is

snow is a Python framework for radio-interferometric **self-calibration** on top of CASA (`casatasks`/`casatools`), with pluggable imagers (`tclean`, GPUvmem, and Pyralysis as the target). Its purpose is to let astronomers describe a self-cal procedure **briefly and reproducibly** instead of rewriting a one-off script for every dataset. Every change should make that easier, and stay correct by self-cal theory.

Source lives in `src/snow/`: `selfcalibration/` (self-cal loop), `imaging/` (imager backends) and `utils/`. Example scripts are in `main_files/`.

## Canonical conventions — read these first

The conventions live in **`.cursor/rules/*.mdc`** (see `.cursor/rules/README.md`) and are mirrored in **`CONVENTIONS.md`**. This file doesn't repeat them; it only adds operational notes. If the two disagree, the `.mdc` files win; regenerate the mirror with `python tools/sync_conventions.py`.

Highlights to keep in mind:

- **Architecture** (`design.mdc`):
  - One self-cal loop, in `Selfcal`.
  - Calibration types are frozen-dataclass subclasses of `CalibrationStep` (`Phasecal`, `Ampcal`, `AmpPhasecal`), in any valid order.
  - `CalStack` is the only place that builds `gaincal`/`applycal` table lists.
  - Theory rules are `multimethod` overloads of `check_transition`; an invalid schedule raises `ScheduleError` before any I/O.
- **Parameter budget:** at most ~8 constructor parameters, all of them used. Group parameters into frozen config objects (`Selection`, `SolveOptions`, `ApplyOptions`, `ImageGeometry`, `Weighting`). Prefer presets and data-derived defaults over new parameters.
- **Pluggable imagers:** subclass `Imager` and implement `_reconstruct` (or wrap a function with `CallableImager`). The final `image()` method enforces the contract: image `CORRECTED_DATA`, write `MODEL_DATA`.
- **Self-cal invariants** (`selfcal-domain.mdc`):
  - `gaincal` solves against `DATA`.
  - `MODEL_DATA` must hold the accepted model.
  - Pre-apply the accepted table stack; `spwmap`/`interp` belong to each table.
  - One selection for solve, apply and imaging.
  - Never use PSNR alone to accept a round.
- No I/O in constructors, no module-level CASA tool singletons, no `os.system`, no files with fixed names in the working directory (`casa-io.mdc`).
- Imagers pass backend options through with the backend's own parameter names; Pyralysis and GPUvmem are optional extras imported lazily.
- Legacy code predates these rules. Follow them in new or changed code; don't refactor beyond the scope of the issue you're working on.

## Roadmap and issues

The redesign is tracked on GitHub (`miguelcarcamov/snow`), milestone **"v1.0 redesign"**:

- Epic: #22.
- Imager interface: #23.
- Pyralysis backend: #24 (the end goal).
- Science bugs: #25, #26, #27.
- Stopping rule: #28.
- Smaller bugs: #29.
- Packaging and Python version: #30.
- Docs and examples: #31.
- Tests: #32.
- CI: #33, #34.

Check the relevant issue before starting work, and reference it in commits and PRs (`Closes #NN`).

Known blocker: snow pins Python `<3.11`, while Pyralysis needs `≥3.11` (#30).

## Environment

```bash
pip install -e .            # pulls casatools/casatasks wheels (Linux/macOS only)
pre-commit install
```

- CASA wheels exist only for specific Python versions; check `requires-python` in `pyproject.toml`.
- CASA needs measures data (`casaconfig`) the first time `casatools` is imported. Without network access, provide a measures directory.
- `import snow.imaging` currently imports Pyralysis at module level (#29/#30); install Pyralysis or import `snow.imaging.tclean` directly.
- The full Docker image (`ghcr.io/miguelcarcamov/snow:base`) has CUDA, casacore and GPUvmem. It's only needed for GPUvmem.

## Tests and linting

```bash
pytest tests/unit/                     # CASA-free unit tests (suite being built, #32)
pytest -m casa tests/integration/      # needs CASA
pre-commit run --all-files             # yapf, YAML/TOML formatting, markdownlint, workflow schema
python tools/sync_conventions.py --check
```

No suite exists yet (#32). When you add behaviour, add the test under `tests/unit/` mirroring `src/snow/`.

## Commits and branches

- `development` is the integration branch; `master` is the release branch.
- Work on `feature/<topic>` or `fix/<topic>` branches with PRs into `development`.
- Commit prefixes: `CHG:` and `FIX:`.
- Don't push to `development` or `master` directly unless asked.
