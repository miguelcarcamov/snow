"""
The self-calibration loop (#22): the only one in snow. Steps provide data and hooks, the imager
follows the contract of snow.imaging.contract, the acceptance policies decide, and CalStack builds
every table list.

    result = Selfcal(
        ms="target.ms",
        imager=Tclean(ImageGeometry(imsize=1024, cell="0.05arcsec"), niter=1000, nsigma=3.0),
        steps=[Phasecal("inf"), Phasecal("60s"), AmpPhasecal("inf")],
        solve=SolveOptions(refant="DA51"),
    ).run()

The loop:

    validate the schedule                         (construction, before any I/O)
    CORRECTED_DATA = DATA, image round 0
    for each step:
        solve against DATA with the accepted stack pre-applied
        apply the stack plus the new table
        image CORRECTED_DATA (the imager writes MODEL_DATA)
        accept: push the table and snapshot the flags
        reject: restore flags and calibration, re-image the accepted model, stop
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence, Tuple

from casatasks import applycal, gaincal

from ..config import ApplyOptions, Selection, SolveOptions
from ..imaging.contract import Imager, ImagingRequest, ImagingResult
from ..utils.selfcal_utils import (
    columns_match, count_spectral_windows, get_spw_ids, is_column_in_ms, solution_summary
)
from .acceptance import AcceptancePolicy, RoundMetrics
from .caltable import CalStack, CalTable, derive_spwmap
from .result import RoundRecord, SelfcalResult
from .schedule import validate_schedule
from .steps import CalibrationStep
from .workspace import MSWorkspace

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Selfcal:
    """
    A self-calibration run. Construction validates everything and touches no file; run() does the
    work.

    Parameters
    ----------
    ms :
        Measurement set with the calibrated data in DATA. gaincal solves against DATA, so
        CORRECTED_DATA is reset to DATA when the run starts
    imager :
        The imager; it must follow the self-calibration contract
    steps :
        The calibration steps, in the order they run. Invalid orders raise ScheduleError here
    selection :
        One data selection for the solve, the apply and the imager
    solve :
        gaincal options shared by every step; each step may override them
    apply :
        applycal options
    acceptance :
        Policy for every step that doesn't set its own; None uses each step type's default
    output :
        Directory for the calibration tables, the images and the result record; by default
        <ms name>_selfcal next to the measurement set
    """
    ms: str
    imager: Imager
    steps: Tuple[CalibrationStep, ...]
    selection: Selection = Selection()
    solve: SolveOptions = SolveOptions()
    apply: ApplyOptions = ApplyOptions()
    acceptance: Optional[AcceptancePolicy] = None
    output: str = ""
    # Tables accepted before this run, e.g. by an earlier run (set by resume, #22c)
    initial_stack: CalStack = field(default=CalStack(), repr=False)

    def __post_init__(self):
        if not isinstance(self.imager, Imager):
            raise TypeError("imager must be an Imager, not {0!r}".format(self.imager))
        object.__setattr__(self, "steps", tuple(self.steps))
        validate_schedule(self.steps, self.initial_stack)
        for step in self.steps:
            # Unknown imaging overrides raise here, before anything runs
            self.imager.with_overrides(**step.imaging)
        if not self.output:
            object.__setattr__(self, "output", str(Path(self.ms).with_suffix("")) + "_selfcal")

    def run(self) -> SelfcalResult:
        """
        Runs the self-calibration and writes the result record to <output>/selfcal.toml
        """
        os.makedirs(self.output, exist_ok=True)
        workspace = MSWorkspace(self.ms, self.selection, self.apply)
        self._reset_corrected(workspace)
        stack = self.initial_stack

        result = self._image(self.imager, "round0")
        accepted = RoundMetrics(result.statistics())
        rounds = [self._record(0, None, "", result, accepted, True, ("initial image", ))]
        logger.info("Round 0 - PSNR: %.3f", accepted.statistics.psnr)
        workspace.snapshot("round0")
        snapshots = ["round0"]
        accepted_imager, accepted_name = self.imager, "round0"

        for index, step in enumerate(self.steps, start=1):
            name = "round{0}".format(index)
            imager = self.imager.with_overrides(**step.imaging)
            table = self._solve(step, stack, name)
            if table is None:
                # Nothing was applied or imaged, so there is nothing to restore
                rounds.append(self._failed_record(index, step, "gaincal found no solutions"))
                break
            self._apply(stack, table)
            result = self._image(imager, name)
            flagged, antennas = solution_summary(table.path)
            current = RoundMetrics(result.statistics(), step, flagged, antennas)
            decision = step.acceptance_policy(self.acceptance).decide(accepted, current)
            rounds.append(
                self._record(
                    index, step, table.path, result, current, decision.accepted, decision.reasons
                )
            )
            logger.info(
                "Round %d %r - PSNR: %.3f, %s: %s", index, step.solint, current.statistics.psnr,
                "accepted" if decision.accepted else "rejected", "; ".join(decision.reasons)
            )
            if not decision.accepted:
                workspace.restore(snapshots[-1], stack)
                # MODEL_DATA must hold the accepted model again
                self._image(accepted_imager, accepted_name + "_model")
                break
            stack = stack.push(table)
            workspace.snapshot(name)
            snapshots.append(name)
            accepted, accepted_imager, accepted_name = current, imager, name

        workspace.discard(*snapshots)
        selfcal_result = SelfcalResult(self.ms, stack, rounds)
        selfcal_result.to_toml(os.path.join(self.output, "selfcal.toml"))
        return selfcal_result

    def _reset_corrected(self, workspace: MSWorkspace) -> None:
        """
        gaincal solves against DATA, so CORRECTED_DATA must be DATA calibrated with the initial
        stack and nothing else
        """
        if len(self.initial_stack) == 0 and is_column_in_ms(self.ms, "CORRECTED_DATA") \
                and not columns_match(self.ms, "DATA", "CORRECTED_DATA"):
            logger.warning(
                "CORRECTED_DATA of %s differs from DATA and is reset. gaincal solves against "
                "DATA: split the calibrated data into DATA first if CORRECTED_DATA held them.",
                self.ms
            )
        workspace.apply_stack(self.initial_stack)

    def _image(self, imager: Imager, name: str) -> ImagingResult:
        return imager.image(
            ImagingRequest(self.ms, os.path.join(self.output, name), self.selection)
        )

    def _solve(self, step: CalibrationStep, stack: CalStack, name: str) -> Optional[CalTable]:
        """
        Solves the step's gains against DATA with the accepted stack pre-applied. Returns None if
        gaincal wrote no table (no solution reached minsnr)
        """
        path = os.path.join(self.output, "{0}_{1}.cal".format(name, step.calmode))
        gaincal(
            vis=self.ms,
            caltable=path,
            **step.gaincal_kwargs(self.solve, self.selection),
            **stack.gaincal_arguments()
        )
        if not os.path.exists(path):
            return None
        spwmap = derive_spwmap(
            step.solve_options(self.solve).combine, get_spw_ids(self.ms, self.selection.spw),
            count_spectral_windows(self.ms)
        )
        return CalTable(path, step.calmode, tuple(spwmap), self.apply.interp)

    def _apply(self, stack: CalStack, table: CalTable) -> None:
        applycal(
            vis=self.ms,
            field=self.selection.field,
            spw=self.selection.spw,
            uvrange=self.selection.uvrange,
            gainfield="",
            calwt=False,
            flagbackup=False,
            applymode=self.apply.applymode,
            **stack.applycal_arguments(table)
        )

    @staticmethod
    def _record(
        index: int, step: Optional[CalibrationStep], caltable: str, result: ImagingResult,
        metrics: RoundMetrics, accepted: bool, reasons: Sequence[str]
    ) -> RoundRecord:
        statistics = metrics.statistics
        return RoundRecord(
            index=index,
            step=_describe(step),
            caltable=caltable,
            image=result.restored,
            psnr=statistics.psnr,
            peak=statistics.peak,
            rms=statistics.rms,
            model_flux=statistics.model_flux,
            beam=statistics.beam,
            flagged_fraction=metrics.flagged_fraction,
            accepted=accepted,
            reasons=tuple(reasons),
        )

    @staticmethod
    def _failed_record(index: int, step: CalibrationStep, reason: str) -> RoundRecord:
        logger.info("Round %d %r - rejected: %s", index, step.solint, reason)
        nan = float("nan")
        return RoundRecord(
            index, _describe(step), "", "", nan, nan, nan, nan, None, 1.0, False, (reason, )
        )


def _describe(step: Optional[CalibrationStep]) -> str:
    if step is None:
        return ""
    details = [repr(step.solint)]
    details += ["solve={0}".format(dict(step.solve))] if step.solve else []
    details += ["imaging={0}".format(dict(step.imaging))] if step.imaging else []
    return "{0}({1})".format(type(step).__name__, ", ".join(details))
